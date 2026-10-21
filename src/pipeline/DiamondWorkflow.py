from __future__ import annotations

import hashlib
import json
import logging
import time
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Collection, Iterable

import pandas as pd

from pipeline.FastaUtils import prepare_query_fasta
from pipeline.RunSteps import ProgressCallback, report_progress
from pipeline.ToolEnv import tool_cmd
from pipeline.Utils import CATEGORIES


def dmnd_stamp_path(dmnd: Path) -> Path:
    return Path(dmnd).parent / f"{Path(dmnd).name}.source.json"


def _fasta_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_dmnd_stamp(fasta: Path, dmnd: Path) -> None:
    """Record which FASTA content a .dmnd was built from."""
    fasta = Path(fasta)
    stat = fasta.stat()
    stamp = {
        "fasta": str(fasta),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": _fasta_sha256(fasta),
    }
    dmnd_stamp_path(dmnd).write_text(json.dumps(stamp, indent=2) + "\n", encoding="utf-8")


def dmnd_is_current(fasta: Path, dmnd: Path) -> bool:
    """True when ``dmnd`` exists and was built from the current content of ``fasta``.

    With a stamp (written by build_diamond_db_for_files) the FASTA content is compared; size
    and mtime unchanged short-circuit the hash. Without one (older runs) the .dmnd must be
    at least as new as the FASTA, so a marker database rebuilt in DB/seq/markers after the
    run makes the run's .dmnd stale.
    """
    fasta, dmnd = Path(fasta), Path(dmnd)
    try:
        if not dmnd.is_file() or dmnd.stat().st_size == 0 or not fasta.is_file():
            return False
        stat = fasta.stat()
        stamp_path = dmnd_stamp_path(dmnd)
        if stamp_path.is_file():
            stamp = json.loads(stamp_path.read_text(encoding="utf-8"))
            if stamp.get("size") == stat.st_size and stamp.get("mtime_ns") == stat.st_mtime_ns:
                return True
            return stamp.get("size") == stat.st_size and stamp.get("sha256") == _fasta_sha256(fasta)
        return dmnd.stat().st_mtime_ns >= stat.st_mtime_ns
    except (OSError, ValueError):
        return False


def hit_tsv_is_current(tsv: Path, db: Path) -> bool:
    """A DIAMOND hit TSV is reusable only when it is not older than the database it searched."""
    try:
        return Path(tsv).is_file() and Path(tsv).stat().st_mtime_ns >= Path(db).stat().st_mtime_ns
    except OSError:
        return False


def run_diamond_annotations(
    genome_files: Iterable[Path],
    db_files: Iterable[Path],
    output_path: Path,
    logger: logging.Logger,
    threads: int = 1,
    max_target_seqs: int = 1,
    sanitize: bool = True,
    force: bool = False,
    progress: ProgressCallback = None,
):
    """Search each genome/query FASTA against each DB, one TSV per pair.

    progress(done, total) is called after each genome/query FASTA (all its DBs).

    An existing TSV is kept when it is not older than its DB, unless force is set:
    the caller then knows the query changed (for example new AnnotationIDs) and the
    old hits must not be reused.
    """
    annotation_dir = output_path
    annotation_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=output_path) as tdir:
        tdir = Path(tdir)

        genome_files = list(genome_files)
        db_files = list(db_files)
        for done, genome in enumerate(genome_files):
            report_progress(progress, done, len(genome_files))
            if not genome.exists() or not genome.is_file():
                logger.warning(f"Skipping missing genome file: {genome}")
                continue

            query_fa = prepare_query_fasta(genome, tdir, logger, write_map=False) if sanitize else genome

            for db in db_files:
                if not db.exists() or not db.is_file():
                    logger.warning(f"Skipping missing DB file: {db}")
                    continue

                out_file = annotation_dir / f"{genome.stem}-vs-{db.stem}.tsv"
                if not force and hit_tsv_is_current(out_file, db):
                    logger.info(f"Skipping existing annotation: {out_file.name}")
                    continue
                if force and out_file.exists():
                    logger.info(f"Recomputing {out_file.name}: the query changed since the hits were written.")
                elif out_file.exists():
                    logger.info(f"Recomputing {out_file.name}: {db.name} is newer than the hits.")

                with tempfile.TemporaryDirectory(dir=annotation_dir) as tmpdir:
                    tmp_out = Path(tmpdir) / out_file.name
                    cmd = [
                        *tool_cmd("diamond"), "blastx",
                        "--query", str(query_fa),
                        "--db", str(db),
                        "--out", str(tmp_out),
                        "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                        "qlen", "slen",
                        "--threads", str(threads),
                    "--id", "80",
                    "--subject-cover", "80",
                    "--max-target-seqs", "0",
                    "--max-hsps", "7",
                    "--evalue", "1e-5",
                    "--masking", "0",
                    "--quiet",
                    "--fast",
                    "--block-size", "16",
                    "--index-chunks", "1",
                ]

                    logger.info(f"Running DIAMOND blastx: {genome.name} vs {db.name} with command: {cmd}")
                    try:
                        subprocess.run(cmd, check=True)
                        if not tmp_out.exists():
                            # DIAMOND can succeed without writing any rows when there are no hits.
                            tmp_out.touch()
                            logger.info(
                                "DIAMOND produced no hits for %s vs %s; writing empty TSV.",
                                genome.name,
                                db.name,
                            )
                        shutil.move(str(tmp_out), out_file)
                        logger.info(f"Annotation written: {out_file}")
                    except subprocess.CalledProcessError as e:
                        logger.error(f"DIAMOND failed: {genome.name} vs {db.name}: {e}")
                        # The old TSV was judged stale above; never leave it to be read or
                        # certified by the caller's query stamp as if this search had run.
                        if out_file.exists():
                            out_file.unlink()
                            logger.error(f"Removed stale {out_file.name}; rerun the step to search it again.")
        report_progress(progress, len(genome_files), len(genome_files))


def build_diamond_db_for_files(
    fasta_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    threads: int,
) -> list[Path]:
    db_output_dir = output_path
    db_output_dir.mkdir(parents=True, exist_ok=True)

    if not fasta_files:
        logger.warning("No FASTA files provided.")
        return []

    created: list[Path] = []
    for fasta in fasta_files:
        if not fasta.is_file():
            logger.warning(f"Skipping missing file: {fasta}")
            continue

        db_name = fasta.stem
        final_db = db_output_dir / f"{db_name}.dmnd"

        if dmnd_is_current(fasta, final_db):
            logger.info(f"Skipping {db_name}: DB already exists.")
            created.append(final_db)
            continue
        if final_db.exists():
            logger.info(f"Rebuilding {db_name}: {fasta.name} changed since the DB was built.")

        with tempfile.TemporaryDirectory(dir=db_output_dir) as tmpdir:
            tmp_db = Path(tmpdir) / f"{db_name}.dmnd"
            cmd = [
                *tool_cmd("diamond"), "makedb",
                "--in", str(fasta),
                "--db", str(tmp_db),
                "--threads", str(threads),
            ]
            logger.info(f"Building DIAMOND DB for {fasta.name}")
            try:
                subprocess.run(cmd, check=True)
                shutil.move(str(tmp_db), final_db)
                write_dmnd_stamp(fasta, final_db)
                logger.info(f"Created DIAMOND DB: {final_db}")
                created.append(final_db)
            except subprocess.CalledProcessError as e:
                logger.error(f"DIAMOND failed for {fasta.name}: {e}")
    return created


def hit_method(db_name: str, full_names: Collection[str] = ()) -> str:
    """The Method of a database's hits: the name before its first '_' for the built-in
    catalog names of older runs ('CARD_AMR' -> 'CARD', 'AMRF_core-AMR' -> 'AMRF'), the
    whole name for a user database (``full_names``, pipeline.UserDatabases)."""
    return db_name if db_name in full_names else db_name.split("_")[0]


def parse_diamond_hit_file(
    path: Path,
    category: str,
    logger: logging.Logger | None = None,
    full_names: Collection[str] = (),
) -> pd.DataFrame:
    empty = pd.DataFrame(columns=[
        "Sample", "Contig", "Start", "Stop", "Strand",
        "Identity", "Evalue", "Bitscore", "Method", "MarkerCategory", "Product",
    ])
    try:
        if not path.exists():
            raise FileNotFoundError(path)
        if path.stat().st_size == 0:
            return empty
        t0 = time.perf_counter()
        print(f"[file] reading {path}")
        df = pd.read_csv(path, sep="\t", header=None)
        print(f"[file] read {path} in {time.perf_counter() - t0:.2f}s")
        if not isinstance(df, pd.DataFrame):
            raise ValueError("Parsed object is not a DataFrame")
        if df.shape[1] < 14:
            raise ValueError(f"Expected ≥14 columns (with qlen/slen), got {df.shape[1]}")

        base_cols = [
            "Query", "Product", "Identity", "Align_len", "Mismatches", "Gap_openings",
            "Q_start", "Q_end", "S_start", "S_end", "Evalue", "Bitscore",
        ]
        df = df.iloc[:, :14]
        df.columns = base_cols + ["Q_len", "S_len"]

        name_parts = path.stem.split("-vs-")
        sample = name_parts[0]
        method = name_parts[1] if len(name_parts) > 1 else "unknown"

        df["Sample"] = sample
        df["Method"] = hit_method(method, full_names)
        df["MarkerCategory"] = category if category in CATEGORIES else None

        df["Contig"] = df["Query"].astype(str).str.extract(r"contig_(\d+)", expand=False)
        df["Start"] = df[["Q_start", "Q_end"]].min(axis=1)
        df["Stop"] = df[["Q_start", "Q_end"]].max(axis=1)
        df["Strand"] = df.apply(lambda row: "+" if row["Q_start"] <= row["Q_end"] else "-", axis=1)
        df["Product"] = df["Product"]

        return df[
            ["Sample", "Contig", "Start", "Stop", "Strand",
             "Identity", "Evalue", "Bitscore", "Method", "MarkerCategory", "Product"]
        ]

    except Exception as e:
        if logger:
            logger.info(f"Failed to parse {path.name}: {e}")
        return empty


def load_diamond_tsvs(
    paths: list[Path],
    categories,
    logger: logging.Logger | None = None,
    full_names: Collection[str] = (),
) -> pd.DataFrame:
    cols = ["Sample", "Contig", "Start", "Stop", "Strand",
            "Identity", "Evalue", "Bitscore", "Method", "MarkerCategory", "Product"]

    frames = []
    for p in paths:
        try:
            df = parse_diamond_hit_file(p, categories[p.stem.split("-vs-")[1]], logger, full_names)
            if df is not None and len(df):
                frames.append(df)
        except Exception:
            if logger:
                logger.info("Failed parsing %s", p)

    if not frames:
        return pd.DataFrame(columns=cols).astype({
            "Sample": "string", "Contig": "Int64", "Start": "Int64", "Stop": "Int64",
            "Strand": "string", "Identity": "Float64",
            "Evalue": "Float64", "Bitscore": "Float64", "Method": "string",
            "MarkerCategory": "string", "Product": "string",
        })

    out = pd.concat(frames, ignore_index=True, sort=False)

    for c in cols:
        if c not in out.columns:
            out[c] = pd.NA

    out["Contig"] = pd.to_numeric(out["Contig"], errors="coerce").astype("Int64")
    out["Start"] = pd.to_numeric(out["Start"], errors="coerce").astype("Int64")
    out["Stop"] = pd.to_numeric(out["Stop"], errors="coerce").astype("Int64")
    out["Identity"] = pd.to_numeric(out["Identity"], errors="coerce").astype("Float64")
    out["Evalue"] = pd.to_numeric(out["Evalue"], errors="coerce").astype("Float64")
    out["Bitscore"] = pd.to_numeric(out["Bitscore"], errors="coerce").astype("Float64")
    out["Sample"] = out["Sample"].astype("string")
    out["Strand"] = out["Strand"].astype("string")
    out["Method"] = out["Method"].astype("string")
    out["MarkerCategory"] = out["MarkerCategory"].astype("string")
    out["Product"] = out["Product"].astype("string")

    out = out.dropna(subset=["Start", "Stop", "Contig", "Sample"])

    def _keep_best_by_overlap(g: pd.DataFrame) -> pd.DataFrame:
        if g.empty:
            return g
        g = g.sort_values(["Start", "Stop", "Identity", "Bitscore"], ascending=[True, True, False, False]).copy()
        keep = []
        cur_idx = []
        cur_start = cur_end = None

        for idx, row in g.iterrows():
            s, e = int(row["Start"]), int(row["Stop"])
            if cur_start is None:
                cur_start, cur_end = s, e
                cur_idx = [idx]
                continue

            if s <= cur_end:
                cur_end = max(cur_end, e)
                cur_idx.append(idx)
            else:
                block = g.loc[cur_idx]
                best = block.sort_values(
                    ["Identity", "Bitscore", "Stop"], ascending=[False, False, False]
                ).index[0]
                keep.append(best)
                cur_start, cur_end = s, e
                cur_idx = [idx]

        if cur_idx:
            block = g.loc[cur_idx]
            best = block.sort_values(
                ["Identity", "Bitscore", "Stop"], ascending=[False, False, False]
            ).index[0]
            keep.append(best)

        return g.loc[keep]

    parts = []
    for (s, c), g in out.groupby(["Sample", "Contig"], sort=False):
        parts.append(_keep_best_by_overlap(g))

    out = pd.concat(parts, ignore_index=True) if parts else out.iloc[0:0]
    return out[cols]


def read_annotations(path: Path, full_names: Collection[str] = ()) -> pd.DataFrame:
    """
    Read DIAMOND annotation TSVs produced for Bakta annotations vs MarkerSequence DBs
    (Method: `hit_method`, the whole name for the user databases in ``full_names``).
    """
    path = Path(path)

    if path.is_dir():
        tsv_files = sorted(path.rglob("*.tsv"))
    else:
        tsv_files = [path] if path.suffix == ".tsv" and path.exists() else []

    cols = [
        "AnnotationID", "Product",
        "Identity", "Align_len", "Mismatches", "Gap_openings",
        "Q_start", "Q_end", "S_start", "S_end",
        "Evalue", "Bitscore", "Q_len", "S_len", "Subject_coverage",
        "Method", "MarkerCategory", "Category", "DB",
    ]

    if not tsv_files:
        return pd.DataFrame(columns=cols).astype({
            "AnnotationID": "string",
            "Product": "string",
            "Identity": "Float64",
            "Align_len": "Int64",
            "Mismatches": "Int64",
            "Gap_openings": "Int64",
            "Q_start": "Int64",
            "Q_end": "Int64",
            "S_start": "Int64",
            "S_end": "Int64",
            "Evalue": "Float64",
            "Bitscore": "Float64",
            "Q_len": "Int64",
            "S_len": "Int64",
            "Subject_coverage": "Float64",
            "Method": "string",
            "MarkerCategory": "string",
            "Category": "string",
            "DB": "string",
        })

    frames: list[pd.DataFrame] = []

    for p in tsv_files:
        try:
            if p.stat().st_size == 0:
                continue
            t0 = time.perf_counter()
            print(f"[file] reading {p}")
            df = pd.read_csv(p, sep="\t", header=None, engine="python")
            print(f"[file] read {p} in {time.perf_counter() - t0:.2f}s")
            if df.shape[1] < 14:
                if logger := logging.getLogger("diamond_prep"):
                    logger.warning(
                        f"Skipping {p.name}: expected ≥14 columns (with qlen/slen), got {df.shape[1]}"
                    )
                continue

            base_cols = [
                "Query", "Product", "Identity", "Align_len", "Mismatches", "Gap_openings",
                "Q_start", "Q_end", "S_start", "S_end", "Evalue", "Bitscore",
            ]
            df = df.iloc[:, :14]
            df.columns = base_cols + ["Q_len", "S_len"]

            df["AnnotationID"] = df["Query"].astype(str)

            stem = p.stem
            after_vs = stem.split("-vs-", 1)[1] if "-vs-" in stem else stem

            parts = after_vs.split("_")
            category = parts[1] if len(parts) >= 2 and after_vs not in full_names else None
            method = hit_method(after_vs, full_names)

            marker_cat = category if category in CATEGORIES else None

            df["Method"] = str(method)
            df["MarkerCategory"] = marker_cat if marker_cat in CATEGORIES else pd.NA
            df["Category"] = df["MarkerCategory"]
            df["DB"] = after_vs

            df["Identity"] = pd.to_numeric(df["Identity"], errors="coerce")
            df["Align_len"] = pd.to_numeric(df["Align_len"], errors="coerce").astype("Int64")
            df["Mismatches"] = pd.to_numeric(df["Mismatches"], errors="coerce").astype("Int64")
            df["Gap_openings"] = pd.to_numeric(df["Gap_openings"], errors="coerce").astype("Int64")
            df["Q_start"] = pd.to_numeric(df["Q_start"], errors="coerce").astype("Int64")
            df["Q_end"] = pd.to_numeric(df["Q_end"], errors="coerce").astype("Int64")
            df["S_start"] = pd.to_numeric(df["S_start"], errors="coerce").astype("Int64")
            df["S_end"] = pd.to_numeric(df["S_end"], errors="coerce").astype("Int64")
            df["Evalue"] = pd.to_numeric(df["Evalue"], errors="coerce")
            df["Bitscore"] = pd.to_numeric(df["Bitscore"], errors="coerce")
            df["Q_len"] = pd.to_numeric(df["Q_len"], errors="coerce").astype("Int64")
            df["S_len"] = pd.to_numeric(df["S_len"], errors="coerce").astype("Int64")

            denom = pd.to_numeric(df["S_len"], errors="coerce")
            numer = pd.to_numeric(df["Align_len"], errors="coerce")
            cov = (numer / denom) * 100
            cov = cov.where(denom.notna() & (denom != 0))
            df["Subject_coverage"] = pd.to_numeric(cov, errors="coerce").astype("Float64")

            df["AnnotationID"] = df["AnnotationID"].astype("string")
            df["Product"] = df["Product"].astype("string")
            df["Method"] = df["Method"].astype("string")
            df["MarkerCategory"] = df["MarkerCategory"].astype("string")
            df["Category"] = df["Category"].astype("string")
            df["DB"] = df["DB"].astype("string")

            frames.append(df[cols])

        except Exception as e:
            if logger := logging.getLogger("diamond_prep"):
                logger.warning(f"Failed to parse DIAMOND annotations from {p.name}: {e}")

    if not frames:
        return pd.DataFrame(columns=cols)

    out = pd.concat(frames, ignore_index=True, sort=False)
    return out[cols]


def collapse_diamond_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each (AnnotationID, DB) keep a single best hit based on Identity, Bitscore, Align_len, and Evalue.

    Grouping is by DB (the full ``source_category`` database name), not by Method (the source
    prefix alone). A single source such as AMRF ships several category databases
    (AMRF_core-AMR, AMRF_plus-AMR, AMRF_plus-BMRG, AMRF_plus-VIRULENCE, ...) that all share
    Method="AMRF". Collapsing by Method would keep only the single best AMRF hit per gene and
    silently drop the others, merging the core/plus and AMR/BMRG/VIRULENCE distinctions and
    dropping resistance-metadata-bearing hits before enrichment. Collapsing by DB keeps one
    best representative per database, preserving every category.
    """
    if df.empty:
        return df

    df = df.copy()
    df["Evalue_tmp"] = pd.to_numeric(df["Evalue"], errors="coerce")
    df["Evalue_tmp"] = df["Evalue_tmp"].fillna(1e9)

    group_key = "DB" if "DB" in df.columns else "Method"

    df = df.sort_values(
        ["AnnotationID", group_key, "Identity", "Bitscore", "Align_len", "Evalue_tmp"],
        ascending=[True, True, False, False, False, True],
    )

    best = df.groupby(["AnnotationID", group_key], sort=False).head(1).copy()
    best = best.drop(columns=["Evalue_tmp"])
    best = best.reset_index(drop=True)
    return best


__all__ = [
    "dmnd_is_current",
    "hit_tsv_is_current",
    "write_dmnd_stamp",
    "run_diamond_annotations",
    "build_diamond_db_for_files",
    "hit_method",
    "parse_diamond_hit_file",
    "load_diamond_tsvs",
    "read_annotations",
    "collapse_diamond_annotations",
]
