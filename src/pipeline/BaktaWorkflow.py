from __future__ import annotations

import json
import logging
import time
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from pipeline.FastaUtils import prepare_query_fasta, subseq
from pipeline.ToolEnv import tool_cmd


def run_bakta_on_fastas(
    fasta_paths: list[Path],
    output_path: Path,
    logger: logging.Logger,
    bakta_db: str = "",
    bakta_threads: int = 1,
):
    """
    Run Bakta on sanitized copies of the input FASTAs.
    """
    bakta_final_dir = output_path / "bakta"
    bakta_final_dir.mkdir(parents=True, exist_ok=True)

    tmp_root = output_path / "bakta_tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)

    completed: list[str] = []
    failures: list[str] = []

    for fasta in fasta_paths:
        if not fasta.exists() or not fasta.is_file():
            msg = f"non-existent input FASTA: {fasta}"
            logger.warning("Skipping %s", msg)
            failures.append(msg)
            continue

        try:
            query_fa = prepare_query_fasta(fasta, tmp_root, logger, write_map=False)
        except Exception as e:
            failures.append(f"{fasta.name}: failed to prepare FASTA ({e})")
            logger.error(f"Failed to prepare FASTA for Bakta: {fasta} — {e}")
            continue

        prefix = query_fa.stem
        sample_final_out = bakta_final_dir / prefix
        result_file = sample_final_out / f"{prefix}.json"

        if result_file.exists():
            logger.info(f"Skipping {prefix}: already processed.")
            completed.append(prefix)
            continue

        tmp_output = tmp_root / f"{prefix}__bakta_out"
        if tmp_output.exists():
            shutil.rmtree(tmp_output, ignore_errors=True)
        tmp_output.mkdir(parents=True, exist_ok=True)

        cmd = [*tool_cmd("bakta")]
        if bakta_db:
            cmd.extend(["--db", bakta_db])
        cmd.extend([
            "--force",
            "--output", str(tmp_output),
            "--threads", str(bakta_threads),
            "--prefix", prefix,
            "--keep-contig-headers",
            "--skip-trna",
            "--skip-tmrna",
            "--skip-rrna",
            "--skip-ncrna",
            "--skip-ncrna-region",
            "--skip-crispr",
            "--skip-sorf",
            "--skip-gap",
            "--skip-ori",
            "--skip-plot",
            "--skip-pseudo",
            str(query_fa),
        ])

        logger.info(f"Running Bakta on {fasta.name} (sanitized: {query_fa.name}) → temp: {tmp_output}")
        logger.info("CMD: " + " ".join(cmd))

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            failures.append(f"{prefix}: bakta exited with code {e.returncode}")
            logger.error(f"Bakta failed on {fasta.name}: {e}")
            continue

        try:
            shutil.copytree(tmp_output, sample_final_out, dirs_exist_ok=True)
            shutil.copy2(query_fa, sample_final_out / f"{prefix}.fna")
            logger.info(f"Finished: {fasta.name} → {sample_final_out}")
            maybe_json = sample_final_out / f"{prefix}.json"
            if not maybe_json.exists():
                msg = f"{prefix}: expected JSON not found in {sample_final_out}"
                failures.append(msg)
                logger.warning(f"Expected JSON not found: {maybe_json.name} in {sample_final_out}")
                continue
            completed.append(prefix)
        except Exception as e:
            failures.append(f"{prefix}: failed to copy Bakta results ({e})")
            logger.error(f"Failed to copy Bakta results for {fasta.name}: {e}")
        finally:
            shutil.rmtree(tmp_output, ignore_errors=True)

    if failures:
        preview = "; ".join(failures[:5])
        if len(failures) > 5:
            preview += f"; ... ({len(failures)} failures total)"
        raise RuntimeError(
            f"Bakta did not produce usable outputs for {len(failures)} genome(s). "
            f"Successful: {len(completed)}. Failures: {preview}"
        )


def parse_json_file(path: Path, logger: logging.Logger):
    try:
        with path.open("r", encoding="utf-8") as fh:
            obj = json.load(fh)

        sample_name = path.stem
        sample_fa = path.parent / f"{sample_name}.fna"

        feats = obj.get("features") or []
        rows = []
        for feat in feats:
            if feat.get("type") != "cds":
                continue

            contig = feat.get("sequence")
            m = re.search(r"contig_(\d+)", str(contig)) if contig else None
            contig_val = int(m.group(1)) if m else None

            start = feat.get("start")
            stop = feat.get("stop")
            strand = feat.get("strand")

            nt_seq = subseq(sample_fa, str(contig), start, stop, strand)

            rows.append({
                "Sample": sample_name,
                "Contig": contig_val,
                "Id": feat.get("id") or None,
                "Gene": feat.get("gene") or None,
                "Start": start,
                "Stop": stop,
                "Strand": strand or None,
                "Product": feat.get("product") or None,
                "Locus": feat.get("locus") or None,
                "StartType": feat.get("start_type") or None,
                "NCBI_NRP_ID": feat.get("ups", {}).get("ncbi_nrp_id"),
                "Refseq_ID": feat.get("ups", {}).get("RefSeq"),
                "RBS": feat.get("rbs_motif") or None,
                "Protein_MW": feat.get("seq_stats", {}).get("molecular_weight"),
                "pI": feat.get("seq_stats", {}).get("isoelectric_point"),
                "Hypothetical": bool(feat.get("hypothetical", False)),
                "Nucleotide_Sequence": nt_seq,
                "Amino_Acid_Sequence": feat.get("aa") or None,
            })

        df = pd.DataFrame(rows)
        if df.empty:
            return df

        df["Contig"] = pd.to_numeric(df["Contig"], errors="coerce")
        df["Protein_MW"] = pd.to_numeric(df["Protein_MW"], errors="coerce")
        df["pI"] = pd.to_numeric(df["pI"], errors="coerce")
        return df

    except Exception as e:
        logger.warning(f"Failed to process {path.name}: {e}")
        return pd.DataFrame()


def parse_multiple_bakta_json(json_paths: list[Path], logger: logging.Logger) -> pd.DataFrame:
    frames = []
    for p in json_paths:
        try:
            t0 = time.perf_counter()
            print(f"[file] reading {p}")
            df = parse_json_file(p, logger)
            print(f"[file] read {p} in {time.perf_counter() - t0:.2f}s")
            if df is not None and len(df):
                frames.append(df)
        except Exception as e:
            logger.exception("Failed parsing %s", p)
    if not frames:
        cols = ["Sample", "Contig", "Id", "Gene", "Start", "Stop", "Strand", "Product",
                "Locus", "StartType", "NCBI_NRP_ID", "Refseq_ID", "RBS", "Protein_MW", "pI", "Hypothetical",
                "Nucleotide_Sequence", "Amino_Acid_Sequence"]
        return pd.DataFrame(columns=cols).astype({
            "Sample": "string", "Contig": "Int64", "Id": "string", "Gene": "string",
            "Start": "Int64", "Stop": "Int64", "Strand": "string", "Product": "string",
            "Locus": "string", "StartType": "string", "NCBI_NRP_ID": "string",
            "Refseq_ID": "string", "RBS": "string", "Protein_MW": "Float64", "pI": "Float64",
            "Hypothetical": "boolean", "Nucleotide_Sequence": "string", "Amino_Acid_Sequence": "string",
        })

    out = pd.concat(frames, ignore_index=True, sort=False)

    for col in ["Sample", "Contig", "Id", "Gene", "Start", "Stop", "Strand", "Product",
                "Locus", "StartType", "NCBI_NRP_ID", "Refseq_ID", "RBS", "Protein_MW", "pI", "Hypothetical",
                "Nucleotide_Sequence", "Amino_Acid_Sequence"]:
        if col not in out.columns:
            out[col] = pd.NA

    out["Contig"] = pd.to_numeric(out["Contig"], errors="coerce").astype("Int64")
    out["Start"] = pd.to_numeric(out["Start"], errors="coerce").astype("Int64")
    out["Stop"] = pd.to_numeric(out["Stop"], errors="coerce").astype("Int64")
    out["Protein_MW"] = pd.to_numeric(out["Protein_MW"], errors="coerce").astype("Float64")
    out["pI"] = pd.to_numeric(out["pI"], errors="coerce").astype("Float64")
    out["Hypothetical"] = out["Hypothetical"].astype("boolean")

    for col in ["Sample", "Id", "Gene", "Strand", "Product", "Locus", "StartType",
                "NCBI_NRP_ID", "Refseq_ID", "RBS", "Nucleotide_Sequence", "Amino_Acid_Sequence"]:
        out[col] = out[col].astype("string")

    return out[["Sample", "Contig", "Id", "Gene", "Start", "Stop", "Strand", "Product",
                "Locus", "StartType", "NCBI_NRP_ID", "Refseq_ID", "RBS", "Protein_MW", "pI", "Hypothetical",
                "Nucleotide_Sequence", "Amino_Acid_Sequence"]]


__all__ = [
    "run_bakta_on_fastas",
    "parse_json_file",
    "parse_multiple_bakta_json",
]
