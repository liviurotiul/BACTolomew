from __future__ import annotations

import gzip
import logging
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from pipeline.CoreTreeWorkflow import tree_stamp_path
from pipeline.RunSteps import ProgressCallback, report_progress
from pipeline.ToolEnv import tool_cmd
from pipeline.TreesIndex import LEGACY_SNIPPY_TREE_REL, SNIPPY_RAW_REL, write_snippy_stamp

# Panaroo output folders are recognized by their core alignment (see _find_panaroo_output_dir).
_CORE_ALIGNMENT_CANDIDATES = (
    "core_gene_alignment.aln",
    "core_gene_alignment.fa",
    "core_gene_alignment.fasta",
    "core_gene_alignment.aln.gz",
    "core_gene_alignment.fa.gz",
    "core_gene_alignment.fasta.gz",
)
_SNIPPY_CORE_ALIGNMENT_CANDIDATES = (
    "core.full.aln",
    "core.aln",
)


@dataclass(frozen=True)
class PhylogenyResult:
    reference_genome: Path
    snippy_tree: Path
    warnings: list[dict[str, str]] = field(default_factory=list)


MIN_RAXML_SEQUENCES = 4


def _cmd_available(tool: str, probe_args: list[str]) -> bool:
    for cmd in ([*tool_cmd(tool), *probe_args], [tool, *probe_args]):
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        except FileNotFoundError:
            continue
        if res.returncode == 0:
            return True
    return False


def _open_text_auto(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def _copy_alignment_plain(source: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if str(source).lower().endswith(".gz"):
        with gzip.open(source, "rb") as fin, dest.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(source, dest)


def _fasta_headers(path: Path) -> list[str]:
    headers: list[str] = []
    with _open_text_auto(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if header:
                headers.append(header)
    return headers


def _write_trivial_tree(headers: list[str], output_tree: Path, lengths: list[float] | None = None) -> Path:
    output_tree.parent.mkdir(parents=True, exist_ok=True)
    if not headers:
        raise RuntimeError("Cannot write tree: no taxa available.")
    if len(headers) == 1:
        newick = f"{headers[0]};\n"
    elif lengths is not None:
        newick = "(" + ",".join(f"{h}:{max(0.0, x):.8g}" for h, x in zip(headers, lengths)) + ");\n"
    else:
        newick = "(" + ",".join(f"{h}:0" for h in headers) + ");\n"
    output_tree.write_text(newick, encoding="utf-8")
    return output_tree


def _three_taxon_lengths(msa_path: Path) -> list[float]:
    """Branch lengths of the only unrooted 3-taxon tree from pairwise p-distances
    (columns called in both sequences): a = (d_ab + d_ac - d_bc) / 2, and so on."""
    seqs: list[np.ndarray] = []
    current: list[str] = []
    with _open_text_auto(msa_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    seqs.append(np.frombuffer("".join(current).upper().encode("ascii", "replace"), dtype=np.uint8))
                current = []
            else:
                current.append(line.strip())
    if current:
        seqs.append(np.frombuffer("".join(current).upper().encode("ascii", "replace"), dtype=np.uint8))
    called = [np.isin(seq, np.frombuffer(b"ACGT", dtype=np.uint8)) for seq in seqs]

    def p(i: int, j: int) -> float:
        both = called[i] & called[j]
        n = int(both.sum())
        return float((seqs[i][both] != seqs[j][both]).sum()) / n if n else 0.0

    ab, ac, bc = p(0, 1), p(0, 2), p(1, 2)
    return [(ab + ac - bc) / 2, (ab + bc - ac) / 2, (ac + bc - ab) / 2]


def _find_alignment_file(root: Path, candidates: tuple[str, ...]) -> Path | None:
    for name in candidates:
        path = root / name
        if path.is_file():
            return path
    return None


def _run_raxml_tree(
    msa_path: Path,
    output_tree: Path,
    *,
    logger: logging.Logger,
    threads: int,
    prefix_name: str,
) -> Path:
    headers = _fasta_headers(msa_path)
    if len(headers) < MIN_RAXML_SEQUENCES:
        # RAxML-NG refuses fewer than 4 sequences; up to 3 taxa there is one unrooted
        # topology, so the tree is written directly (3 taxa: lengths from p-distances).
        logger.info("Writing the only possible tree for %s (%d taxa, no RAxML-NG).", msa_path.name, len(headers))
        lengths = _three_taxon_lengths(msa_path) if len(headers) == 3 else None
        return _write_trivial_tree(headers, output_tree, lengths)

    output_tree.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output_tree.parent) as tmpdir:
        tmpdir = Path(tmpdir)
        plain_msa = tmpdir / msa_path.name.removesuffix(".gz")
        _copy_alignment_plain(msa_path, plain_msa)
        prefix = tmpdir / prefix_name
        cmd = [
            *tool_cmd("raxml-ng"),
            "--search",
            "--msa",
            str(plain_msa),
            "--model",
            "GTR+G",
            "--threads",
            str(max(1, int(threads))),
            "--seed",
            "1",
            "--prefix",
            str(prefix),
            "--redo",
        ]
        logger.info("Running RAxML-NG tree inference on %s", msa_path.name)
        logger.info("CMD: " + " ".join(cmd))
        subprocess.run(cmd, check=True)
        best_tree = Path(str(prefix) + ".raxml.bestTree")
        if not best_tree.is_file():
            raise RuntimeError(f"RAxML-NG did not produce bestTree output for {msa_path.name}")
        shutil.move(str(best_tree), str(output_tree))
    return output_tree


def _find_fasta_for_gff(gff: Path) -> Path | None:
    base = gff.name
    base_no_ext = base[:-5] if base.endswith(".gff3") else (base[:-4] if base.endswith(".gff") else gff.stem)
    cands = [
        gff.with_suffix(".fa"),
        gff.with_suffix(".fna"),
        gff.with_suffix(".fasta"),
        gff.parent / f"{gff.stem}.fa",
        gff.parent / f"{gff.stem}.fna",
        gff.parent / f"{gff.stem}.fasta",
        gff.parent / f"{base_no_ext}.fa",
        gff.parent / f"{base_no_ext}.fna",
        gff.parent / f"{base_no_ext}.fasta",
    ]
    for c in cands:
        if c.exists() and c.is_file():
            return c
    return None


def _sanitize_gff3(text: str) -> str:
    lines = text.splitlines()
    out = []
    for ln in lines:
        if ln.startswith("##FASTA"):
            break
        if (
            ln.startswith("##feature-ontology")
            or ln.startswith("##attribute-ontology")
            or ln.startswith("##source-ontology")
        ):
            continue
        out.append(ln)
    return "\n".join(out).rstrip("\n")


def _sanitize_fasta_comments(text: str) -> str:
    """
    Remove FASTA comment lines that Biopython now deprecates when present
    before the first sequence header, and ignore comment-style lines overall.
    """
    lines = text.splitlines()
    out: list[str] = []
    seen_header = False
    for ln in lines:
        s = ln.strip()
        if not s:
            # keep blank lines only after FASTA has started
            if seen_header:
                out.append("")
            continue
        if s.startswith(">"):
            seen_header = True
            out.append(ln)
            continue
        if s.startswith((";", "#", "!")):
            # drop comment lines both before and after first header
            continue
        if seen_header:
            out.append(ln)
        else:
            # Non-comment before first header is malformed FASTA; skip it.
            continue
    return "\n".join(out).rstrip("\n")


def _prepare_gff_inputs(
    gff_files: list[Path],
    *,
    prep_dir: Path,
    logger: logging.Logger,
    consumer_name: str,
    sanitize_fasta_comments: bool = False,
) -> list[Path]:
    prep_dir.mkdir(parents=True, exist_ok=True)
    prepared: list[Path] = []
    for gff in gff_files:
        dest = prep_dir / (gff.stem + ".gff")
        fasta = _find_fasta_for_gff(gff)

        try:
            with gff.open("r", encoding="utf-8", errors="ignore") as gin:
                cleaned = _sanitize_gff3(gin.read())
            with dest.open("w", encoding="utf-8") as gout:
                gout.write(cleaned)
                if fasta and fasta.exists():
                    with fasta.open("r", encoding="utf-8", errors="ignore") as fin:
                        fasta_txt = fin.read()
                        if sanitize_fasta_comments:
                            fasta_txt = _sanitize_fasta_comments(fasta_txt)
                        gout.write("\n##FASTA\n")
                        gout.write(fasta_txt)
                else:
                    logger.warning(
                        f"No matching FASTA for {gff.name}; writing GFF without FASTA ({consumer_name} may fail)."
                    )
            prepared.append(dest)
            logger.info(f"Prepared for {consumer_name}: {dest.name}")
        except Exception as e:
            logger.error(f"Failed preparing {gff.name} for {consumer_name}: {e}")

    return [p for p in prepared if p.exists()]


def _find_panaroo_output_dir(search_root: Path, *, preferred_dir: Path | None = None) -> Path | None:
    candidates = {p.parent for p in search_root.rglob("gene_presence_absence.csv") if p.is_file()}
    if not candidates:
        return None

    def _rank(path: Path) -> tuple[int, int, str]:
        pref = 1
        if preferred_dir is not None and path == preferred_dir:
            pref = 0
        try:
            depth = len(path.relative_to(search_root).parts)
        except ValueError:
            depth = 999
        has_core = any((path / name).is_file() for name in _CORE_ALIGNMENT_CANDIDATES)
        return (pref, 0 if has_core else 1, depth, str(path))

    return sorted(candidates, key=_rank)[0]


def ensure_panaroo_out(output_path: Path, logger: logging.Logger) -> Path | None:
    panaroo_root = output_path / "panaroo"
    out_dir = panaroo_root / "out"
    if not panaroo_root.exists():
        return None

    source_dir = _find_panaroo_output_dir(panaroo_root, preferred_dir=out_dir)
    if source_dir is None:
        return None

    if source_dir != out_dir:
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        for item in source_dir.iterdir():
            target = out_dir / item.name
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target, ignore_errors=True)
                else:
                    target.unlink()
            shutil.move(str(item), str(target))
        logger.info(f"Panaroo outputs consolidated: {source_dir} -> {out_dir}")
    else:
        out_dir.mkdir(parents=True, exist_ok=True)

    if not (out_dir / "gene_presence_absence.csv").is_file():
        return None
    return out_dir


def run_panaroo(gff_files: list[Path], output_path: Path, logger: logging.Logger, threads: int = 1) -> Path | None:
    if not gff_files:
        logger.error("No GFF files provided to run_panaroo.")
        return None
    gff_files = [p for p in gff_files if p.exists() and p.is_file()]
    if not gff_files:
        logger.error("All provided GFF paths are missing.")
        return None
    if not (_cmd_available("panaroo", ["--help"]) or _cmd_available("panaroo", ["--version"])):
        logger.error("panaroo not available (configured env or PATH).")
        return None

    panaroo_root = output_path / "panaroo"
    if panaroo_root.exists() and panaroo_root.is_dir():
        for item in panaroo_root.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()

    prep_dir = panaroo_root / "gff_prepared"
    prepared = _prepare_gff_inputs(
        gff_files,
        prep_dir=prep_dir,
        logger=logger,
        consumer_name="Panaroo",
        sanitize_fasta_comments=True,
    )
    if not prepared:
        logger.error("No valid prepared GFFs for Panaroo.")
        return None

    with tempfile.TemporaryDirectory(dir=panaroo_root) as tmpdir:
        tmp_out = Path(tmpdir) / "out"
        cmd = [
            *tool_cmd("panaroo"),
            "-i",
            *map(str, prepared),
            "-o",
            str(tmp_out),
            "-t",
            str(threads),
            "--clean-mode",
            "strict",
            "--remove-invalid-genes",
            "-a",
            "core",
        ]
        logger.info(f"Running Panaroo on {len(prepared)} files -> {tmp_out}")
        logger.info("CMD: " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
            logger.info("Panaroo completed.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Panaroo failed: {e}")
            return None

        out_dir = ensure_panaroo_out(output_path, logger)
        if out_dir is None:
            logger.error("Panaroo completed but no gene_presence_absence.csv was found in Panaroo output.")
            return None
        return out_dir


def _default_reference_genome(genome_files: list[Path]) -> Path:
    if not genome_files:
        raise RuntimeError("No genome files available for phylogeny.")
    return sorted(genome_files, key=lambda p: (p.stem, str(p)))[0]


def _run_snippy_per_sample(
    genome_files: list[Path],
    *,
    reference_genome: Path,
    per_sample_dir: Path,
    logger: logging.Logger,
    threads: int,
    progress: ProgressCallback = None,
) -> list[Path]:
    sample_dirs: list[Path] = []
    per_sample_dir.mkdir(parents=True, exist_ok=True)
    queries = [
        genome
        for genome in sorted(genome_files, key=lambda p: (p.stem, str(p)))
        if genome.resolve() != reference_genome.resolve()
    ]
    for done, genome in enumerate(queries):
        report_progress(progress, done, len(queries))
        sample_dir = per_sample_dir / genome.stem
        if sample_dir.exists():
            shutil.rmtree(sample_dir, ignore_errors=True)
        cmd = [
            *tool_cmd("snippy"),
            "--cpus",
            str(max(1, int(threads))),
            "--outdir",
            str(sample_dir),
            "--ref",
            str(reference_genome),
            "--ctgs",
            str(genome),
        ]
        logger.info("Running Snippy on %s vs %s", genome.name, reference_genome.name)
        logger.info("CMD: " + " ".join(cmd))
        completed = subprocess.run(cmd, capture_output=True, text=True)
        if completed.stdout:
            logger.info("Snippy stdout for %s:\n%s", genome.name, completed.stdout.rstrip())
        if completed.stderr:
            logger.info("Snippy stderr for %s:\n%s", genome.name, completed.stderr.rstrip())
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout).strip()
            if detail:
                raise RuntimeError(
                    f"Snippy failed for {genome.name} against {reference_genome.name}: {detail}"
                )
            raise subprocess.CalledProcessError(completed.returncode, cmd)
        if not (sample_dir / "snps.vcf").is_file():
            raise RuntimeError(f"Snippy did not produce snps.vcf for {genome.name}")
        sample_dirs.append(sample_dir)
    report_progress(progress, len(queries), len(queries))
    return sample_dirs


def _run_snippy_core(
    *,
    reference_genome: Path,
    sample_dirs: list[Path],
    core_dir: Path,
    logger: logging.Logger,
) -> Path:
    core_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        *tool_cmd("snippy-core", workdir=str(core_dir)),
        "--ref",
        str(reference_genome),
        "--prefix",
        "core",
        *map(str, sample_dirs),
    ]
    logger.info("Running snippy-core on %d sample(s).", len(sample_dirs))
    logger.info("CMD: " + " ".join(cmd))
    subprocess.run(cmd, cwd=core_dir, check=True)
    alignment = _find_alignment_file(core_dir, _SNIPPY_CORE_ALIGNMENT_CANDIDATES)
    if alignment is None:
        raise RuntimeError("snippy-core did not produce a core alignment.")
    return alignment


def raxml_ng_version() -> str | None:
    try:
        result = subprocess.run([*tool_cmd("raxml-ng"), "--version"], capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return None
    found = re.search(r"RAxML-NG v\.\s*(\S+)", (result.stdout or "") + (result.stderr or ""))
    return found.group(1) if found else None


def run_phylogeny(
    genome_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    snippy_threads: int = 1,
    raxml_threads: int = 1,
    reference_genome: Path | None = None,
    progress: ProgressCallback = None,
) -> PhylogenyResult:
    """Snippy against the reference, snippy-core, then RAxML-NG GTR+G on core.full.aln.

    The raw (unrooted) tree goes to trees/raw/snippy_core.raxml.nwk with its stamp
    trees/snippy_core.json; pipeline.TreesIndex roots it, renames its Reference tip and
    lists it as the secondary, not recombination-corrected tree. The core genome tree is
    the core_tree step (IQ-TREE on the Gubbins output), not this one.
    """
    if not (_cmd_available("raxml-ng", ["--help"]) or _cmd_available("raxml-ng", ["--version"])):
        raise RuntimeError("raxml-ng not available (configured env or PATH).")
    if not (_cmd_available("snippy", ["--help"]) or _cmd_available("snippy", ["--version"])):
        raise RuntimeError("snippy not available (configured env or PATH).")
    if not (_cmd_available("snippy-core", ["--help"]) or _cmd_available("snippy-core", ["--version"])):
        raise RuntimeError("snippy-core not available (configured env or PATH).")

    reference_genome = (reference_genome or _default_reference_genome(genome_files)).resolve()
    snippy_tree = output_path / SNIPPY_RAW_REL
    snippy_tree.parent.mkdir(parents=True, exist_ok=True)
    tree_stamp_path(output_path, "snippy_core").unlink(missing_ok=True)

    snippy_root = output_path / "snippy"
    with tempfile.TemporaryDirectory(dir=output_path) as tmpdir:
        tmpdir = Path(tmpdir)
        tmp_snippy = tmpdir / "snippy"
        tmp_per_sample = tmp_snippy / "per_sample"
        tmp_core = tmp_snippy / "core"

        sample_dirs = _run_snippy_per_sample(
            genome_files,
            reference_genome=reference_genome,
            per_sample_dir=tmp_per_sample,
            logger=logger,
            threads=snippy_threads,
            progress=progress,
        )

        taxa = 1
        if not sample_dirs:
            _write_trivial_tree([reference_genome.stem], snippy_tree)
            tmp_core.mkdir(parents=True, exist_ok=True)
        else:
            core_alignment = _run_snippy_core(
                reference_genome=reference_genome,
                sample_dirs=sample_dirs,
                core_dir=tmp_core,
                logger=logger,
            )
            taxa = len(_fasta_headers(core_alignment))
            _run_raxml_tree(
                core_alignment,
                snippy_tree,
                logger=logger,
                threads=raxml_threads,
                prefix_name="snippy_raxml",
            )

        ref_txt = tmp_snippy / "reference.txt"
        ref_txt.parent.mkdir(parents=True, exist_ok=True)
        ref_txt.write_text(str(reference_genome) + "\n", encoding="utf-8")

        if snippy_root.exists():
            shutil.rmtree(snippy_root, ignore_errors=True)
        shutil.move(str(tmp_snippy), str(snippy_root))

    (output_path / LEGACY_SNIPPY_TREE_REL).unlink(missing_ok=True)
    warnings: list[dict[str, str]] = []
    if taxa < MIN_RAXML_SEQUENCES:
        note = (
            f"The Snippy alignment has {taxa} sequence(s); RAxML-NG needs at least "
            f"{MIN_RAXML_SEQUENCES}, so the only possible tree was written"
            + (" with branch lengths from pairwise p-distances." if taxa == 3 else ".")
        )
        logger.warning(note)
        warnings.append({"step": "phylogeny", "code": "snippy_tree_few_genomes", "message": note})
        write_snippy_stamp(output_path, tool_version=None, threads=None, small_tree_note=note)
    else:
        write_snippy_stamp(output_path, tool_version=raxml_ng_version(), threads=raxml_threads)
    return PhylogenyResult(
        reference_genome=reference_genome,
        snippy_tree=snippy_tree,
        warnings=warnings,
    )


__all__ = [
    "run_panaroo",
    "ensure_panaroo_out",
    "run_phylogeny",
    "PhylogenyResult",
]
