from __future__ import annotations

import csv
import json
import logging
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from pipeline.ToolEnv import configured_assembly_qc_settings, tool_cmd


ASSEMBLY_QC_OUT_FILES: tuple[str, ...] = (
    "summary.tsv",
    "summary.json",
    "quast/report.tsv",
    "checkm2/quality_report.tsv",
)

ASSEMBLY_QC_COLUMNS: list[str] = [
    "Sample",
    "Quality",
    "QuastTotalLength",
    "QuastContigs",
    "QuastLargestContig",
    "QuastN50",
    "QuastL50",
    "QuastGC",
    "QuastNs",
    "QuastNsPer100kbp",
    "CheckM2Completeness",
    "CheckM2Contamination",
    "CheckM2CompletenessModel",
    "QualityNotes",
]


@dataclass(frozen=True)
class AssemblyQualityResult:
    output_dir: Path
    summary_tsv: Path
    summary_json: Path
    quast_report_tsv: Path
    checkm2_report_tsv: Path
    sample_count: int


def _clean_text(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"} else text


def _to_float(value: Any) -> float | None:
    text = _clean_text(value).replace(",", "")
    if not text:
        return None
    try:
        val = float(text)
    except Exception:
        return None
    return val if math.isfinite(val) else None


def _to_int(value: Any) -> int | None:
    val = _to_float(value)
    if val is None:
        return None
    return int(round(val))


def _fasta_metrics(path: Path) -> dict[str, float | int | None]:
    lengths: list[int] = []
    total_gc = 0
    total_n = 0
    total = 0
    current_len = 0
    current_gc = 0
    current_n = 0

    def flush_current() -> None:
        nonlocal current_len, current_gc, current_n, total, total_gc, total_n
        if current_len <= 0:
            return
        lengths.append(current_len)
        total += current_len
        total_gc += current_gc
        total_n += current_n
        current_len = 0
        current_gc = 0
        current_n = 0

    with Path(path).open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_current()
                continue
            seq = line.upper()
            current_len += len(seq)
            current_gc += seq.count("G") + seq.count("C")
            current_n += seq.count("N")
    flush_current()

    if not lengths:
        return {
            "contigs": None,
            "largest_contig": None,
            "total_length": None,
            "gc": None,
            "n50": None,
            "l50": None,
            "ns": None,
            "ns_per_100kbp": None,
        }

    lengths_sorted = sorted(lengths, reverse=True)
    half = total / 2.0
    cumulative = 0
    n50 = lengths_sorted[-1]
    l50 = len(lengths_sorted)
    for idx, length in enumerate(lengths_sorted, start=1):
        cumulative += length
        if cumulative >= half:
            n50 = length
            l50 = idx
            break

    return {
        "contigs": len(lengths),
        "largest_contig": max(lengths),
        "total_length": total,
        "gc": (100.0 * total_gc / total) if total else None,
        "n50": n50,
        "l50": l50,
        "ns": total_n,
        "ns_per_100kbp": (100000.0 * total_n / total) if total else None,
    }


def _metric_lookup(metrics: Mapping[str, str], *names: str) -> str:
    normalized = {
        str(key).strip().casefold(): str(value).strip()
        for key, value in metrics.items()
    }
    for name in names:
        key = str(name).strip().casefold()
        if key in normalized:
            return normalized[key]
    return ""


def _load_quast_metrics(report_tsv: Path) -> dict[str, dict[str, str]]:
    if not report_tsv.is_file():
        return {}
    try:
        df = pd.read_csv(report_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}
    if df.empty or len(df.columns) < 2:
        return {}
    metric_col = df.columns[0]
    out: dict[str, dict[str, str]] = {}
    for sample_col in df.columns[1:]:
        sample = str(sample_col).strip()
        if not sample:
            continue
        metrics: dict[str, str] = {}
        for _idx, row in df.iterrows():
            metric = str(row.get(metric_col) or "").strip()
            if metric:
                metrics[metric] = str(row.get(sample_col) or "").strip()
        out[sample] = metrics
    return out


def _load_checkm2_metrics(report_tsv: Path) -> dict[str, dict[str, str]]:
    if not report_tsv.is_file():
        return {}
    try:
        df = pd.read_csv(report_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}
    if df.empty:
        return {}
    name_col = next((col for col in df.columns if col.strip().casefold() in {"name", "sample", "bin id", "bin"}), None)
    if name_col is None:
        name_col = df.columns[0]
    out: dict[str, dict[str, str]] = {}
    for row in df.to_dict(orient="records"):
        sample = _clean_text(row.get(name_col))
        if sample:
            metrics = {str(k): str(v) for k, v in row.items()}
            for alias in _sample_aliases(sample):
                out[alias] = metrics
    return out


def _sample_aliases(sample: str) -> list[str]:
    raw = _clean_text(sample)
    if not raw:
        return []
    candidates = [raw, Path(raw).name, Path(raw).stem]
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if raw.endswith(suffix):
            candidates.append(raw[: -len(suffix)])
        name = Path(raw).name
        if name.endswith(suffix):
            candidates.append(name[: -len(suffix)])

    out: list[str] = []
    seen: set[str] = set()
    for value in candidates:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def _checkm2_value(metrics: Mapping[str, str], *names: str) -> str:
    lower = {str(k).strip().casefold(): str(v).strip() for k, v in metrics.items()}
    for name in names:
        key = str(name).strip().casefold()
        if key in lower:
            return lower[key]
    return ""


def _median(values: list[float]) -> float | None:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return (clean[mid - 1] + clean[mid]) / 2.0


def _quality_score(row: Mapping[str, Any], *, median_total: float | None, median_gc: float | None) -> tuple[int, list[str]]:
    notes: list[str] = []
    completeness = _to_float(row.get("CheckM2Completeness"))
    contamination = _to_float(row.get("CheckM2Contamination")) or 0.0
    contigs = _to_float(row.get("QuastContigs"))
    n50 = _to_float(row.get("QuastN50"))
    total = _to_float(row.get("QuastTotalLength"))
    gc = _to_float(row.get("QuastGC"))
    ns_per_100kbp = _to_float(row.get("QuastNsPer100kbp"))

    score = completeness if completeness is not None else 0.0
    if completeness is None:
        notes.append("missing CheckM2 completeness")

    contamination_penalty = min(45.0, max(0.0, contamination) * 5.0)
    if contamination_penalty:
        notes.append(f"contamination penalty {contamination_penalty:.1f}")
    score -= contamination_penalty

    if contigs is not None and contigs > 80:
        penalty = min(18.0, (contigs - 80.0) / 28.0)
        score -= penalty
        notes.append(f"fragmentation penalty {penalty:.1f}")

    if total and n50 and total > 0:
        n50_ratio = n50 / total
        if n50_ratio < 0.02:
            penalty = min(12.0, ((0.02 - n50_ratio) / 0.02) * 12.0)
            score -= penalty
            notes.append(f"low N50 penalty {penalty:.1f}")

    if ns_per_100kbp is not None and ns_per_100kbp > 0:
        penalty = min(12.0, ns_per_100kbp / 40.0)
        score -= penalty
        notes.append(f"N content penalty {penalty:.1f}")

    if total and median_total and median_total > 0:
        rel_delta = abs(total - median_total) / median_total
        if rel_delta > 0.20:
            penalty = min(10.0, ((rel_delta - 0.20) / 0.50) * 10.0)
            score -= penalty
            notes.append(f"length outlier penalty {penalty:.1f}")

    if gc is not None and median_gc is not None:
        delta = abs(gc - median_gc)
        if delta > 5.0:
            penalty = min(8.0, ((delta - 5.0) / 10.0) * 8.0)
            score -= penalty
            notes.append(f"GC outlier penalty {penalty:.1f}")

    return int(round(max(0.0, min(100.0, score)))), notes


def _write_summary(root: Path, rows: list[dict[str, Any]]) -> tuple[Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    tsv_path = root / "summary.tsv"
    json_path = root / "summary.json"
    df = pd.DataFrame(rows, columns=ASSEMBLY_QC_COLUMNS)
    df.to_csv(tsv_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    json_path.write_text(
        json.dumps({"sample_count": len(rows), "samples": rows}, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return tsv_path, json_path


def _run_command(cmd: list[str], *, stdout_log: Path, stderr_log: Path, logger: logging.Logger) -> None:
    logger.info("CMD: %s", " ".join(str(part) for part in cmd))
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    with stdout_log.open("w", encoding="utf-8") as out_fh, stderr_log.open("w", encoding="utf-8") as err_fh:
        subprocess.run(cmd, check=True, stdout=out_fh, stderr=err_fh)


def run_assembly_quality(
    genome_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int,
    config: Mapping[str, Any] | None = None,
) -> AssemblyQualityResult:
    genomes = [Path(path) for path in genome_files if Path(path).is_file()]
    if not genomes:
        raise RuntimeError("No genome FASTA files provided for assembly QC.")

    settings = configured_assembly_qc_settings(config)
    root = output_path / "assembly_qc"
    quast_dir = root / "quast"
    checkm2_dir = root / "checkm2"
    checkm2_input_dir = root / "checkm2_input"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True, exist_ok=True)
    checkm2_input_dir.mkdir(parents=True, exist_ok=True)
    for genome in genomes:
        shutil.copy2(genome, checkm2_input_dir / genome.name)

    labels = [genome.stem for genome in genomes]
    quast_cmd = [
        *tool_cmd("quast", executable="quast.py"),
        *[str(genome) for genome in genomes],
        "-o",
        str(quast_dir),
        "--threads",
        str(max(1, int(threads))),
        "--min-contig",
        str(int(settings.get("quast_min_contig") or 500)),
        "--labels",
        ",".join(labels),
        "--no-plots",
        "--no-html",
    ]
    quast_cmd.extend(str(arg) for arg in settings.get("quast_extra_args") or [] if str(arg).strip())
    logger.info("Running QUAST assembly QC on %d genome(s).", len(genomes))
    _run_command(
        quast_cmd,
        stdout_log=root / "quast.command.stdout.log",
        stderr_log=root / "quast.command.stderr.log",
        logger=logger,
    )

    checkm2_cmd = [
        *tool_cmd("checkm2"),
        "predict",
        "--threads",
        str(max(1, int(threads))),
        "--input",
        str(checkm2_input_dir),
        "--output-directory",
        str(checkm2_dir),
    ]
    db_path = _clean_text(settings.get("checkm2_database_path"))
    if db_path:
        checkm2_cmd.extend(["--database_path", db_path])
    if bool(settings.get("checkm2_lowmem")):
        checkm2_cmd.append("--lowmem")
    checkm2_cmd.extend(str(arg) for arg in settings.get("checkm2_extra_args") or [] if str(arg).strip())
    logger.info("Running CheckM2 genome quality prediction on %d genome(s).", len(genomes))
    _run_command(
        checkm2_cmd,
        stdout_log=root / "checkm2.command.stdout.log",
        stderr_log=root / "checkm2.command.stderr.log",
        logger=logger,
    )

    quast_metrics = _load_quast_metrics(quast_dir / "report.tsv")
    checkm2_metrics = _load_checkm2_metrics(checkm2_dir / "quality_report.tsv")
    fasta_metrics = {genome.stem: _fasta_metrics(genome) for genome in genomes}

    rows: list[dict[str, Any]] = []
    for genome in genomes:
        sample = genome.stem
        q = quast_metrics.get(sample, {})
        c = next((checkm2_metrics[alias] for alias in _sample_aliases(sample) if alias in checkm2_metrics), {})
        f = fasta_metrics.get(sample, {})
        ns = _to_int(_metric_lookup(q, "# N's")) if q else None
        total = _to_int(_metric_lookup(q, "Total length", "Total length (>= 0 bp)")) if q else None
        ns_per_100kbp = _to_float(_metric_lookup(q, "# N's per 100 kbp", "# N's per 100 kbp")) if q else None
        if ns_per_100kbp is None and ns is not None and total:
            ns_per_100kbp = 100000.0 * float(ns) / float(total)
        row = {
            "Sample": sample,
            "Quality": "",
            "QuastTotalLength": total if total is not None else f.get("total_length"),
            "QuastContigs": (
                _to_int(_metric_lookup(q, "# contigs", "# contigs (>= 0 bp)")) if q else None
            ) or f.get("contigs"),
            "QuastLargestContig": (
                _to_int(_metric_lookup(q, "Largest contig")) if q else None
            ) or f.get("largest_contig"),
            "QuastN50": (_to_int(_metric_lookup(q, "N50")) if q else None) or f.get("n50"),
            "QuastL50": (_to_int(_metric_lookup(q, "L50")) if q else None) or f.get("l50"),
            "QuastGC": (_to_float(_metric_lookup(q, "GC (%)")) if q else None) or f.get("gc"),
            "QuastNs": ns if ns is not None else f.get("ns"),
            "QuastNsPer100kbp": ns_per_100kbp if ns_per_100kbp is not None else f.get("ns_per_100kbp"),
            "CheckM2Completeness": _to_float(_checkm2_value(c, "Completeness")),
            "CheckM2Contamination": _to_float(_checkm2_value(c, "Contamination")),
            "CheckM2CompletenessModel": _checkm2_value(c, "Completeness_Model_Used", "Completeness Model Used", "Model"),
            "QualityNotes": "",
        }
        rows.append(row)

    median_total = _median([float(v) for row in rows if (v := _to_float(row.get("QuastTotalLength"))) is not None])
    median_gc = _median([float(v) for row in rows if (v := _to_float(row.get("QuastGC"))) is not None])
    for row in rows:
        quality, notes = _quality_score(row, median_total=median_total, median_gc=median_gc)
        row["Quality"] = quality
        row["QualityNotes"] = "; ".join(notes)

    summary_tsv, summary_json = _write_summary(root, rows)
    return AssemblyQualityResult(
        output_dir=root,
        summary_tsv=summary_tsv,
        summary_json=summary_json,
        quast_report_tsv=quast_dir / "report.tsv",
        checkm2_report_tsv=checkm2_dir / "quality_report.tsv",
        sample_count=len(rows),
    )


def load_assembly_quality_summary(output_path: Path) -> pd.DataFrame:
    path = Path(output_path) / "assembly_qc" / "summary.tsv"
    if not path.is_file():
        return pd.DataFrame(columns=ASSEMBLY_QC_COLUMNS)
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return pd.DataFrame(columns=ASSEMBLY_QC_COLUMNS)


__all__ = [
    "ASSEMBLY_QC_COLUMNS",
    "ASSEMBLY_QC_OUT_FILES",
    "AssemblyQualityResult",
    "load_assembly_quality_summary",
    "run_assembly_quality",
]
