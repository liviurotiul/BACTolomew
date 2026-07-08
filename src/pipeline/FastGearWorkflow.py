from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Sequence

import pandas as pd

from pipeline.ToolEnv import configured_fastgear_settings

_PANAROO_ALIGNMENT_CANDIDATES = (
    "core_gene_alignment.aln",
    "core_gene_alignment.fa",
    "core_gene_alignment.fasta",
    "core_gene_alignment.aln.gz",
    "core_gene_alignment.fa.gz",
    "core_gene_alignment.fasta.gz",
)

_SAMPLE_CLADE_COLUMNS = [
    "Sample",
    "FastGearLineage",
    "FastGearCluster",
    "FastGearStrainIndex",
    "FastGearMethod",
]


@dataclass
class FastGearResult:
    sample_clades_tsv: Path
    lineage_information: Path
    result_mat: Path


def _as_settings(config: Mapping[str, Any] | None) -> dict[str, Any]:
    if config:
        return configured_fastgear_settings(config)
    return configured_fastgear_settings()


def _ensure_fastgear_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    settings = _as_settings(config)
    image = str(settings.get("image") or "").strip()
    if not image:
        raise RuntimeError(
            "fastGEAR Docker image is not configured. Set bactolomew_config.json fastgear.image."
        )
    return settings


def _resolve_panaroo_alignment(output_path: Path) -> Path:
    panaroo_out = output_path / "panaroo" / "out"
    if not panaroo_out.is_dir():
        raise RuntimeError(f"Panaroo output directory not found: {panaroo_out}")
    for name in _PANAROO_ALIGNMENT_CANDIDATES:
        path = panaroo_out / name
        if path.is_file():
            return path
    raise RuntimeError(
        f"Panaroo core alignment not found under {panaroo_out} "
        f"(looked for: {', '.join(_PANAROO_ALIGNMENT_CANDIDATES)})."
    )


def _stage_alignment(src: Path, fastgear_dir: Path, logger: logging.Logger) -> Path:
    input_dir = fastgear_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    staged = input_dir / "core_gene_alignment.aln"
    if src.suffix == ".gz":
        with gzip.open(src, "rt", encoding="utf-8", errors="ignore") as fin, staged.open(
            "w", encoding="utf-8"
        ) as fout:
            for line in fin:
                fout.write(line)
        logger.info(f"Prepared fastGEAR input alignment: {src.name} -> {staged}")
    else:
        shutil.copy2(src, staged)
        logger.info(f"Prepared fastGEAR input alignment: {src.name} -> {staged}")
    return staged


def _docker_client_available(docker_bin: str) -> bool:
    try:
        res = subprocess.run(
            [docker_bin, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return False
    return res.returncode == 0


def _docker_image_available(docker_bin: str, image: str) -> bool:
    try:
        res = subprocess.run(
            [docker_bin, "image", "inspect", image],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return False
    return res.returncode == 0


def _format_subprocess_error(result: subprocess.CompletedProcess[str]) -> str:
    stderr = (result.stderr or "").strip()
    stdout = (result.stdout or "").strip()
    parts: list[str] = []
    if stdout:
        parts.append(stdout)
    if stderr and stderr != stdout:
        parts.append(stderr)
    return "\n".join(parts[:2]).strip()


def _sample_clade_column_map(columns: Sequence[object]) -> dict[str, str]:
    normalized = {
        str(col): "".join(ch for ch in str(col).strip().lower() if ch.isalnum())
        for col in columns
    }
    wanted = {
        "sample": {"sample", "samples", "name", "strain", "isolate"},
        "lineage": {"lineage", "fastgearlineage"},
        "cluster": {"cluster", "fastgearcluster"},
        "strainindex": {"strainindex", "fastgearstrainindex", "index", "strainidx"},
    }

    mapped: dict[str, str] = {}
    for target, aliases in wanted.items():
        for original, key in normalized.items():
            if key in aliases:
                mapped[target] = original
                break
    return mapped


def _normalize_sample_clades(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
    if df.empty:
        raise RuntimeError(f"{context} is empty.")

    col_map = _sample_clade_column_map(df.columns)
    required = {"sample", "lineage", "cluster", "strainindex"}
    if not required.issubset(col_map):
        raise RuntimeError(
            f"{context} is missing required columns. Found columns: {', '.join(map(str, df.columns))}"
        )

    out = pd.DataFrame(
        {
            "Sample": df[col_map["sample"]].astype(str).str.strip(),
            "FastGearLineage": df[col_map["lineage"]],
            "FastGearCluster": df[col_map["cluster"]],
            "FastGearStrainIndex": df[col_map["strainindex"]],
        }
    )
    out = out[out["Sample"] != ""].copy()
    if out.empty:
        raise RuntimeError(f"{context} does not contain any sample names.")

    if out["Sample"].duplicated().any():
        dupes = sorted(set(out.loc[out["Sample"].duplicated(), "Sample"].tolist()))
        raise RuntimeError(f"{context} contains duplicate sample rows: {', '.join(dupes)}")

    sort_index = pd.to_numeric(out["FastGearStrainIndex"], errors="coerce")
    out["FastGearMethod"] = "fastgear"
    out = out[_SAMPLE_CLADE_COLUMNS]
    out = out.assign(_sort_index=sort_index).sort_values(
        ["_sort_index", "Sample"],
        kind="stable",
        na_position="last",
    )
    out = out.drop(columns=["_sort_index"]).reset_index(drop=True)
    return out


def _validate_expected_samples(df: pd.DataFrame, expected_samples: Sequence[str], *, context: str) -> None:
    expected = {str(sample).strip() for sample in expected_samples if str(sample).strip()}
    observed = {str(sample).strip() for sample in df["Sample"].tolist() if str(sample).strip()}
    if expected != observed:
        missing = sorted(expected - observed)
        unexpected = sorted(observed - expected)
        details: list[str] = []
        if missing:
            details.append(f"missing samples: {', '.join(missing)}")
        if unexpected:
            details.append(f"unexpected samples: {', '.join(unexpected)}")
        raise RuntimeError(f"{context} sample set mismatch ({'; '.join(details)})")


def parse_lineage_information(
    lineage_path: Path,
    *,
    expected_samples: Sequence[str],
) -> pd.DataFrame:
    try:
        df = pd.read_csv(
            lineage_path,
            sep=r"\s+",
            engine="python",
            comment="#",
        )
    except Exception as exc:
        raise RuntimeError(f"Failed to parse {lineage_path}: {exc}") from exc

    normalized = _normalize_sample_clades(df, context=str(lineage_path))
    _validate_expected_samples(normalized, expected_samples, context=str(lineage_path))
    return normalized


def validate_sample_clades_tsv(
    sample_clades_path: Path,
    *,
    expected_samples: Sequence[str],
) -> pd.DataFrame:
    try:
        df = pd.read_csv(sample_clades_path, sep="\t", dtype=str)
    except Exception as exc:
        raise RuntimeError(f"Failed to parse {sample_clades_path}: {exc}") from exc

    normalized = _normalize_sample_clades(df, context=str(sample_clades_path))
    _validate_expected_samples(normalized, expected_samples, context=str(sample_clades_path))
    return normalized


def run_fastgear_clades(
    output_path: Path,
    *,
    expected_samples: Sequence[str],
    logger: logging.Logger,
    fastgear_config: Mapping[str, Any] | None = None,
) -> FastGearResult:
    settings = _ensure_fastgear_config(fastgear_config)
    launcher = str(settings.get("launcher") or "/opt/run_fastGEAR.sh")
    mcr_root = str(settings.get("mcr_root") or "/opt/v901")
    specs_file = str(settings.get("specs_file") or "/opt/fG_input_specs.txt")
    extra_args = [str(token) for token in settings.get("extra_args") or []]

    from pipeline.ToolEnv import (
        container_image_cmd,
        current_container_context,
        make_container_context,
    )

    fastgear_dir = output_path / "fastgear"
    fastgear_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = fastgear_dir / "mcr_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    alignment_src = _resolve_panaroo_alignment(output_path)
    staged_alignment = _stage_alignment(alignment_src, fastgear_dir, logger)
    result_mat_out = fastgear_dir / "fastgear_res.mat"

    # fastGEAR (MATLAB MCR) now lives in the unified tools image at /opt. Run it with the
    # run's container runtime and identity-mounted host paths (native falls back to docker,
    # since fastGEAR cannot run outside a container).
    ctx = current_container_context() or make_container_context(output_dir=output_path)
    cmd = container_image_cmd(
        ctx,
        [launcher, mcr_root, str(staged_alignment), str(result_mat_out), specs_file, *extra_args],
        workdir=str(fastgear_dir),
        extra_env={"HOME": str(cache_dir), "MCR_CACHE_ROOT": str(cache_dir)},
        extra_binds=[fastgear_dir],
    )

    logger.info("Running fastGEAR (%s).", ctx.runtime if ctx.runtime != "native" else "docker")
    logger.info("CMD: " + " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        details = _format_subprocess_error(result)
        raise RuntimeError(
            f"fastGEAR Docker execution failed with exit code {result.returncode}"
            + (f": {details}" if details else ".")
        )

    result_mat = fastgear_dir / "fastgear_res.mat"
    lineage_path = fastgear_dir / "output" / "lineage_information.txt"
    if not result_mat.is_file():
        raise RuntimeError(f"fastGEAR finished but result file is missing: {result_mat}")
    if not lineage_path.is_file():
        raise RuntimeError(
            f"fastGEAR finished but lineage information is missing: {lineage_path}"
        )

    normalized = parse_lineage_information(lineage_path, expected_samples=expected_samples)
    sample_clades_tsv = fastgear_dir / "sample_clades.tsv"
    normalized.to_csv(sample_clades_tsv, sep="\t", index=False)
    logger.info(f"Normalized fastGEAR sample clades: {sample_clades_tsv}")

    return FastGearResult(
        sample_clades_tsv=sample_clades_tsv,
        lineage_information=lineage_path,
        result_mat=result_mat,
    )


__all__ = [
    "FastGearResult",
    "parse_lineage_information",
    "run_fastgear_clades",
    "validate_sample_clades_tsv",
]
