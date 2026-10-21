from __future__ import annotations

import json
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from pipeline.CoreAlignment import ORDERED_ALIGNMENT_REL, file_sha256, read_ordering_summary
from pipeline.ToolEnv import configured_fastgear_settings

# fastGEAR reads the genome-ordered core alignment (the same file as Gubbins and the
# core tree); fastgear/input.json records what it ran on, so it reruns when that changes.
ALGORITHM_VERSION = 1
INPUT_STAMP = "input.json"

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


# The MATLAB R2016a runtime links libncurses.so.5, which the tools image (Ubuntu 24.04)
# lacks, so fastGEAR could not start. The launcher puts its working directory first on
# LD_LIBRARY_PATH, so a link there to the image's ncurses 6 lets it load. Images built
# from the current docker/Dockerfile carry the same link in the MCR folder.
_NCURSES_COMPAT = ("libncurses.so.5", "/usr/lib/x86_64-linux-gnu/libncursesw.so.6")


def _link_ncurses_compat(workdir: Path) -> None:
    link = workdir / _NCURSES_COMPAT[0]
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(_NCURSES_COMPAT[1])


def _image_identity(ctx) -> str | None:
    """Docker (or Podman) image id, or the Singularity image's path, size and mtime."""
    if ctx.runtime == "singularity" and ctx.sif_path:
        try:
            st = Path(ctx.sif_path).stat()
        except OSError:
            return None
        return f"{ctx.sif_path}:{st.st_size}:{st.st_mtime_ns}"
    from pipeline.ContainerRuntime import docker_image_info
    from pipeline.ToolEnv import oci_exe

    info = docker_image_info(ctx.image, oci_exe(ctx))
    return (info or {}).get("id") or None


def _fastgear_stamp(output_path: Path, settings: Mapping[str, Any], ctx) -> dict[str, Any]:
    alignment = output_path / ORDERED_ALIGNMENT_REL
    ordering = read_ordering_summary(output_path) or {}
    specs = Path(str(settings.get("specs_file") or ""))
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "input_path": str(alignment),
        "input_sha256": (ordering.get("output") or {}).get("sha256"),
        "settings": {
            key: settings.get(key) for key in ("launcher", "mcr_root", "specs_file", "extra_args")
        },
        # The default specs file lives in the image, whose identity covers it.
        "specs_sha256": file_sha256(specs) if specs.is_file() else None,
        "image": _image_identity(ctx),
    }


def fastgear_engine(ctx) -> str:
    """The container runtime fastGEAR runs under: the run's, or for a native run the
    docker or podman it falls back to."""
    from pipeline.ToolEnv import native_container_runtime

    return ctx.runtime if ctx.runtime != "native" else native_container_runtime()


def fastgear_reuse_issue(output_path: Path, fastgear_config: Mapping[str, Any] | None = None) -> str | None:
    """Why the existing fastGEAR output cannot be reused (None when it can)."""
    from pipeline.ToolEnv import current_container_context, make_container_context

    output_path = Path(output_path)
    try:
        recorded = json.loads((output_path / "fastgear" / INPUT_STAMP).read_text(encoding="utf-8"))
    except Exception:
        return "no input stamp"
    settings = _as_settings(fastgear_config)
    ctx = current_container_context() or make_container_context(output_dir=output_path)
    current = _fastgear_stamp(output_path, settings, ctx)
    if recorded.get("input_sha256") != current["input_sha256"] or not current["input_sha256"]:
        return "input alignment changed"
    for key in ("algorithm_version", "settings", "specs_sha256", "image"):
        if recorded.get(key) != current[key]:
            return f"{key.replace('_', ' ')} changed"
    return None


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

    alignment = output_path / ORDERED_ALIGNMENT_REL
    if not alignment.is_file():
        raise RuntimeError("The ordered core alignment (step core_alignment) is required before fastGEAR.")
    # Outputs of an earlier run (and the old staged input copy) must not survive a rerun.
    for name in ("input", "output", INPUT_STAMP, "fastgear_res.mat", "fastgear_res_snpData.mat", "sample_clades.tsv"):
        stale = fastgear_dir / name
        if stale.is_dir():
            shutil.rmtree(stale, ignore_errors=True)
        else:
            stale.unlink(missing_ok=True)
    result_mat_out = fastgear_dir / "fastgear_res.mat"
    _link_ncurses_compat(fastgear_dir)

    # fastGEAR (MATLAB MCR) now lives in the unified tools image at /opt. Run it with the
    # run's container runtime and identity-mounted host paths (native falls back to docker or podman,
    # since fastGEAR cannot run outside a container).
    ctx = current_container_context() or make_container_context(output_dir=output_path)
    cmd = container_image_cmd(
        ctx,
        [launcher, mcr_root, str(alignment), str(result_mat_out), specs_file, *extra_args],
        workdir=str(fastgear_dir),
        extra_env={"HOME": str(cache_dir), "MCR_CACHE_ROOT": str(cache_dir)},
        extra_binds=[fastgear_dir],
    )
    engine = fastgear_engine(ctx)

    logger.info("Running fastGEAR (%s).", engine)
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
            f"fastGEAR {engine} execution failed with exit code {result.returncode}"
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
    stamp_path = fastgear_dir / INPUT_STAMP
    stamp_path.write_text(json.dumps(_fastgear_stamp(output_path, settings, ctx), indent=2), encoding="utf-8")

    return FastGearResult(
        sample_clades_tsv=sample_clades_tsv,
        lineage_information=lineage_path,
        result_mat=result_mat,
    )


__all__ = [
    "FastGearResult",
    "fastgear_engine",
    "fastgear_reuse_issue",
    "parse_lineage_information",
    "run_fastgear_clades",
    "validate_sample_clades_tsv",
]
