"""Population structure of the cohort (step ``fastgear_clades``): fastGEAR or FastBAPS.

The step finds the lineages Gubbins runs per (``phylo.gubbins_mode``) and the clades and
clusters the explorer colors the trees by. Two methods can produce them:

- ``fastgear``: fastGEAR (MATLAB runtime, x86_64 only) on the ordered core alignment,
  output in ``fastgear/`` (``pipeline.FastGearWorkflow``);
- ``fastbaps``: FastBAPS (R package, ``fastbaps`` conda env) on the same alignment,
  output in ``fastbaps/`` (``pipeline.FastBapsWorkflow``). Lineage = its Level 1 cluster,
  cluster = its Level 2 cluster (nested in Level 1).

Both write ``<folder>/sample_clades.tsv`` in one format (Sample, FastGearLineage,
FastGearCluster, FastGearStrainIndex, FastGearMethod; the column names stay for
compatibility, FastGearMethod says which method made the file), so Gubbins, the trees
and every panel read either unchanged.

``phylo.population_method`` (fastbaps | fastgear | auto, default fastbaps; a run can set
its own) picks the method. The tools image records the methods it holds (verify_image.sh:
``/opt/bactolomew/image_contents.json`` "population_methods", and the label
``org.bactolomew.population_methods``); an older image without that record has fastGEAR
when its fastGEAR files are present. auto takes FastBAPS when available, else fastGEAR,
else none (the step is skipped and Gubbins runs once on all genomes); a requested method
the image lacks stops the run before it starts. The run records the choice in
run_config.json ``inputs.population_structure`` ({requested, method, available}); a run
without that record (made before FastBAPS existed) used fastGEAR, and every reader
(``sample_clades_path``) finds its file in ``fastgear/`` as before.

No heavy imports: the explorer reads the record and the paths too.
"""
from __future__ import annotations

import json
import logging
import subprocess
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional

POPULATION_METHODS: tuple[str, ...] = ("fastgear", "fastbaps")
POPULATION_METHOD_CHOICES: tuple[str, ...] = ("auto", *POPULATION_METHODS)
# The order auto tries the methods in.
AUTO_PREFERENCE: tuple[str, ...] = ("fastbaps", "fastgear")
DEFAULT_METHOD = "fastbaps"
METHOD_LABELS = {"fastgear": "fastGEAR", "fastbaps": "FastBAPS"}
METHOD_DIRS = {"fastgear": "fastgear", "fastbaps": "fastbaps"}
SAMPLE_CLADES_NAME = "sample_clades.tsv"
# The method of a run without a population_structure record (all runs before FastBAPS).
LEGACY_METHOD = "fastgear"

IMAGE_CONTENTS_PATH = "/opt/bactolomew/image_contents.json"
METHODS_LABEL = "org.bactolomew.population_methods"
RUN_CONFIG_KEY = "population_structure"

_log = logging.getLogger(__name__)


class PopulationMethodUnavailable(RuntimeError):
    """phylo.population_method names a method the configured tools image lacks."""


def method_label(method: Optional[str]) -> str:
    return METHOD_LABELS.get(str(method or ""), str(method or "none"))


def method_dir(output_dir: Path, method: str) -> Path:
    return Path(output_dir) / METHOD_DIRS[method]


@dataclass(frozen=True)
class PopulationChoice:
    requested: str                        # auto | fastgear | fastbaps
    method: Optional[str]                 # the method the run uses; None: none available
    available: Optional[tuple[str, ...]]  # None: the image could not be checked

    @property
    def note(self) -> str:
        if self.method is not None:
            return ""
        return ("no population structure method (fastGEAR or FastBAPS) is available in the "
                "configured tools image")

    def record(self) -> dict[str, Any]:
        return {
            "requested": self.requested,
            "method": self.method,
            "available": list(self.available) if self.available is not None else None,
        }


def resolve_population_method(
    requested: str, available: Optional[tuple[str, ...]], *, explicit: bool = True
) -> PopulationChoice:
    """The method a run uses (see the module docstring). ``available`` None (the image
    could not be checked) keeps the behavior from before the check: auto is fastGEAR and
    a requested method is tried as asked. A method nobody set (``explicit`` False: the
    default) that the image lacks is replaced by auto's choice instead of stopping the run,
    so an image without FastBAPS (0.2.0 and older) keeps running fastGEAR."""
    requested = str(requested or DEFAULT_METHOD).strip().lower()
    if requested not in POPULATION_METHOD_CHOICES:
        raise ValueError(
            f"phylo.population_method: {requested!r} is not one of {', '.join(POPULATION_METHOD_CHOICES)}"
        )
    if available is None:
        return PopulationChoice(requested, LEGACY_METHOD if requested == "auto" else requested, None)
    if requested == "auto":
        method = next((m for m in AUTO_PREFERENCE if m in available), None)
        return PopulationChoice(requested, method, tuple(available))
    if requested not in available and not explicit:
        method = next((m for m in AUTO_PREFERENCE if m in available), None)
        return PopulationChoice(requested, method, tuple(available))
    if requested not in available:
        have = ", ".join(method_label(m) for m in available) or "neither fastGEAR nor FastBAPS"
        raise PopulationMethodUnavailable(
            f"phylo.population_method is {requested}, but the configured tools image has "
            f"{have}. Choose another population structure method, or configure an image that "
            f"includes {method_label(requested)}."
        )
    return PopulationChoice(requested, requested, tuple(available))


# ------------------------------------------------------------- what the image holds
_PROBE = (
    f"cat {IMAGE_CONTENTS_PATH} 2>/dev/null; echo; "
    "if [ -x /opt/run_fastGEAR.sh ] && [ -f /opt/fastGEAR ]; then echo @@fastgear-files; fi"
)
# Per image identity, the methods or None (the check failed: not retried for that image).
_cache: dict[str, Optional[tuple[str, ...]]] = {}
_MISSING = object()
_cache_lock = threading.Lock()


def parse_methods(text: str) -> tuple[str, ...]:
    """The known methods in a comma or space separated list, in POPULATION_METHODS order."""
    names = {part.strip().lower() for part in str(text or "").replace(",", " ").split()}
    return tuple(m for m in POPULATION_METHODS if m in names)


def _methods_from_probe(stdout: str) -> tuple[str, ...]:
    """The probe's output: image_contents.json (its population_methods when recorded),
    then a marker line when the fastGEAR files exist (images before the record)."""
    body, _, _ = stdout.partition("@@fastgear-files")
    try:
        contents = json.loads(body.strip() or "{}")
    except ValueError:
        contents = {}
    recorded = contents.get("population_methods") if isinstance(contents, dict) else None
    if isinstance(recorded, list):
        return parse_methods(" ".join(map(str, recorded)))
    return ("fastgear",) if "@@fastgear-files" in stdout else ()


def _native_fastbaps() -> bool:
    from pipeline.ToolEnv import get_tool_env_prefix

    return bool(get_tool_env_prefix("fastbaps"))


def image_population_methods(ctx=None) -> Optional[tuple[str, ...]]:
    """The population methods the run's tools image holds, cached per image identity
    (None when the image cannot be checked: no engine, image missing, probe failed; a
    failed probe is not repeated for the same image).
    A native context runs fastGEAR from the tools image and FastBAPS from a local
    ``fastbaps`` env (config envs.fastbaps or envs/fastbaps)."""
    from pipeline.ContainerRuntime import docker_image_info
    from pipeline.FastGearWorkflow import _image_identity
    from pipeline.ToolEnv import (
        container_image_cmd,
        current_container_context,
        make_container_context,
        oci_exe,
    )

    ctx = ctx or current_container_context() or make_container_context()
    try:
        identity = _image_identity(ctx)
    except Exception:
        identity = None
    methods: Optional[tuple[str, ...]]
    if identity is None:
        methods = None
    else:
        with _cache_lock:
            cached = _cache.get(identity, _MISSING)
        if cached is not _MISSING:
            methods = cached
        else:
            methods = None
            if ctx.runtime != "singularity":
                info = docker_image_info(ctx.image, oci_exe(ctx)) or {}
                label = (info.get("labels") or {}).get(METHODS_LABEL)
                if label is not None:
                    methods = parse_methods(label)
            if methods is None:
                cmd = container_image_cmd(ctx, ["sh", "-c", _PROBE])
                try:
                    res = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=False)
                except Exception as exc:
                    _log.warning("Could not check the population methods of the tools image: %s", exc)
                    res = None
                if res is not None and res.returncode == 0:
                    methods = _methods_from_probe(res.stdout)
                elif res is not None:
                    _log.warning("Could not check the population methods of the tools image (exit %s).",
                                 res.returncode)
            with _cache_lock:
                _cache[identity] = methods
    if ctx.runtime == "native":
        image_fastgear = methods is None or "fastgear" in methods
        return tuple(m for m in POPULATION_METHODS
                     if (m == "fastgear" and image_fastgear) or (m == "fastbaps" and _native_fastbaps()))
    return methods


def population_method_explicit(run_phylo: Optional[Mapping[str, Any]] = None) -> bool:
    """True when phylo.population_method is set for the run (DataInput.phylo_config: the
    New Run choice, --population-method, a resume) or in the config, not just the default."""
    from pipeline.ToolEnv import _load_app_config

    config_phylo = _load_app_config().get("phylo")
    return any(isinstance(src, Mapping) and src.get("population_method")
               for src in (run_phylo, config_phylo))


def configured_population_choice(
    phylo_settings: Mapping[str, Any], ctx=None, *, explicit: bool = True
) -> PopulationChoice:
    """The method of a run with these phylo settings on the run's tools image (raises
    PopulationMethodUnavailable for an explicitly requested method the image lacks)."""
    requested = str(phylo_settings.get("population_method") or DEFAULT_METHOD)
    available = image_population_methods(ctx)
    return resolve_population_method(requested, available, explicit=explicit)


# ------------------------------------------------------------- the run's record
def run_population_record(output_dir: Path) -> Optional[dict[str, Any]]:
    """run_config.json inputs.population_structure (None for runs before it)."""
    try:
        cfg = json.loads((Path(output_dir) / "run_config.json").read_text(encoding="utf-8"))
    except Exception:
        return None
    inputs = cfg.get("inputs") if isinstance(cfg, dict) else None
    record = inputs.get(RUN_CONFIG_KEY) if isinstance(inputs, dict) else None
    return record if isinstance(record, dict) else None


def run_population_method(output_dir: Path) -> Optional[str]:
    """The method the run used: its record, fastGEAR for a run before the record, None
    when the run had none."""
    record = run_population_record(output_dir)
    if record is None:
        return LEGACY_METHOD
    method = record.get("method")
    return method if method in POPULATION_METHODS else None


def resume_phylo_overrides(output_dir: Path, requested: Optional[str] = None) -> dict[str, str]:
    """The phylo settings a resume passes (DataInput.phylo_config): ``requested`` when
    given, else the method the run used, so a run keeps its population method across
    resumes: fastGEAR for a run from before the record. Nothing for a run that had no
    method, or a folder without run_config.json (the config's setting applies)."""
    if requested:
        return {"population_method": str(requested)}
    if not (Path(output_dir) / "run_config.json").is_file():
        return {}
    method = run_population_method(output_dir)
    return {"population_method": method} if method in POPULATION_METHODS else {}


def run_population_label(output_dir: Path) -> str:
    """The name the explorer gives the run's clades and clusters: fastGEAR or FastBAPS
    (fastGEAR for a run without a method, whose clade legends stay empty)."""
    return method_label(run_population_method(output_dir) or LEGACY_METHOD)


def sample_clades_path(output_dir: Path, method: Optional[str] = "") -> Optional[Path]:
    """<method folder>/sample_clades.tsv of the run's method (or of ``method``); None
    when the run had no method. The file may not exist yet."""
    method = run_population_method(output_dir) if method == "" else method
    if method not in POPULATION_METHODS:
        return None
    return method_dir(output_dir, method) / SAMPLE_CLADES_NAME


__all__ = [
    "AUTO_PREFERENCE",
    "DEFAULT_METHOD",
    "LEGACY_METHOD",
    "METHOD_DIRS",
    "METHOD_LABELS",
    "METHODS_LABEL",
    "POPULATION_METHODS",
    "POPULATION_METHOD_CHOICES",
    "PopulationChoice",
    "PopulationMethodUnavailable",
    "configured_population_choice",
    "image_population_methods",
    "method_dir",
    "method_label",
    "population_method_explicit",
    "parse_methods",
    "resolve_population_method",
    "resume_phylo_overrides",
    "run_population_label",
    "run_population_method",
    "run_population_record",
    "sample_clades_path",
]
