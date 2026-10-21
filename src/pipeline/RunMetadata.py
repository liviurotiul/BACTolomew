"""The run's own sample table, for a run opened without a metadata CSV.

One row per genome with what the run itself knows about it: the species (from its
MLST scheme), the MLST sequence type, the fastGEAR clade and cluster and the assembly
QC summary. Every value is text, as if read from a CSV; the explorer types the
columns itself. No Qt.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

import pandas as pd

from pipeline.AmrFinderWorkflow import MLST_SCHEME_ORGANISMS
from pipeline.AssemblyQualityWorkflow import load_assembly_quality_summary
from pipeline.ToolEnv import project_root

# Columns other than the assembly QC ones (MLSTScheme only when the genomes were typed
# with different schemes); ST, Clade and Cluster are labels, not numbers.
RUN_METADATA_COLUMNS = ("Sample", "Species", "MLSTScheme", "ST", "Clade", "Cluster")
RUN_METADATA_LABEL_COLUMNS = ("ST", "Clade", "Cluster")
SCHEME_CATALOG_NAME = "scheme_species_catalog.json"
_MISSING = {"", "-", "nan", "none", "na", "n/a"}


def _clean(value: Any) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text.lower() in _MISSING else text


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _mlst_rows(output_dir: Path) -> list[Dict[str, str]]:
    """(sample, scheme, sequence_type) per genome from mlst/results.json, else results.csv."""
    raw = _read_json(output_dir / "mlst" / "results.json")
    if not isinstance(raw, list):
        try:
            with (output_dir / "mlst" / "results.csv").open(newline="", encoding="utf-8") as handle:
                raw = list(csv.DictReader(handle))
        except OSError:
            raw = []
    return [
        {
            "sample": _clean(row.get("sample")),
            "scheme": _clean(row.get("scheme")),
            "sequence_type": _clean(row.get("sequence_type")),
        }
        for row in raw
        if isinstance(row, Mapping) and _clean(row.get("sample"))
    ]


def _scheme_catalog_path(output_dir: Path) -> Path:
    config = _read_json(output_dir / "run_config.json")
    inputs = config.get("inputs") if isinstance(config, Mapping) else None
    chewbbaca = inputs.get("chewbbaca") if isinstance(inputs, Mapping) else None
    configured = str((chewbbaca or {}).get("catalog_path") or "").strip() if isinstance(chewbbaca, Mapping) else ""
    if configured and Path(configured).expanduser().is_file():
        return Path(configured).expanduser()
    return project_root() / "DB" / "cgmlst" / SCHEME_CATALOG_NAME


def mlst_scheme_species(output_dir: Path) -> Dict[str, str]:
    """MLST scheme (lower case) -> species name: the cgMLST scheme catalog first, then
    the schemes AMRFinderPlus knows ('Acinetobacter_baumannii' -> 'Acinetobacter baumannii')."""
    out = {scheme: organism.replace("_", " ") for scheme, organism in MLST_SCHEME_ORGANISMS.items()}
    raw = _read_json(_scheme_catalog_path(output_dir))
    entries = raw.get("schemes") if isinstance(raw, Mapping) else raw
    for entry in entries if isinstance(entries, list) else []:
        if not isinstance(entry, Mapping):
            continue
        species = [_clean(name) for name in entry.get("species") or [] if _clean(name)]
        if not species:
            continue
        for scheme in entry.get("mlst_schemes") or []:
            if _clean(scheme):
                out[_clean(scheme).lower()] = species[0]
    return out


def _fastgear_rows(output_dir: Path) -> Dict[str, tuple[str, str]]:
    """sample -> (clade, cluster) from fastgear/sample_clades.tsv."""
    out: Dict[str, tuple[str, str]] = {}
    try:
        with (output_dir / "fastgear" / "sample_clades.tsv").open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = _clean(row.get("Sample"))
                if sample and sample not in out:
                    out[sample] = (_clean(row.get("FastGearLineage")), _clean(row.get("FastGearCluster")))
    except OSError:
        pass
    return out


def build_run_metadata(output_dir: Path, samples: Iterable[str] = ()) -> pd.DataFrame:
    """The run's sample table (``RUN_METADATA_COLUMNS`` plus the assembly QC columns).

    Rows are ``samples`` when given (the genomes the explorer shows), else every genome
    any of the run's outputs names. Columns the run has no value for are left out, so
    a run without fastGEAR has no Clade column. Missing values are empty strings.
    """
    output_dir = Path(output_dir)
    mlst = {row["sample"]: row for row in _mlst_rows(output_dir)}
    species_by_scheme = mlst_scheme_species(output_dir) if mlst else {}
    fastgear = _fastgear_rows(output_dir)
    qc = load_assembly_quality_summary(output_dir)
    qc = qc.astype(str) if not qc.empty and "Sample" in qc.columns else pd.DataFrame()
    if not qc.empty:
        qc["Sample"] = qc["Sample"].map(_clean)
        qc = qc[qc["Sample"] != ""].drop_duplicates(subset=["Sample"], keep="first")

    wanted = [_clean(sample) for sample in samples]
    names = sorted({name for name in wanted if name}) if wanted else sorted(
        set(mlst) | set(fastgear) | (set(qc["Sample"]) if not qc.empty else set())
    )
    columns: Dict[str, list[str]] = {
        "Sample": names,
        "Species": [species_by_scheme.get(mlst.get(name, {}).get("scheme", "").lower(), "") for name in names],
        "MLSTScheme": [mlst.get(name, {}).get("scheme", "") for name in names],
        "ST": [mlst.get(name, {}).get("sequence_type", "") for name in names],
        "Clade": [fastgear.get(name, ("", ""))[0] for name in names],
        "Cluster": [fastgear.get(name, ("", ""))[1] for name in names],
    }
    if len({scheme for scheme in columns["MLSTScheme"] if scheme}) < 2:
        del columns["MLSTScheme"]
    df = pd.DataFrame({name: values for name, values in columns.items() if name == "Sample" or any(values)})
    if not qc.empty:
        df = df.merge(qc, on="Sample", how="left").fillna("")
    return df.astype(str)


__all__ = [
    "RUN_METADATA_COLUMNS",
    "RUN_METADATA_LABEL_COLUMNS",
    "build_run_metadata",
    "mlst_scheme_species",
]
