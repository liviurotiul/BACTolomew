from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from pipeline.FastaUtils import prepare_query_fasta
from pipeline.RunSteps import ProgressCallback, report_progress
from pipeline.ToolEnv import configured_mobsuite_settings, tool_cmd


MOBSUITE_OUT_FILES = (
    "contigs.tsv",
    "plasmids.tsv",
    "summary.json",
)

MOBSUITE_CONTIG_COLUMNS = [
    "Sample",
    "Contig",
    "SourceContig",
    "PlasmidID",
    "IsPlasmid",
    "Prediction",
    "MoleculeType",
    "PrimaryClusterID",
    "SecondaryClusterID",
    "Size",
    "GC",
    "MD5",
    "CircularityStatus",
    "RepliconTypes",
    "RepliconAccessions",
    "RelaxaseTypes",
    "RelaxaseAccessions",
    "MPFType",
    "MPFAccessions",
    "OriTTypes",
    "OriTAccessions",
    "PredictedMobility",
    "MashNearestNeighbor",
    "MashNeighborDistance",
    "MashNeighborIdentification",
    "RepetitiveDNAID",
    "RepetitiveDNAType",
    "FilteringReason",
    "RawColumns",
]

MOBSUITE_PLASMID_COLUMNS = [
    "Sample",
    "PlasmidID",
    "NumContigs",
    "Size",
    "GC",
    "MD5",
    "RepliconTypes",
    "RepliconAccessions",
    "RelaxaseTypes",
    "RelaxaseAccessions",
    "MPFType",
    "MPFAccessions",
    "OriTTypes",
    "OriTAccessions",
    "PredictedMobility",
    "MashNearestNeighbor",
    "MashNeighborDistance",
    "MashNeighborIdentification",
    "PrimaryClusterID",
    "SecondaryClusterID",
    "PredictedHostRangeRank",
    "PredictedHostRangeName",
    "ObservedHostRangeRank",
    "ObservedHostRangeName",
    "ReportedHostRangeRank",
    "ReportedHostRangeName",
    "AssociatedPMIDs",
    "RawColumns",
]


@dataclass(frozen=True)
class MobSuiteResult:
    output_dir: Path
    contigs_tsv: Path
    plasmids_tsv: Path
    summary_json: Path
    sample_count: int
    plasmid_contig_count: int
    plasmid_count: int


def _normalize_contig(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    text = text.split()[0]
    match = re.match(r"^contig[_-]?(\d+)$", text, flags=re.IGNORECASE)
    if match:
        return str(int(match.group(1)))
    try:
        fval = float(text)
    except Exception:
        return text
    if fval.is_integer():
        return str(int(fval))
    return text


def _load_tsv(path: Path | None) -> pd.DataFrame:
    if path is None or not path.is_file():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def _find_file(output_dir: Path, name: str) -> Path | None:
    direct = output_dir / name
    if direct.is_file():
        return direct
    candidates = sorted(output_dir.rglob(name))
    return candidates[0] if candidates else None


def _find_mobtyper_file(output_dir: Path) -> Path | None:
    candidates = [
        path
        for path in output_dir.rglob("*")
        if path.is_file() and "mobtyper" in path.name.lower()
    ]
    return sorted(candidates, key=lambda path: (path.name != "mobtyper_results.txt", path.name))[0] if candidates else None


def _value(record: Mapping[str, Any], *names: str) -> str:
    lower_to_key = {str(key).lower().strip(): key for key in record.keys()}
    for name in names:
        key = lower_to_key.get(name.lower().strip())
        if key is not None:
            return str(record.get(key) or "").strip()
    return ""


def _is_plasmid_molecule(value: Any) -> bool:
    text = str(value or "").strip().lower()
    return text == "plasmid" or text.startswith("plasmid")


_MOBSUITE_EMPTY_VALUES = {
    "",
    "-",
    ".",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "<na>",
    "not tested",
}


def _has_mobsuite_value(value: Any) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    tokens = [token.strip() for token in re.split(r"[,;|]\s*", text) if token.strip()]
    if not tokens:
        tokens = [text]
    return any(token.casefold() not in _MOBSUITE_EMPTY_VALUES for token in tokens)


def _has_supported_mobility(value: Any) -> bool:
    text = str(value or "").strip().casefold()
    if not text or text in _MOBSUITE_EMPTY_VALUES:
        return False
    normalized = re.sub(r"[\s_]+", "-", text)
    if "non-mobil" in normalized or "nonmobil" in normalized:
        return False
    return normalized in {"mobilizable", "mobilisable", "conjugative"} or any(
        token in normalized for token in ("mobilizable", "mobilisable", "conjugative")
    )


def mobsuite_row_has_direct_plasmid_evidence(record: Mapping[str, Any]) -> bool:
    """Return true for direct plasmid evidence, not only MOB cluster/Mash similarity."""
    for names in (
        ("RepliconTypes", "rep_type(s)"),
        ("RepliconAccessions", "rep_type_accession(s)"),
        ("RelaxaseTypes", "relaxase_type(s)"),
        ("RelaxaseAccessions", "relaxase_type_accession(s)"),
        ("MPFType", "mpf_type"),
        ("MPFAccessions", "mpf_type_accession(s)"),
        ("OriTTypes", "orit_type(s)"),
        ("OriTAccessions", "orit_accession(s)"),
    ):
        if _has_mobsuite_value(_value(record, *names)):
            return True
    return _has_supported_mobility(_value(record, "PredictedMobility", "predicted_mobility"))


def mobsuite_row_is_supported_plasmid(record: Mapping[str, Any]) -> bool:
    raw_plasmid_call = (
        _is_plasmid_molecule(_value(record, "MoleculeType", "molecule_type"))
        or _value(record, "IsPlasmid", "is_plasmid").casefold() in {"true", "1", "yes", "y", "p", "plasmid"}
        or _value(record, "Prediction", "prediction").casefold() == "plasmid"
    )
    return raw_plasmid_call and mobsuite_row_has_direct_plasmid_evidence(record)


def _supported_plasmid_keys(
    plasmid_rows: list[dict[str, Any]],
) -> tuple[set[tuple[str, str]], set[tuple[str, str]], set[tuple[str, str]]]:
    plasmid_ids: set[tuple[str, str]] = set()
    primary_ids: set[tuple[str, str]] = set()
    secondary_ids: set[tuple[str, str]] = set()
    for row in plasmid_rows:
        if not mobsuite_row_has_direct_plasmid_evidence(row):
            continue
        sample = _value(row, "Sample", "sample")
        plasmid_id = _value(row, "PlasmidID", "sample_id")
        primary = _value(row, "PrimaryClusterID", "primary_cluster_id")
        secondary = _value(row, "SecondaryClusterID", "secondary_cluster_id")
        if plasmid_id:
            plasmid_ids.add((sample, plasmid_id))
        if primary:
            primary_ids.add((sample, primary))
        if secondary:
            secondary_ids.add((sample, secondary))
    return plasmid_ids, primary_ids, secondary_ids


def _apply_supported_plasmid_classification(
    contig_rows: list[dict[str, Any]],
    plasmid_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    plasmid_ids, primary_ids, secondary_ids = _supported_plasmid_keys(plasmid_rows)
    for row in contig_rows:
        raw_plasmid_call = (
            _is_plasmid_molecule(_value(row, "MoleculeType", "molecule_type"))
            or _value(row, "IsPlasmid", "is_plasmid").casefold() in {"true", "1", "yes", "y", "p", "plasmid"}
            or _value(row, "Prediction", "prediction").casefold() == "plasmid"
        )
        row_supported = mobsuite_row_has_direct_plasmid_evidence(row)
        sample = _value(row, "Sample", "sample_id")
        plasmid_supported = (
            (sample, _value(row, "PlasmidID", "sample_id")) in plasmid_ids
            or (sample, _value(row, "PrimaryClusterID", "primary_cluster_id")) in primary_ids
            or (sample, _value(row, "SecondaryClusterID", "secondary_cluster_id")) in secondary_ids
        )
        supported = raw_plasmid_call and (row_supported or plasmid_supported)
        row["IsPlasmid"] = "true" if supported else "false"
        if supported:
            row["Prediction"] = "plasmid"
            primary_cluster_id = _value(row, "PrimaryClusterID", "primary_cluster_id")
            sample = _value(row, "Sample", "sample_id")
            if not _value(row, "PlasmidID", "sample_id") and sample and primary_cluster_id:
                row["PlasmidID"] = f"{sample}:{primary_cluster_id}"
        elif raw_plasmid_call:
            row["Prediction"] = "plasmid_similarity_only"
            row["PlasmidID"] = ""
        else:
            row["Prediction"] = "chromosome"
            row["PlasmidID"] = ""
    return contig_rows


def _normalize_contig_rows(sample: str, raw_df: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if raw_df.empty:
        return rows
    for record in raw_df.to_dict(orient="records"):
        source_contig = _value(record, "contig_id", "contig", "sequence", "seq_id", "id", "name")
        contig = _normalize_contig(source_contig)
        if not contig:
            continue
        molecule_type = _value(record, "molecule_type")
        is_plasmid = _is_plasmid_molecule(molecule_type)
        primary_cluster_id = _value(record, "primary_cluster_id")
        plasmid_id = f"{sample}:{primary_cluster_id}" if is_plasmid and primary_cluster_id else ""
        rows.append(
            {
                "Sample": sample,
                "Contig": contig,
                "SourceContig": source_contig,
                "PlasmidID": plasmid_id,
                "IsPlasmid": "true" if is_plasmid else "false",
                "Prediction": "plasmid" if is_plasmid else "chromosome",
                "MoleculeType": molecule_type,
                "PrimaryClusterID": primary_cluster_id,
                "SecondaryClusterID": _value(record, "secondary_cluster_id"),
                "Size": _value(record, "size"),
                "GC": _value(record, "gc"),
                "MD5": _value(record, "md5"),
                "CircularityStatus": _value(record, "circularity_status"),
                "RepliconTypes": _value(record, "rep_type(s)"),
                "RepliconAccessions": _value(record, "rep_type_accession(s)"),
                "RelaxaseTypes": _value(record, "relaxase_type(s)"),
                "RelaxaseAccessions": _value(record, "relaxase_type_accession(s)"),
                "MPFType": _value(record, "mpf_type"),
                "MPFAccessions": _value(record, "mpf_type_accession(s)"),
                "OriTTypes": _value(record, "orit_type(s)"),
                "OriTAccessions": _value(record, "orit_accession(s)"),
                "PredictedMobility": _value(record, "predicted_mobility"),
                "MashNearestNeighbor": _value(record, "mash_nearest_neighbor"),
                "MashNeighborDistance": _value(record, "mash_neighbor_distance"),
                "MashNeighborIdentification": _value(record, "mash_neighbor_identification"),
                "RepetitiveDNAID": _value(record, "repetitive_dna_id"),
                "RepetitiveDNAType": _value(record, "repetitive_dna_type"),
                "FilteringReason": _value(record, "filtering_reason"),
                "RawColumns": json.dumps(record, sort_keys=True),
            }
        )
    return rows


def _normalize_plasmid_rows(sample: str, raw_df: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if raw_df.empty:
        return rows
    for index, record in enumerate(raw_df.to_dict(orient="records"), start=1):
        plasmid_id = _value(record, "sample_id") or f"{sample}_plasmid_{index}"
        rows.append(
            {
                "Sample": sample,
                "PlasmidID": plasmid_id,
                "NumContigs": _value(record, "num_contigs"),
                "Size": _value(record, "size"),
                "GC": _value(record, "gc"),
                "MD5": _value(record, "md5"),
                "RepliconTypes": _value(record, "rep_type(s)"),
                "RepliconAccessions": _value(record, "rep_type_accession(s)"),
                "RelaxaseTypes": _value(record, "relaxase_type(s)"),
                "RelaxaseAccessions": _value(record, "relaxase_type_accession(s)"),
                "MPFType": _value(record, "mpf_type"),
                "MPFAccessions": _value(record, "mpf_type_accession(s)"),
                "OriTTypes": _value(record, "orit_type(s)"),
                "OriTAccessions": _value(record, "orit_accession(s)"),
                "PredictedMobility": _value(record, "predicted_mobility"),
                "MashNearestNeighbor": _value(record, "mash_nearest_neighbor"),
                "MashNeighborDistance": _value(record, "mash_neighbor_distance"),
                "MashNeighborIdentification": _value(record, "mash_neighbor_identification"),
                "PrimaryClusterID": _value(record, "primary_cluster_id"),
                "SecondaryClusterID": _value(record, "secondary_cluster_id"),
                "PredictedHostRangeRank": _value(record, "predicted_host_range_overall_rank"),
                "PredictedHostRangeName": _value(record, "predicted_host_range_overall_name"),
                "ObservedHostRangeRank": _value(record, "observed_host_range_ncbi_rank"),
                "ObservedHostRangeName": _value(record, "observed_host_range_ncbi_name"),
                "ReportedHostRangeRank": _value(record, "reported_host_range_lit_rank"),
                "ReportedHostRangeName": _value(record, "reported_host_range_lit_name"),
                "AssociatedPMIDs": _value(record, "associated_pmid(s)"),
                "RawColumns": json.dumps(record, sort_keys=True),
            }
        )
    return rows


def _write_table(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        df = pd.DataFrame(columns=columns)
    df.to_csv(path, sep="\t", index=False)


def validate_mobsuite_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    settings = configured_mobsuite_settings(raw_config)
    db_dir = str(settings.get("database_directory") or "").strip()
    if db_dir:
        db_path = Path(db_dir).expanduser()
        db_path.mkdir(parents=True, exist_ok=True)
        settings["database_directory"] = str(db_path.resolve())
    return settings


def run_mobsuite(
    genome_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int,
    config: Mapping[str, Any] | None = None,
    progress: ProgressCallback = None,
) -> MobSuiteResult:
    settings = validate_mobsuite_settings(config)
    db_dir = str(settings.get("database_directory") or "").strip()
    extra_args = [str(arg) for arg in settings.get("extra_args") or [] if str(arg).strip()]
    mobsuite_root = output_path / "mobsuite"
    per_sample_root = mobsuite_root / "per_sample"
    mobsuite_root.mkdir(parents=True, exist_ok=True)
    per_sample_root.mkdir(parents=True, exist_ok=True)

    all_contig_rows: list[dict[str, Any]] = []
    all_plasmid_rows: list[dict[str, Any]] = []
    summary_samples: list[dict[str, Any]] = []
    failures: list[str] = []

    with tempfile.TemporaryDirectory(prefix="mobsuite_tmp_", dir=output_path) as tmp_text:
        tmp_root = Path(tmp_text)
        sanitized_root = tmp_root / "input"
        sanitized_root.mkdir(parents=True, exist_ok=True)

        genome_files = list(genome_files)
        for done, genome in enumerate(genome_files):
            report_progress(progress, done, len(genome_files))
            sample = genome.stem
            sample_tmp = tmp_root / f"{sample}__mobrecon_out"
            sample_out = per_sample_root / sample
            if sample_out.exists():
                shutil.rmtree(sample_out, ignore_errors=True)

            stdout_log = tmp_root / f"{sample}.command.stdout.log"
            stderr_log = tmp_root / f"{sample}.command.stderr.log"
            try:
                query_fa = prepare_query_fasta(genome, sanitized_root, logger, write_map=True)
                cmd = [
                    *tool_cmd("mob_recon"),
                    "--infile",
                    str(query_fa),
                    "--outdir",
                    str(sample_tmp),
                    "--num_threads",
                    str(max(1, int(threads))),
                    "--sample_id",
                    sample,
                    "--force",
                ]
                if db_dir:
                    cmd.extend(["--database_directory", db_dir])
                cmd.extend(extra_args)

                logger.info("Running MOB-suite MOB-recon on %s -> %s", genome.name, sample_out)
                logger.info("CMD: %s", " ".join(cmd))
                with stdout_log.open("w", encoding="utf-8") as out_fh, stderr_log.open("w", encoding="utf-8") as err_fh:
                    subprocess.run(cmd, check=True, stdout=out_fh, stderr=err_fh)

                contig_report = _find_file(sample_tmp, "contig_report.txt")
                if contig_report is None:
                    raise RuntimeError("MOB-recon did not produce contig_report.txt")
                mobtyper_report = _find_mobtyper_file(sample_tmp)

                contig_rows = _normalize_contig_rows(sample, _load_tsv(contig_report))
                raw_plasmid_rows = _normalize_plasmid_rows(sample, _load_tsv(mobtyper_report))
                contig_rows = _apply_supported_plasmid_classification(contig_rows, raw_plasmid_rows)
                plasmid_rows = [
                    row for row in raw_plasmid_rows if mobsuite_row_has_direct_plasmid_evidence(row)
                ]
                all_contig_rows.extend(contig_rows)
                all_plasmid_rows.extend(plasmid_rows)

                shutil.copytree(sample_tmp, sample_out, dirs_exist_ok=True)
                shutil.copy2(query_fa, sample_out / f"{sample}.fna")
                if stdout_log.is_file():
                    shutil.copy2(stdout_log, sample_out / "command.stdout.log")
                if stderr_log.is_file():
                    shutil.copy2(stderr_log, sample_out / "command.stderr.log")

                plasmid_contigs = sum(1 for row in contig_rows if str(row.get("IsPlasmid")).lower() == "true")
                chromosome_contigs = max(0, len(contig_rows) - plasmid_contigs)
                summary_samples.append(
                    {
                        "sample": sample,
                        "input_fasta": str(genome),
                        "sanitized_fasta": str(sample_out / f"{sample}.fna"),
                        "output_dir": str(sample_out),
                        "raw_contig_report": str(sample_out / contig_report.relative_to(sample_tmp)),
                        "raw_mobtyper_report": str(sample_out / mobtyper_report.relative_to(sample_tmp)) if mobtyper_report else "",
                        "plasmid_contigs": plasmid_contigs,
                        "chromosome_contigs": chromosome_contigs,
                        "plasmids": len(plasmid_rows),
                        "status": "success",
                    }
                )
            except subprocess.CalledProcessError as exc:
                failures.append(f"{sample}: mob_recon exited with code {exc.returncode}")
                logger.error("MOB-suite failed on %s: %s", genome.name, exc)
            except Exception as exc:
                failures.append(f"{sample}: {exc}")
                logger.error("MOB-suite failed on %s: %s", genome.name, exc)
        report_progress(progress, len(genome_files), len(genome_files))

    if failures:
        preview = "; ".join(failures[:5])
        if len(failures) > 5:
            preview += f"; ... ({len(failures)} failures total)"
        raise RuntimeError(f"MOB-suite failed for {len(failures)} genome(s): {preview}")

    contigs_tsv = mobsuite_root / "contigs.tsv"
    plasmids_tsv = mobsuite_root / "plasmids.tsv"
    summary_json = mobsuite_root / "summary.json"
    _write_table(contigs_tsv, all_contig_rows, MOBSUITE_CONTIG_COLUMNS)
    _write_table(plasmids_tsv, all_plasmid_rows, MOBSUITE_PLASMID_COLUMNS)
    summary_payload = {
        "database_directory": db_dir,
        "threads": int(threads),
        "sample_count": len(summary_samples),
        "plasmid_contig_count": sum(1 for row in all_contig_rows if str(row.get("IsPlasmid")).lower() == "true"),
        "plasmid_count": len(all_plasmid_rows),
        "samples": summary_samples,
    }
    summary_json.write_text(json.dumps(summary_payload, indent=2, sort_keys=True), encoding="utf-8")

    return MobSuiteResult(
        output_dir=mobsuite_root,
        contigs_tsv=contigs_tsv,
        plasmids_tsv=plasmids_tsv,
        summary_json=summary_json,
        sample_count=len(summary_samples),
        plasmid_contig_count=int(summary_payload["plasmid_contig_count"]),
        plasmid_count=len(all_plasmid_rows),
    )


def load_mobsuite_contig_predictions(output_path: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = Path(output_path) / "mobsuite" / "contigs.tsv"
    if not path.is_file():
        return {}
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}
    plasmids_path = Path(output_path) / "mobsuite" / "plasmids.tsv"
    try:
        plasmid_rows = (
            pd.read_csv(plasmids_path, sep="\t", dtype=str, keep_default_na=False).to_dict(orient="records")
            if plasmids_path.is_file()
            else []
        )
    except Exception:
        plasmid_rows = []
    contig_rows = _apply_supported_plasmid_classification(df.to_dict(orient="records"), plasmid_rows)
    out: dict[tuple[str, str], dict[str, str]] = {}
    for row in contig_rows:
        sample = str(row.get("Sample") or "").strip()
        contig = _normalize_contig(row.get("Contig"))
        if not sample or not contig:
            continue
        out[(sample, contig)] = {str(k): str(v) for k, v in row.items()}
    return out


__all__ = [
    "MOBSUITE_CONTIG_COLUMNS",
    "MOBSUITE_OUT_FILES",
    "MOBSUITE_PLASMID_COLUMNS",
    "MobSuiteResult",
    "load_mobsuite_contig_predictions",
    "mobsuite_row_has_direct_plasmid_evidence",
    "mobsuite_row_is_supported_plasmid",
    "run_mobsuite",
    "validate_mobsuite_settings",
]
