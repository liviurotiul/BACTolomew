from __future__ import annotations

import csv
import json
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from pipeline.ToolEnv import conda_tool_cmd

MLST_OUT_FILES: tuple[str, ...] = (
    "results.csv",
    "results.json",
)


def _candidate_cmds(*args: str) -> list[list[str]]:
    out: list[list[str]] = []
    env_cmd = conda_tool_cmd("mlst")
    if env_cmd:
        out.append([*env_cmd, *args])
    out.append(["mlst", *args])
    return out


def _resolve_mlst_base_cmd() -> list[str] | None:
    for probe_args in (("--version",), ("--help",)):
        for cmd in _candidate_cmds(*probe_args):
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=False)
            except FileNotFoundError:
                continue
            if res.returncode == 0:
                return cmd[:-len(probe_args)]
    return None


def _parse_mlst_stdout(stdout: str, genome_path: Path) -> dict[str, object]:
    line = ""
    for raw in stdout.splitlines():
        text = raw.strip()
        if text:
            line = text
            break
    if not line:
        raise RuntimeError(f"mlst returned no output for {genome_path.name}")

    fields = next(csv.reader([line], delimiter="\t"))
    if len(fields) < 3:
        fields = line.split()
    if len(fields) < 3:
        raise RuntimeError(f"mlst output could not be parsed for {genome_path.name}: {line}")

    return {
        "sample": genome_path.stem,
        "input_file": str(genome_path.resolve()),
        "scheme": fields[1],
        "sequence_type": fields[2],
        "alleles": fields[3:],
    }


def _run_mlst_one(base_cmd: list[str], genome_path: Path) -> dict[str, object]:
    cmd = [*base_cmd, str(genome_path)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(f"mlst executable not available for {genome_path.name}") from exc
    if res.returncode != 0:
        stderr = (res.stderr or "").strip()
        stdout = (res.stdout or "").strip()
        detail = stderr or stdout or f"exit code {res.returncode}"
        raise RuntimeError(f"mlst failed for {genome_path.name}: {detail}")
    return _parse_mlst_stdout(res.stdout or "", genome_path)


def run_mlst(
    genome_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int = 1,
) -> Path:
    genomes = [Path(p) for p in genome_files if Path(p).exists() and Path(p).is_file()]
    if not genomes:
        raise RuntimeError("No genome FASTA files provided to mlst.")

    base_cmd = _resolve_mlst_base_cmd()
    if base_cmd is None:
        raise RuntimeError("mlst executable not available (configured env or PATH).")

    workers = max(1, min(int(threads), len(genomes), 64))
    logger.info("Running mlst on %d genomes with %d worker(s).", len(genomes), workers)

    results: list[dict[str, object]] = []
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_map = {executor.submit(_run_mlst_one, base_cmd, genome): genome for genome in genomes}
        for future in as_completed(future_map):
            genome = future_map[future]
            try:
                results.append(future.result())
            except Exception as exc:
                logger.error("mlst failed for %s: %s", genome.name, exc)
                errors.append(str(exc))

    if errors:
        preview = "; ".join(errors[:5])
        if len(errors) > 5:
            preview += f"; ... ({len(errors)} failures total)"
        raise RuntimeError(preview)

    results.sort(key=lambda row: (str(row.get("sample") or ""), str(row.get("input_file") or "")))

    out_dir = output_path / "mlst"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "results.csv"
    json_path = out_dir / "results.json"

    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample", "input_file", "scheme", "sequence_type", "alleles"],
        )
        writer.writeheader()
        for row in results:
            writer.writerow(
                {
                    "sample": row["sample"],
                    "input_file": row["input_file"],
                    "scheme": row["scheme"],
                    "sequence_type": row["sequence_type"],
                    "alleles": ";".join(str(v) for v in (row.get("alleles") or [])),
                }
            )

    json_path.write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    logger.info("mlst completed: %s", csv_path)
    return out_dir


__all__ = [
    "MLST_OUT_FILES",
    "run_mlst",
]
