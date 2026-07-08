from __future__ import annotations

import gzip
import json
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from pipeline.ToolEnv import tool_cmd

PANAROO_ALIGNMENT_CANDIDATES: tuple[str, ...] = (
    "core_gene_alignment.aln",
    "core_gene_alignment.fa",
    "core_gene_alignment.fasta",
    "core_gene_alignment.aln.gz",
    "core_gene_alignment.fa.gz",
    "core_gene_alignment.fasta.gz",
)

GUBBINS_PREFIX = "core_gene"
GUBBINS_OUT_FILES: tuple[Path, ...] = (
    Path(f"{GUBBINS_PREFIX}.final_tree.tre"),
    Path(f"{GUBBINS_PREFIX}.filtered_polymorphic_sites.fasta"),
    Path(f"{GUBBINS_PREFIX}.recombination_predictions.gff"),
    Path("summary.json"),
)

GUBBINS_ALLOWED_BASES = set("ACGTN-")


@dataclass(frozen=True)
class GubbinsResult:
    output_dir: Path
    input_alignment: Path
    final_tree: Path
    filtered_alignment: Path
    recombination_gff: Path
    summary_json: Path


@dataclass(frozen=True)
class StagedAlignment:
    path: Path
    sequence_count: int
    alignment_length: int
    replaced_characters: int


def _resolve_panaroo_alignment(output_path: Path) -> Path:
    panaroo_out = Path(output_path) / "panaroo" / "out"
    if not panaroo_out.is_dir():
        raise RuntimeError(f"Panaroo output directory not found: {panaroo_out}")

    for name in PANAROO_ALIGNMENT_CANDIDATES:
        path = panaroo_out / name
        if path.is_file() and path.stat().st_size > 0:
            return path

    raise RuntimeError(
        f"Panaroo core alignment not found under {panaroo_out} "
        f"(looked for: {', '.join(PANAROO_ALIGNMENT_CANDIDATES)})."
    )


def _open_alignment(src: Path):
    if src.suffix == ".gz":
        return gzip.open(src, "rt", encoding="utf-8", errors="ignore")
    return src.open("r", encoding="utf-8", errors="ignore")


def _sanitize_alignment_line(line: str) -> tuple[str, int]:
    out: list[str] = []
    replaced = 0
    for char in line.strip():
        base = char.upper()
        if base in GUBBINS_ALLOWED_BASES:
            out.append(base)
        else:
            out.append("N")
            replaced += 1
    return "".join(out), replaced


def _stage_alignment(src: Path, gubbins_dir: Path, logger: logging.Logger) -> StagedAlignment:
    input_dir = gubbins_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    staged = input_dir / "core_gene_alignment.aln"

    sequence_count = 0
    replaced_characters = 0
    sequence_lengths: list[int] = []
    empty_sequences: list[str] = []
    current_name: str | None = None
    current_length = 0

    def finish_sequence() -> None:
        nonlocal current_name, current_length
        if current_name is None:
            return
        if current_length <= 0:
            empty_sequences.append(current_name)
        sequence_lengths.append(current_length)

    with _open_alignment(src) as fin, staged.open("w", encoding="utf-8") as fout:
        for raw_line in fin:
            if raw_line.startswith(">"):
                finish_sequence()
                header = raw_line[1:].strip()
                if not header:
                    raise RuntimeError(f"Invalid FASTA header in Panaroo core alignment: {src}")
                current_name = header
                current_length = 0
                sequence_count += 1
                fout.write(f">{header}\n")
                continue

            if current_name is None:
                if raw_line.strip():
                    raise RuntimeError(f"FASTA sequence data appears before the first header: {src}")
                continue

            sanitized, replaced = _sanitize_alignment_line(raw_line)
            if not sanitized:
                continue
            current_length += len(sanitized)
            replaced_characters += replaced
            fout.write(sanitized)
            fout.write("\n")
        finish_sequence()

    if sequence_count <= 0:
        raise RuntimeError(f"Gubbins input staging failed: no FASTA sequences found in {src}")
    if empty_sequences:
        preview = ", ".join(empty_sequences[:10])
        raise RuntimeError(f"Gubbins input staging found empty sequence(s): {preview}")
    distinct_lengths = sorted(set(sequence_lengths))
    if len(distinct_lengths) != 1:
        raise RuntimeError(
            "Gubbins input staging found non-uniform alignment lengths: "
            + ", ".join(str(length) for length in distinct_lengths[:10])
        )

    stats_path = input_dir / "core_gene_alignment.sanitization.json"
    _write_summary(
        stats_path,
        {
            "source_alignment": str(src),
            "staged_alignment": str(staged),
            "sequence_count": sequence_count,
            "alignment_length": distinct_lengths[0],
            "replaced_characters": replaced_characters,
            "allowed_bases": "".join(sorted(GUBBINS_ALLOWED_BASES)),
            "replacement": "unsupported alignment characters converted to N",
        },
    )
    logger.info(
        "Prepared Gubbins input alignment: %s -> %s (%d unsupported character(s) converted to N).",
        src.name,
        staged,
        replaced_characters,
    )
    return StagedAlignment(
        path=staged,
        sequence_count=sequence_count,
        alignment_length=distinct_lengths[0],
        replaced_characters=replaced_characters,
    )


def _format_subprocess_error(result: subprocess.CompletedProcess[str]) -> str:
    stderr = (result.stderr or "").strip()
    stdout = (result.stdout or "").strip()
    parts: list[str] = []
    if stdout:
        parts.append(stdout)
    if stderr and stderr != stdout:
        parts.append(stderr)
    if not parts:
        return f"exit code {result.returncode}"
    text = "\n".join(parts)
    lines = text.splitlines()
    return "\n".join(lines[-30:])


def _write_command_logs(gubbins_dir: Path, result: subprocess.CompletedProcess[str]) -> None:
    (gubbins_dir / "command.stdout.log").write_text(result.stdout or "", encoding="utf-8")
    (gubbins_dir / "command.stderr.log").write_text(result.stderr or "", encoding="utf-8")


def _write_summary(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def _ensure_outputs(paths: Sequence[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file() or path.stat().st_size <= 0]
    if missing:
        raise RuntimeError("Gubbins did not produce required output(s): " + ", ".join(missing))


def run_gubbins(
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int,
) -> GubbinsResult:
    output_path = Path(output_path).expanduser().resolve()
    gubbins_dir = output_path / "gubbins"
    gubbins_dir.mkdir(parents=True, exist_ok=True)

    panaroo_alignment = _resolve_panaroo_alignment(output_path)
    staged_alignment = _stage_alignment(panaroo_alignment, gubbins_dir, logger)

    cmd = [
        *tool_cmd("run_gubbins.py", workdir=str(gubbins_dir)),
        "--prefix",
        GUBBINS_PREFIX,
        "--threads",
        str(max(1, int(threads))),
        str(staged_alignment.path),
    ]
    logger.info("Running Gubbins on Panaroo core alignment.")
    result = subprocess.run(
        cmd,
        cwd=gubbins_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    _write_command_logs(gubbins_dir, result)
    if result.returncode != 0:
        raise RuntimeError(f"Gubbins failed: {_format_subprocess_error(result)}")

    final_tree = gubbins_dir / f"{GUBBINS_PREFIX}.final_tree.tre"
    filtered_alignment = gubbins_dir / f"{GUBBINS_PREFIX}.filtered_polymorphic_sites.fasta"
    recombination_gff = gubbins_dir / f"{GUBBINS_PREFIX}.recombination_predictions.gff"
    summary_json = gubbins_dir / "summary.json"
    _ensure_outputs((final_tree, filtered_alignment, recombination_gff))

    _write_summary(
        summary_json,
        {
            "tool": "gubbins",
            "input_alignment": str(panaroo_alignment),
            "staged_alignment": str(staged_alignment.path),
            "sequence_count": int(staged_alignment.sequence_count),
            "alignment_length": int(staged_alignment.alignment_length),
            "sanitized_characters": int(staged_alignment.replaced_characters),
            "threads": int(max(1, threads)),
            "outputs": {
                "final_tree": str(final_tree),
                "filtered_alignment": str(filtered_alignment),
                "recombination_gff": str(recombination_gff),
            },
        },
    )

    logger.info("Gubbins completed: %s", gubbins_dir)
    return GubbinsResult(
        output_dir=gubbins_dir,
        input_alignment=staged_alignment.path,
        final_tree=final_tree,
        filtered_alignment=filtered_alignment,
        recombination_gff=recombination_gff,
        summary_json=summary_json,
    )


__all__ = [
    "GUBBINS_OUT_FILES",
    "GubbinsResult",
    "run_gubbins",
]
