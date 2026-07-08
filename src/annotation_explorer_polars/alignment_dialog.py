from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from PyQt6 import QtCore, QtWidgets


class AlignmentDialog(QtWidgets.QDialog):
    def __init__(self, parent: Optional[QtWidgets.QWidget], ann_df: pd.DataFrame) -> None:
        super().__init__(parent)
        self.ann_df = ann_df
        self.reference_path: Optional[str] = None
        self.alignment_path: Optional[Path] = None
        self.mutations_csv_path: Optional[Path] = None
        self.current_tmpdir: Optional[Path] = None
        self.seq_type_group: Optional[QtWidgets.QButtonGroup] = None
        self.setWindowTitle("Align annotations")
        self.resize(620, 420)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        gene_row = QtWidgets.QHBoxLayout()
        gene_label = QtWidgets.QLabel("GeneID")
        self.gene_combo = QtWidgets.QComboBox()
        self.gene_combo.setEditable(True)
        self.gene_combo.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        gene_ids = sorted(
            {
                str(v).strip()
                for v in self.ann_df.get("GeneID", pd.Series()).dropna().unique().tolist()
                if str(v).strip() and str(v).lower() != "nan"
            }
        )
        self.gene_combo.addItems(gene_ids)
        completer = self.gene_combo.completer()
        if completer:
            completer.setFilterMode(QtCore.Qt.MatchFlag.MatchContains)
            completer.setCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
        gene_row.addWidget(gene_label)
        gene_row.addWidget(self.gene_combo, stretch=1)
        layout.addLayout(gene_row)

        seq_row = QtWidgets.QHBoxLayout()
        seq_label = QtWidgets.QLabel("Sequence type")
        self.nt_radio = QtWidgets.QRadioButton("Nucleotide")
        self.aa_radio = QtWidgets.QRadioButton("Amino acid")
        self.nt_radio.setChecked(True)
        self.seq_type_group = QtWidgets.QButtonGroup(self)
        self.seq_type_group.addButton(self.nt_radio)
        self.seq_type_group.addButton(self.aa_radio)
        seq_row.addWidget(seq_label)
        seq_row.addWidget(self.nt_radio)
        seq_row.addWidget(self.aa_radio)
        seq_row.addStretch(1)
        layout.addLayout(seq_row)

        ref_row = QtWidgets.QHBoxLayout()
        ref_label = QtWidgets.QLabel("Reference FASTA (optional)")
        self.ref_path_edit = QtWidgets.QLineEdit()
        self.ref_path_edit.setReadOnly(True)
        ref_btn = QtWidgets.QPushButton("Browse…")
        ref_btn.clicked.connect(self._choose_reference)
        ref_row.addWidget(ref_label)
        ref_row.addWidget(self.ref_path_edit, stretch=1)
        ref_row.addWidget(ref_btn)
        layout.addLayout(ref_row)

        paste_label = QtWidgets.QLabel("Or paste reference sequence (optional)")
        self.ref_text_edit = QtWidgets.QPlainTextEdit()
        self.ref_text_edit.setPlaceholderText("Paste reference FASTA entry or raw sequence")
        self.ref_text_edit.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        self.ref_text_edit.setFixedHeight(110)
        layout.addWidget(paste_label)
        layout.addWidget(self.ref_text_edit)

        self.status_label = QtWidgets.QLabel(
            "Select a GeneID, choose nucleotide or amino acid, and add an optional reference (file or pasted)."
        )
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addStretch(1)
        self.align_btn = QtWidgets.QPushButton("Align")
        self.align_btn.clicked.connect(self._align)
        self.download_btn = QtWidgets.QPushButton("Download result")
        self.download_btn.setEnabled(False)
        self.download_btn.clicked.connect(self._download_alignment)
        self.download_mut_btn = QtWidgets.QPushButton("Download mutations (CSV)")
        self.download_mut_btn.setEnabled(False)
        self.download_mut_btn.clicked.connect(self._download_mutations)
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(self.align_btn)
        btn_row.addWidget(self.download_btn)
        btn_row.addWidget(self.download_mut_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    def _choose_reference(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select reference FASTA",
            "",
            "FASTA files (*.fasta *.fa *.fna *.fa.gz *.fasta.gz);;All files (*)",
        )
        if path:
            self.reference_path = path
            self.ref_path_edit.setText(path)

    def _sequence_column(self) -> str:
        return "Amino_Acid_Sequence" if self.aa_radio.isChecked() else "Nucleotide_Sequence"

    def _sequence_label(self) -> str:
        return "amino acid" if self.aa_radio.isChecked() else "nucleotide"

    def _reference_text(self) -> str:
        return self.ref_text_edit.toPlainText().strip()

    def _read_reference_file(self) -> str:
        if not self.reference_path:
            return ""
        try:
            return Path(self.reference_path).read_text()
        except Exception:
            return ""

    def _collect_sequences(self, gene_id: str, seq_col: str) -> List[Tuple[str, str]]:
        rows = self.ann_df[
            self.ann_df.get("GeneID", pd.Series()).fillna("").astype(str).str.strip() == gene_id
        ]
        records: List[Tuple[str, str]] = []
        for _, row in rows.iterrows():
            seq = (row.get(seq_col) or "").strip()
            if not seq:
                continue
            header_parts = [
                str(row.get("AnnotationID") or "").strip(),
                str(row.get("Sample") or "").strip(),
                str(row.get("Contig") or "").strip(),
                str(row.get("Gene") or "").strip(),
            ]
            header = "|".join([p for p in header_parts if p]) or f"{gene_id}_{len(records) + 1}"
            records.append((header, seq))
        return records

    def _collect_sample_sequences(self, gene_id: str, seq_col: str) -> Dict[str, List[Tuple[str, str]]]:
        rows = self.ann_df[
            self.ann_df.get("GeneID", pd.Series()).fillna("").astype(str).str.strip() == gene_id
        ]
        sample_records: Dict[str, List[Tuple[str, str]]] = {}
        for _, row in rows.iterrows():
            seq = (row.get(seq_col) or "").strip()
            if not seq:
                continue
            sample = str(row.get("Sample") or "").strip() or "sample"
            header_parts = [
                str(row.get("AnnotationID") or "").strip(),
                str(row.get("Sample") or "").strip(),
                str(row.get("Contig") or "").strip(),
                str(row.get("Gene") or "").strip(),
            ]
            header = "|".join([p for p in header_parts if p]) or f"{gene_id}_{len(sample_records.get(sample, [])) + 1}"
            sample_records.setdefault(sample, []).append((header, seq))
        return sample_records

    def _prepare_reference_file(self, manual_reference: str, tmpdir: Path) -> Optional[Path]:
        """
        Choose a reference FASTA for snippy: prefer the selected file, fall back to pasted text.
        Returns a path to a copy under tmpdir, or None if nothing is available.
        """
        if self.reference_path:
            src = Path(self.reference_path)
            try:
                ref_dest = tmpdir / "reference.fasta"
                ref_dest.write_bytes(src.read_bytes())
                return ref_dest
            except Exception:
                pass
        manual = (manual_reference or "").strip()
        if manual:
            ref_dest = tmpdir / "reference_manual.fasta"
            if not manual.lstrip().startswith(">"):
                manual = f">reference\n{manual}"
            if not manual.endswith("\n"):
                manual += "\n"
            ref_dest.write_text(manual)
            return ref_dest
        return None

    def _write_sample_fastas(
        self, sample_records: Dict[str, List[Tuple[str, str]]], tmpdir: Path
    ) -> Dict[str, Path]:
        """
        Write per-sample FASTA files for snippy and return a map of sample -> fasta path.
        """
        out: Dict[str, Path] = {}
        used_names = set()
        input_dir = tmpdir / "snippy_inputs"
        input_dir.mkdir(parents=True, exist_ok=True)
        for sample, recs in sample_records.items():
            if not recs:
                continue
            safe = "".join(ch if (ch.isalnum() or ch in "-_.") else "_" for ch in sample) or "sample"
            name = safe
            i = 1
            while name in used_names:
                i += 1
                name = f"{safe}_{i}"
            used_names.add(name)
            path = input_dir / f"{name}.fasta"
            with path.open("w") as fh:
                for header, seq in recs:
                    fh.write(f">{header}\n{seq}\n")
            out[sample] = path
        return out

    def _parse_vcf_mutations(self, vcf_path: Path) -> str:
        if not vcf_path.exists():
            return ""
        mutations: List[str] = []
        try:
            with vcf_path.open() as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 5:
                        continue
                    chrom, pos, _id, ref, alt = parts[:5]
                    alts = [a for a in alt.split(",") if a and a != "."]
                    for a in alts:
                        mutations.append(f"{chrom}:{pos}{ref}>{a}")
        except Exception:
            return ""
        return "; ".join(mutations)

    def _run_snippy(
        self,
        gene_id: str,
        sample_records: Dict[str, List[Tuple[str, str]]],
        reference_fasta: Optional[Path],
        tmpdir: Path,
    ) -> Tuple[str, Optional[Path]]:
        from pipeline.ToolEnv import gui_tool_cmd, runtime_is_container

        if reference_fasta is None or not reference_fasta.exists():
            return ("Snippy skipped: reference unavailable.", None)
        if not runtime_is_container() and shutil.which("snippy") is None:
            return ("Snippy skipped: snippy not found in PATH.", None)
        sample_fastas = self._write_sample_fastas(sample_records, tmpdir)
        if not sample_fastas:
            return ("Snippy skipped: no sequences available for snippy.", None)

        snippy_root = tmpdir / "snippy_outputs"
        snippy_root.mkdir(parents=True, exist_ok=True)
        mutations_rows: List[Dict[str, str]] = []
        for sample, fasta_path in sample_fastas.items():
            sample_dir = snippy_root / (fasta_path.stem + "_snippy")
            sample_dir.mkdir(parents=True, exist_ok=True)
            cmd = [
                *gui_tool_cmd("snippy", "snippy", mount_paths=[str(tmpdir)]),
                "--cpus",
                "1",
                "--outdir",
                str(sample_dir),
                "--ref",
                str(reference_fasta),
                "--ctgs",
                str(fasta_path),
            ]
            try:
                res = subprocess.run(cmd, capture_output=True, text=True)
                if res.returncode != 0:
                    msg = res.stderr.strip() or res.stdout.strip() or "snippy failed"
                    mutations_rows.append({"Sample": sample, "locus": gene_id, "mutation list": f"error: {msg}"})
                    continue
            except Exception as exc:
                mutations_rows.append({"Sample": sample, "locus": gene_id, "mutation list": f"error: {exc}"})
                continue

            vcf_path = sample_dir / "snps.vcf"
            mut_list = self._parse_vcf_mutations(vcf_path)
            mutations_rows.append({"Sample": sample, "locus": gene_id, "mutation list": mut_list})

        if not mutations_rows:
            return ("Snippy finished but found no mutation records.", None)
        csv_path = tmpdir / "snippy_mutations.csv"
        pd.DataFrame(mutations_rows).to_csv(csv_path, index=False)
        return (f"Snippy completed for {len(mutations_rows)} sample(s).", csv_path)

    def _write_reference_entry(self, fh: Any, text: str, header: str, has_existing: bool) -> bool:
        cleaned = (text or "").strip()
        if not cleaned:
            return has_existing
        if not cleaned.lstrip().startswith(">"):
            cleaned = f">{header}\n{cleaned}"
        if not cleaned.endswith("\n"):
            cleaned += "\n"
        if has_existing:
            fh.write("\n")
        fh.write(cleaned)
        return True

    def _write_input_fasta(self, records: List[Tuple[str, str]], path: Path, manual_reference: str) -> None:
        with path.open("w") as fh:
            wrote_any = False
            for header, seq in records:
                fh.write(f">{header}\n{seq}\n")
                wrote_any = True
            ref_file_text = self._read_reference_file()
            wrote_any = self._write_reference_entry(fh, ref_file_text, "reference_file", wrote_any)
            self._write_reference_entry(fh, manual_reference, "reference_manual", wrote_any)

    def _find_aligner(self) -> Optional[str]:
        from pipeline.ToolEnv import runtime_is_container

        if runtime_is_container():
            return "mafft"  # provided by the tools image
        for cmd in ["mafft", "muscle", "clustalo"]:
            if shutil.which(cmd):
                return cmd
        return None

    def _perform_alignment(self, input_path: Path, output_path: Path) -> str:
        aligner = self._find_aligner()
        if not aligner:
            shutil.copyfile(input_path, output_path)
            return "No aligner (mafft/muscle/clustalo) found. Saved unaligned FASTA."
        try:
            if aligner == "mafft":
                from pipeline.ToolEnv import gui_tool_cmd

                res = subprocess.run(
                    [*gui_tool_cmd("mafft", "mafft", mount_paths=[str(input_path.parent)]), "--auto", str(input_path)],
                    capture_output=True,
                    text=True,
                )
                if res.returncode == 0 and res.stdout.strip():
                    output_path.write_text(res.stdout)
                    return "Alignment completed with mafft."
                msg = res.stderr.strip() or res.stdout.strip() or "Unknown error"
                shutil.copyfile(input_path, output_path)
                return f"mafft failed: {msg}"
            if aligner == "muscle":
                res = subprocess.run(
                    [aligner, "-in", str(input_path), "-out", str(output_path)],
                    capture_output=True,
                    text=True,
                )
                if res.returncode == 0:
                    return "Alignment completed with muscle."
                msg = res.stderr.strip() or res.stdout.strip() or "Unknown error"
                shutil.copyfile(input_path, output_path)
                return f"muscle failed: {msg}"
            if aligner == "clustalo":
                res = subprocess.run(
                    [aligner, "-i", str(input_path), "-o", str(output_path), "--force"],
                    capture_output=True,
                    text=True,
                )
                if res.returncode == 0:
                    return "Alignment completed with clustalo."
                msg = res.stderr.strip() or res.stdout.strip() or "Unknown error"
                shutil.copyfile(input_path, output_path)
                return f"clustalo failed: {msg}"
        except FileNotFoundError:
            shutil.copyfile(input_path, output_path)
            return "Aligner not found. Saved unaligned FASTA."
        except Exception as exc:
            shutil.copyfile(input_path, output_path)
            return f"Alignment failed: {exc}"
        shutil.copyfile(input_path, output_path)
        return "Alignment failed unexpectedly. Saved unaligned FASTA."

    def _align(self) -> None:
        gene_id = self.gene_combo.currentText().strip()
        if not gene_id:
            QtWidgets.QMessageBox.warning(self, "Align", "Please select a GeneID.")
            return
        seq_col = self._sequence_column()
        seq_label = self._sequence_label()
        manual_ref = self._reference_text()
        records = self._collect_sequences(gene_id, seq_col)
        if not records and not self.reference_path and not manual_ref:
            QtWidgets.QMessageBox.warning(
                self,
                "Align",
                f"No {seq_label} sequences found for this GeneID and no reference provided.",
            )
            return

        self.mutations_csv_path = None
        self.download_mut_btn.setEnabled(False)
        self.current_tmpdir = Path(tempfile.mkdtemp(prefix="align_annotations_"))
        tmpdir = self.current_tmpdir
        input_path = tmpdir / "input.fasta"
        output_path = tmpdir / "alignment.fasta"
        self._write_input_fasta(records, input_path, manual_ref)
        msg = self._perform_alignment(input_path, output_path)
        self.alignment_path = output_path
        self.download_btn.setEnabled(True)
        snippy_msg = ""
        if seq_col == "Nucleotide_Sequence" and (self.reference_path or manual_ref):
            sample_records = self._collect_sample_sequences(gene_id, seq_col)
            reference_fasta = self._prepare_reference_file(manual_ref, tmpdir)
            snippy_msg, self.mutations_csv_path = self._run_snippy(
                gene_id, sample_records, reference_fasta, tmpdir
            )
            if self.mutations_csv_path and self.mutations_csv_path.exists():
                self.download_mut_btn.setEnabled(True)
        elif seq_col != "Nucleotide_Sequence" and (self.reference_path or manual_ref):
            snippy_msg = "Snippy skipped: amino acid mode selected (requires nucleotide)."
        else:
            snippy_msg = "Snippy skipped: no reference provided."
        ref_notes = []
        if self.reference_path:
            ref_notes.append("reference file")
        if manual_ref:
            ref_notes.append("pasted reference")
        ref_extra = f" Reference added: {', '.join(ref_notes)}." if ref_notes else ""
        self.status_label.setText(
            f"{msg} Sequences included: {len(records)} {seq_label} entries.{ref_extra}\n{snippy_msg}"
        )

    def _download_alignment(self) -> None:
        if not self.alignment_path or not self.alignment_path.exists():
            QtWidgets.QMessageBox.information(self, "Download", "No alignment is available yet.")
            return
        dest, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save alignment",
            "alignment.fasta",
            "FASTA files (*.fasta *.fa *.fna);;All files (*)",
        )
        if not dest:
            return
        try:
            shutil.copyfile(self.alignment_path, dest)
            QtWidgets.QMessageBox.information(self, "Download", f"Alignment saved to {dest}")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Download", f"Could not save alignment: {exc}")

    def _download_mutations(self) -> None:
        if not self.mutations_csv_path or not self.mutations_csv_path.exists():
            QtWidgets.QMessageBox.information(self, "Download", "No mutations CSV is available yet.")
            return
        dest, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save mutations CSV",
            "mutations.csv",
            "CSV files (*.csv);;All files (*)",
        )
        if not dest:
            return
        try:
            shutil.copyfile(self.mutations_csv_path, dest)
            QtWidgets.QMessageBox.information(self, "Download", f"Mutations saved to {dest}")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Download", f"Could not save mutations: {exc}")
