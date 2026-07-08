from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

@dataclass
class DataInput:
    genomes: List[Path]
    marker_sequences: List[Path]
    reference_locus_sequences: List[Path]
    output: Path
    snippy_reference: Optional[Path]
    bakta_db: Optional[str]
    dmnd_threads: int
    bakta_threads: int
    window_length: int
    diamond: bool
    bakta: bool
    total_cores: int
    marker_categories: Dict[str, str] = field(default_factory=dict)              # stem -> category
    reference_locus_categories: Dict[str, str] = field(default_factory=dict)     # stem -> category
    fastgear_config: Dict[str, Any] = field(default_factory=dict)
    mobsuite_config: Dict[str, Any] = field(default_factory=dict)
    chewbbaca_config: Dict[str, Any] = field(default_factory=dict)
    assembly_qc_config: Dict[str, Any] = field(default_factory=dict)
    run_id: Optional[str] = None

    def __post_init__(self):
        # Always run Bakta and Diamond; require their inputs.
        self.bakta = True
        self.diamond = True

        self.genomes = self._norm_paths(self.genomes, name="genomes", must_exist=True, must_be_file=True, allow_empty=False)
        self.marker_sequences = self._norm_paths(self.marker_sequences, name="marker_sequences", must_exist=True, must_be_file=True, allow_empty=False)
        self.reference_locus_sequences = self._norm_paths(self.reference_locus_sequences, name="reference_locus_sequences", must_exist=True, must_be_file=True, allow_empty=False)
        self.snippy_reference = self._norm_optional_path(self.snippy_reference, name="snippy_reference")

        self.output = self._norm_dir(self.output, name="output", create=True)

        if not (self.bakta_db or "").strip():
            raise ValueError("bakta_db: required")
        p = Path(self.bakta_db).expanduser()
        if not p.is_dir():
            raise ValueError(f"bakta_db: directory not found: {p}")
        self.bakta_db = str(p.resolve())

        self.dmnd_threads = self._pos_int(self.dmnd_threads, "dmnd_threads", 1)
        self.bakta_threads = self._pos_int(self.bakta_threads, "bakta_threads", 1)
        self.window_length = self._pos_int(self.window_length, "window_length", 1)
        self.total_cores = self._pos_int(self.total_cores, "total_cores", 1)

        self.diamond = True
        self.bakta = True

        if not self.marker_sequences:
            raise ValueError("marker_sequences: at least one path required")
        if not self.reference_locus_sequences:
            raise ValueError("reference_locus_sequences: at least one path required")

        m_stems = [p.stem for p in self.marker_sequences]
        r_stems = [p.stem for p in self.reference_locus_sequences]
        self._ensure_unique_stems(m_stems, "marker_sequences")
        self._ensure_unique_stems(r_stems, "reference_locus_sequences")

        if self.snippy_reference is None:
            self.snippy_reference = sorted(self.genomes, key=lambda p: (p.stem, str(p)))[0]
        elif self.snippy_reference not in self.genomes:
            raise ValueError(f"snippy_reference: must be one of the input genomes: {self.snippy_reference}")

        self.marker_categories = self._norm_category_map(self.marker_categories, m_stems, "marker_categories")
        self.reference_locus_categories = self._norm_category_map(self.reference_locus_categories, r_stems, "reference_locus_categories")
        self.fastgear_config = self._norm_json_map(self.fastgear_config, "fastgear_config")
        self.mobsuite_config = self._norm_json_map(self.mobsuite_config, "mobsuite_config")
        self.chewbbaca_config = self._norm_json_map(self.chewbbaca_config, "chewbbaca_config")
        self.assembly_qc_config = self._norm_json_map(self.assembly_qc_config, "assembly_qc_config")

        if self.run_id is not None:
            run_id = str(self.run_id).strip()
            self.run_id = run_id if run_id else None

    # helpers
    def _norm_paths(self, items: Any, *, name: str, must_exist: bool, must_be_file: bool, allow_empty: bool) -> List[Path]:
        if items is None:
            items = []
        if not isinstance(items, (list, tuple, set)):
            items = [items]
        out: List[Path] = []
        seen: set[str] = set()
        for x in items:
            p = Path(x).expanduser()
            if must_exist and not p.exists():
                raise ValueError(f"{name}: not found: {p}")
            if must_exist and must_be_file and not p.is_file():
                raise ValueError(f"{name}: expected file: {p}")
            key = str(p.resolve() if p.exists() else p)
            if key in seen:
                continue
            seen.add(key)
            out.append(p.resolve() if p.exists() else p)
        if not allow_empty and not out:
            raise ValueError(f"{name}: at least one path required")
        return out

    def _norm_optional_path(self, value: Any, *, name: str) -> Optional[Path]:
        if value in (None, "", False):
            return None
        p = Path(value).expanduser()
        if not p.exists():
            raise ValueError(f"{name}: not found: {p}")
        if not p.is_file():
            raise ValueError(f"{name}: expected file: {p}")
        return p.resolve()

    def _norm_dir(self, path_like: Any, *, name: str, create: bool = False) -> Path:
        p = Path(path_like).expanduser()
        if not p.exists() and create:
            try:
                p.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise ValueError(f"{name}: cannot create directory {p}: {e}")
        if not p.exists() or not p.is_dir():
            raise ValueError(f"{name}: directory not found: {p}")
        return p.resolve()

    def _pos_int(self, v: Any, name: str, minimum: int = 1, maximum: Optional[int] = None) -> int:
        try:
            iv = int(v)
        except Exception:
            raise ValueError(f"{name}: must be an integer")
        if iv < minimum:
            raise ValueError(f"{name}: must be ≥ {minimum}")
        if maximum is not None and iv > maximum:
            raise ValueError(f"{name}: must be ≤ {maximum}")
        return iv

    def _ensure_unique_stems(self, stems: List[str], name: str) -> None:
        if len(stems) != len(set(stems)):
            dupes = [s for s in set(stems) if stems.count(s) > 1]
            raise ValueError(f"{name}: file names (stems) must be unique: {', '.join(sorted(dupes))}")

    def _norm_category_map(self, cats: Optional[Mapping[str, str]], stems: List[str], name: str) -> Dict[str, str]:
        stems_set = set(stems)
        if cats is None:
            cats = {}
        if not isinstance(cats, Mapping):
            raise ValueError(f"{name}: must be a mapping of file-stem -> category")
        out: Dict[str, str] = {}
        for k, v in cats.items():
            ks = str(k)
            if ks not in stems_set:
                raise ValueError(f"{name}: unknown key '{ks}' (expected one of: {', '.join(sorted(stems_set))})")
            out[ks] = "" if v is None else str(v)
        for s in stems:
            out.setdefault(s, "")
        return out

    def _norm_json_map(self, value: Any, name: str) -> Dict[str, Any]:
        if value is None:
            return {}
        if not isinstance(value, Mapping):
            raise ValueError(f"{name}: must be a mapping")
        return {str(k): v for k, v in value.items()}
