#!/usr/bin/env python3
from __future__ import annotations

import re
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd


RAW_METADATA = Path("/Data/ST492/metadata/metadata_ALL_new_merged.csv")
TARGET_METADATA = Path(
    "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.pre_standardize_2026-06-11.csv"
)
DETAIL_COLUMN = "isolation_source_detail"


def _clean(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    text = re.sub(r"\s+", " ", text)
    return text


def normalize_isolation_source_detail(value: object, broad_source: object = "") -> str:
    text = _clean(value)
    key = text.lower()
    key = key.replace("_", " ")
    key = re.sub(r"\s+", " ", key).strip()

    broad = _clean(broad_source).lower()

    if not key or key in {"not applicable", "not collected", "unknown", "na", "n/a", "none"}:
        return "not available"

    if broad == "wastewater" or key in {"influent", "sewage", "wastewater"}:
        return "wastewater"

    if "ventilator" in key or "ecg monitor" in key:
        return "clinical sample"

    if "rectal" in key or "rectum" in key:
        return "fecal"

    if "fecal" in key or key == "feces":
        return "fecal"

    if "nostril" in key:
        return "respiratory"

    if "endotracheal" in key or "tracheal" in key or key == "trachea":
        return "respiratory"

    if "aspiration catheter" in key:
        return "respiratory"

    if "sputum" in key:
        return "respiratory"

    if "bronchoalveolar" in key or key == "lavage":
        return "respiratory"

    if "bronchial wash" in key:
        return "respiratory"

    if key == "bronchial":
        return "respiratory"

    if key in {"respiratory", "lower respiratory tract", "low respiratory tract"}:
        return "respiratory"

    if key == "blood":
        return "blood"

    if key == "urine":
        return "urine"

    if "wound" in key or "abscess" in key or key == "right leg":
        return "wound/skin"

    if key == "skin":
        return "wound/skin"

    if "catheter" in key:
        return "catheter"

    if key == "sterile body site":
        return "sterile site/body fluid"

    if key == "bodily fluid":
        return "sterile site/body fluid"

    if key in {"clinical", "clinical sample", "clinical isolate", "clinical settings", "bacterial isolate specimen"}:
        return "clinical sample"

    return key


def main() -> None:
    raw = pd.read_csv(RAW_METADATA)
    target = pd.read_csv(TARGET_METADATA)

    raw_samples = set(raw["Sample"].astype(str))
    target_samples = set(target["Sample"].astype(str))
    missing_from_raw = sorted(target_samples - raw_samples)
    missing_from_target = sorted(raw_samples - target_samples)
    if missing_from_raw or missing_from_target:
        raise SystemExit(
            "Sample mismatch between metadata files. "
            f"Missing from raw={len(missing_from_raw)}, missing from target={len(missing_from_target)}"
        )

    raw_detail_by_sample = raw.set_index("Sample")["isolation_source"]
    target[DETAIL_COLUMN] = [
        normalize_isolation_source_detail(raw_detail_by_sample.loc[sample], broad_source)
        for sample, broad_source in zip(target["Sample"], target["isolation_source"])
    ]

    if DETAIL_COLUMN in target.columns:
        detail = target.pop(DETAIL_COLUMN)
        insert_at = list(target.columns).index("isolation_source") + 1
        target.insert(insert_at, DETAIL_COLUMN, detail)

    backup_path = TARGET_METADATA.with_name(
        f"{TARGET_METADATA.name}.bak_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    shutil.copy2(TARGET_METADATA, backup_path)
    target.to_csv(TARGET_METADATA, index=False)

    summary = target[DETAIL_COLUMN].value_counts(dropna=False).sort_index()
    print(f"Updated: {TARGET_METADATA}")
    print(f"Backup:  {backup_path}")
    print(f"Rows:    {len(target)}")
    print("\nNormalized isolation_source_detail counts:")
    for label, count in summary.items():
        print(f"{count}\t{label}")


if __name__ == "__main__":
    main()
