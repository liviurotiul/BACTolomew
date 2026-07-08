from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd

try:
    import polars as pl
except Exception:
    pl = None  # type: ignore[assignment]


def polars_available() -> bool:
    return pl is not None


def _require_polars() -> None:
    if pl is None:
        raise RuntimeError(
            "Polars backend is required for annotation_explorer_polars, "
            "but 'polars' is not installed in the current environment."
        )


def read_csv_df(
    path: Path | str,
    *,
    dtype: Any | None = None,
    usecols: Sequence[str | int] | None = None,
    nrows: int | None = None,
    sep: str = ",",
    header: int | None | str = "infer",
    names: Sequence[str] | None = None,
    low_memory: bool | None = None,
) -> pd.DataFrame:
    _require_polars()

    # Map common pandas CSV options to polars for faster parsing.
    kwargs: dict[str, Any] = {"separator": sep}
    if nrows is not None:
        kwargs["n_rows"] = int(nrows)
    if usecols is not None:
        kwargs["columns"] = list(usecols)

    if header is None:
        kwargs["has_header"] = False
    elif header in ("infer", 0):
        kwargs["has_header"] = True
    else:
        # Unsupported/rare header mode: fallback for compatibility.
        return pd.read_csv(
            path,
            dtype=dtype,
            usecols=usecols,
            nrows=nrows,
            sep=sep,
            header=header,
            names=names,
            low_memory=low_memory,
        )

    if names is not None:
        kwargs["new_columns"] = list(names)

    df_pl = pl.read_csv(path, **kwargs)
    if dtype is str:
        df_pl = df_pl.with_columns([pl.col(c).cast(pl.Utf8) for c in df_pl.columns])
    return df_pl.to_pandas()


def read_parquet_df(
    path: Path | str,
    *,
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    _require_polars()
    return pl.read_parquet(path, columns=list(columns) if columns is not None else None).to_pandas()
