"""Deterministic fixture helpers shared by the ASV suite."""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

for _name in (
    "BLIS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

SEED = 20260808


def dense_adata(n_obs: int, n_vars: int, *, metadata_columns: int = 0) -> AnnData:
    """Build a deterministic float32 dense AnnData fixture."""
    rng = np.random.default_rng(SEED)
    matrix = rng.standard_normal((n_obs, n_vars), dtype=np.float32)
    obs = metadata_frame(n_obs, metadata_columns, rng=rng)
    obs.index = pd.Index([f"o{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=pd.Index([f"g{i}" for i in range(n_vars)]))
    return AnnData(X=matrix, obs=obs, var=var)


def sparse_adata(n_obs: int, n_vars: int, *, matrix_format: str) -> AnnData:
    """Build a deterministic 0.1%-dense float32 CSR or CSC fixture."""
    rng = np.random.default_rng(SEED)
    matrix = sparse.random(
        n_obs,
        n_vars,
        density=0.001,
        format=matrix_format,
        dtype=np.float32,
        random_state=rng,
    )
    matrix.data = rng.standard_normal(matrix.nnz, dtype=np.float32)
    obs = pd.DataFrame(index=pd.Index([f"o{i}" for i in range(n_obs)]))
    var = pd.DataFrame(index=pd.Index([f"g{i}" for i in range(n_vars)]))
    return AnnData(X=matrix, obs=obs, var=var)


def metadata_frame(n_rows: int, n_columns: int, *, rng: np.random.Generator | None = None) -> pd.DataFrame:
    """Build deterministic wide float32 metadata without object columns."""
    rng = np.random.default_rng(SEED) if rng is None else rng
    return pd.DataFrame({f"m{i}": rng.standard_normal(n_rows, dtype=np.float32) for i in range(n_columns)})


def feature_names(n: int) -> list[str]:
    """Return the first ``n`` deterministic feature names."""
    return [f"g{i}" for i in range(n)]


def realize_frame(frame: pd.DataFrame) -> tuple[int, int, float]:
    """Consume an eager table and return a small structural oracle."""
    numeric = frame.select_dtypes(include="number")
    checksum = 0.0 if numeric.empty or frame.empty else float(numeric.iloc[0].sum())
    return frame.shape[0], frame.shape[1], checksum


def assert_columns(actual: Sequence[object], expected: Sequence[str]) -> None:
    """Check output column realization independently of annplyr internals."""
    assert [str(value) for value in actual] == list(expected)


def unwrap_grouped(value: Any) -> Any:
    """Return the underlying AnnData across v0.2 and persistent v0.3 grouping."""
    ungroup = getattr(value, "ungroup", None)
    return ungroup() if callable(ungroup) else value
