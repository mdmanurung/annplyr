"""Expressions over sparse matrices whose subtype pandas cannot compute with.

Single-cell matrices are almost always `float32`, but pandas implements sparse
binary kernels for `float64`, `int64`, and `bool` only. Combining two projected
matrix columns therefore has to widen the evaluation inputs, without changing
the subtype of exported frames.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap

VALUES = np.array(
    [
        [1.0, 2.0, 0.0],
        [3.0, 0.0, 4.0],
        [0.0, 5.0, 6.0],
        [7.0, 0.0, 0.0],
    ]
)


def sparse_adata(dtype: np.dtype | type) -> AnnData:
    matrix = sparse.csr_matrix(VALUES.astype(dtype))
    return AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"group": ["a", "a", "b", "b"]},
            index=[f"cell_{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32, np.int64])
def test_two_matrix_columns_combine_for_every_sparse_subtype(dtype: np.dtype | type) -> None:
    adata = sparse_adata(dtype)

    mutated = adata.ap.mutate(x={"total": ap.col("g1") + ap.col("g2")})

    assert mutated.obs["total"].tolist() == [3, 3, 5, 7]


@pytest.mark.parametrize("dtype", [np.float32, np.int32])
def test_narrow_sparse_subtypes_combine_in_every_evaluating_verb(dtype: np.dtype | type) -> None:
    adata = sparse_adata(dtype)
    total = ap.col("g1") + ap.col("g2")

    assert adata.ap.filter(x=total > 3).obs_names.tolist() == ["cell_2", "cell_3"]
    assert adata.ap.arrange(x=ap.desc(total)).obs_names.tolist() == [
        "cell_3",
        "cell_2",
        "cell_0",
        "cell_1",
    ]
    assert adata.ap.summarize(x={"mean_total": ap.mean(total)}, by="group")["mean_total"].tolist() == [3.0, 6.0]
    grouped = adata.ap.group_by(obs="group").summarize(x={"mean_total": ap.mean(total)})
    assert grouped["mean_total"].tolist() == [3.0, 6.0]


def test_exports_keep_the_source_sparse_subtype() -> None:
    adata = sparse_adata(np.float32)

    wide = adata.ap.to_df(x=["g1", "g2"])
    long = adata.ap.to_tidy(obs=["group"], x=["g1", "g2"], max_matrix_values=2 * adata.n_obs)

    assert [str(dtype) for dtype in wide.dtypes] == ["Sparse[float32, 0.0]", "Sparse[float32, 0.0]"]
    assert str(long["value"].dtype) == "Sparse[float32, 0.0]"


def test_promotion_does_not_write_through_to_the_source_matrix() -> None:
    adata = sparse_adata(np.float32)

    adata.ap.mutate(x={"total": ap.col("g1") + ap.col("g2")})

    assert adata.X.dtype == np.float32
    assert "total" not in adata.obs
