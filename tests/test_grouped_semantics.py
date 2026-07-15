from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

import annplyr as ap


def test_group_metadata_methods(dense_adata: AnnData) -> None:
    grouped = dense_adata.ap.group_by(obs=["batch", "cell_type"])

    assert grouped.group_vars() == ["batch", "cell_type"]

    keys = grouped.group_keys()
    pd.testing.assert_frame_equal(
        keys.reset_index(drop=True),
        pd.DataFrame(
            {
                "batch": ["A", "A", "B", "B"],
                "cell_type": ["T", "B", "T", "Mono"],
            }
        ),
    )

    data = grouped.group_data()
    assert data.columns.tolist() == ["batch", "cell_type", ".rows"]
    assert data[".rows"].tolist() == [["c0"], ["c1", "c4"], ["c2"], ["c3"]]
    assert grouped.ungroup() is dense_adata


def test_grouped_obs_mutate_uses_group_local_row_number(dense_adata: AnnData) -> None:
    mutated = dense_adata.ap.group_by(obs="batch").mutate(obs={"within_batch": ap.row_number()})

    assert mutated.obs["within_batch"].tolist() == [1, 2, 1, 2, 3]
    assert "within_batch" not in dense_adata.obs


def test_grouped_filter_slice_and_add_count_are_group_local(dense_adata: AnnData) -> None:
    grouped = dense_adata.ap.group_by(obs="batch")

    filtered = grouped.filter(obs=ap.row_number() == 1)
    assert filtered.obs_names.tolist() == ["c0", "c2"]

    sliced = grouped.slice_head(n=1)
    assert sliced.obs_names.tolist() == ["c0", "c2"]

    counted = grouped.add_count(name="batch_n")
    assert counted.obs["batch_n"].tolist() == [3, 3, 2, 2, 3]


def test_grouped_summarize_accepts_matrix_sources(dense_adata: AnnData) -> None:
    summary = dense_adata.ap.group_by(obs="batch").summarize(x={"mean_g0": ap.mean("g0")})

    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "mean_g0": [2.0, 2.5]}),
    )


def test_grouped_var_count_and_keys(dense_adata: AnnData) -> None:
    grouped = dense_adata.ap.group_by(var="feature_type")

    pd.testing.assert_frame_equal(
        grouped.count().reset_index(drop=True),
        pd.DataFrame({"feature_type": ["rna", "protein"], "n": [3, 1]}),
    )
    assert grouped.group_data()[".rows"].tolist() == [["g0", "g1", "g3"], ["g2"]]


def test_grouped_mutate_handles_duplicate_obs_names() -> None:
    # Label-based .loc assignment breaks when obs_names are non-unique; positional
    # indexing must be used instead.
    x = np.zeros((4, 2), dtype=float)
    obs = pd.DataFrame({"batch": ["A", "A", "B", "B"]}, index=["x", "x", "x", "x"])
    adata = AnnData(X=x, obs=obs)

    mutated = adata.ap.group_by(obs="batch").mutate(obs={"n": ap.row_number()})

    assert mutated.obs["n"].tolist() == [1, 2, 1, 2]


def test_grouped_mutate_rejects_backed_before_writing(tmp_path, dense_adata: AnnData) -> None:
    path = tmp_path / "adata.h5ad"
    dense_adata.write_h5ad(path)
    backed = read_h5ad(path, backed="r")

    try:
        with pytest.raises(ap.AnnplyrError, match="backed"):
            backed.ap.group_by(obs="batch").mutate(obs={"within_batch": ap.row_number()})
        assert "within_batch" not in backed.obs
    finally:
        backed.file.close()
