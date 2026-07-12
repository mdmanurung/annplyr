from __future__ import annotations

import pandas as pd
from anndata import AnnData

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


def test_grouped_var_count_and_keys(dense_adata: AnnData) -> None:
    grouped = dense_adata.ap.group_by(var="feature_type")

    pd.testing.assert_frame_equal(
        grouped.count().reset_index(drop=True),
        pd.DataFrame({"feature_type": ["rna", "protein"], "n": [3, 1]}),
    )
    assert grouped.group_data()[".rows"].tolist() == [["g0", "g1", "g3"], ["g2"]]
