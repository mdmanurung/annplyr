from __future__ import annotations

import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap


def test_tidyselect_helpers_resolve_expected_columns(dense_adata: AnnData) -> None:
    selected = dense_adata.ap.select(
        obs=[
            ap.all_of(["batch"]),
            ap.any_of(["missing", "score"]),
            ap.last_col(),
        ],
        x=ap.num_range("g", [3, 0]),
    )

    assert selected.obs.columns.tolist() == ["batch", "score", "cell_type"]
    assert selected.var_names.tolist() == ["g3", "g0"]


def test_where_selector_selects_columns_by_dtype(dense_adata: AnnData) -> None:
    selected = dense_adata.ap.select(obs=ap.where(pd.api.types.is_numeric_dtype))

    assert selected.obs.columns.tolist() == ["score"]


def test_rename_updates_metadata_columns_and_feature_names(dense_adata: AnnData) -> None:
    renamed = dense_adata.ap.rename(
        obs={"condition": "batch"},
        var={"chr": "chrom"},
        x={"Gene0": "g0"},
    )

    assert renamed.obs.columns.tolist() == ["condition", "score", "cell_type"]
    assert renamed.var.columns.tolist() == ["feature_type", "chr", "length"]
    assert renamed.var_names.tolist() == ["Gene0", "g1", "g2", "g3"]
    assert dense_adata.var_names.tolist() == ["g0", "g1", "g2", "g3"]


def test_relocate_moves_metadata_columns(dense_adata: AnnData) -> None:
    relocated = dense_adata.ap.relocate(obs=["score"], before="batch", var=["length"], after="chrom")

    assert relocated.obs.columns.tolist() == ["score", "batch", "cell_type"]
    assert relocated.var.columns.tolist() == ["feature_type", "chrom", "length"]


def test_distinct_subsets_first_rows_for_unique_metadata(dense_adata: AnnData) -> None:
    distinct = dense_adata.ap.distinct(obs="batch", keep_all=True)

    assert distinct.obs_names.tolist() == ["c0", "c2"]
    assert distinct.n_vars == dense_adata.n_vars


def test_transmute_keeps_only_new_metadata_columns(dense_adata: AnnData) -> None:
    transmuted = dense_adata.ap.transmute(
        obs={"score2": ap.col("score") * 2},
        x={"g0_counts": ap.col("g0")},
        layer="counts",
    )

    assert transmuted.obs.columns.tolist() == ["score2", "g0_counts"]
    assert transmuted.obs["score2"].tolist() == [2.0, 4.0, 1.0, 6.0, 5.0]
    assert transmuted.obs["g0_counts"].tolist() == [10.0, 20.0, 0.0, 50.0, 30.0]


def test_tally_and_add_count(dense_adata: AnnData) -> None:
    tally = dense_adata.ap.tally(by="batch")
    pd.testing.assert_frame_equal(
        tally.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "n": [3, 2]}),
    )

    counted = dense_adata.ap.add_count(by="batch", name="batch_n")
    assert counted.obs["batch_n"].tolist() == [3, 3, 2, 2, 3]
    assert "batch_n" not in dense_adata.obs


def test_everything_ignores_virtual_name_columns_outside_select(dense_adata: AnnData) -> None:
    wide = dense_adata.ap.to_df(obs=ap.everything())
    assert wide.columns.tolist() == dense_adata.obs.columns.tolist()

    counted = dense_adata.ap.count(by=ap.everything())
    assert counted.columns.tolist() == [*dense_adata.obs.columns.tolist(), "n"]
    assert "obs_names" not in counted.columns


def test_rename_rejects_duplicate_sources_and_targets(dense_adata: AnnData) -> None:
    with pytest.raises(ap.DuplicateNameError, match="source"):
        dense_adata.ap.rename(obs={"condition": "batch", "group": "batch"})

    with pytest.raises(ap.DuplicateNameError, match="Duplicate obs"):
        dense_adata.ap.rename(obs={"score": "batch"})


def test_rename_with_rejects_transformed_name_collisions(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["Batch"] = adata.obs["batch"]

    with pytest.raises(ap.DuplicateNameError, match="Duplicate obs"):
        adata.ap.rename_with(str.lower, obs=["batch", "Batch"])


def test_relocate_rejects_missing_anchor_for_target_source(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownColumnError, match="missing"):
        dense_adata.ap.relocate(obs=["score"], before="missing")

    relocated = dense_adata.ap.relocate(obs=["score"], before="batch", var=["length"], after="chrom")
    assert relocated.obs.columns.tolist() == ["score", "batch", "cell_type"]
    assert relocated.var.columns.tolist() == ["feature_type", "chrom", "length"]


def test_last_col_rejects_out_of_bounds_offset(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownColumnError, match="last_col"):
        dense_adata.ap.select(obs=ap.last_col(offset=-1))
