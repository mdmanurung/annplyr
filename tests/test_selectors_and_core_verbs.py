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


def test_real_virtual_name_columns_are_not_overwritten(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["obs_names"] = ["real"] * adata.n_obs
    adata.var["var_names"] = ["real_var"] * adata.n_vars

    wide = adata.ap.to_df(obs=["obs_names"])
    assert wide["obs_names"].tolist() == ["real"] * adata.n_obs
    assert adata.ap.filter(obs=ap.obs_names == "c0").obs_names.tolist() == ["c0"]
    assert adata.ap.filter(obs_names=ap.obs_names == "c1").obs_names.tolist() == ["c1"]


def test_select_rejects_computed_columns_for_anndata_alignment(dense_adata: AnnData) -> None:
    with pytest.raises(ap.SelectionError, match="computed"):
        dense_adata.ap.select(obs=(ap.col("score") * 2).alias("score2"))

    with pytest.raises(ap.SelectionError, match="computed"):
        dense_adata.ap.select(x=(ap.col("g0") + ap.col("g1")).alias("sum"))


def test_mutate_assignments_are_sequential_and_recycle_scalars(dense_adata: AnnData) -> None:
    mutated = dense_adata.ap.mutate(
        obs={
            "score2": ap.col("score") * 2,
            "score4": ap.col("score2") * 2,
            "cells": ap.n(),
            "mean_score": ap.mean("score"),
        }
    )

    assert mutated.obs["score4"].tolist() == [4.0, 8.0, 2.0, 12.0, 10.0]
    assert mutated.obs["cells"].tolist() == [5] * dense_adata.n_obs
    assert mutated.obs["mean_score"].tolist() == [1.8] * dense_adata.n_obs


def test_arrange_uses_obs_keys_before_matrix_keys(dense_adata: AnnData) -> None:
    arranged = dense_adata.ap.arrange(obs="batch", x="g0", layer="counts")

    assert arranged.obs_names.tolist() == ["c0", "c1", "c4", "c2", "c3"]


def test_advanced_expression_helpers_work_in_summarize_and_mutate(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["maybe_score"] = [1.0, None, 0.5, None, 2.5]

    summary = adata.ap.summarize(
        obs={
            "cell_types": ap.n_distinct("cell_type"),
            "first_score": ap.first("score"),
            "last_score": ap.last("score"),
        },
        by="batch",
    )
    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame(
            {
                "batch": ["A", "B"],
                "cell_types": [2, 2],
                "first_score": [1.0, 0.5],
                "last_score": [2.5, 3.0],
            }
        ),
    )

    mutated = adata.ap.mutate(
        obs={
            "next_score": ap.lead("score"),
            "previous_score": ap.lag("score", default=-1),
            "score_or_zero": ap.replace_na("maybe_score", 0),
            "score_is_missing": ap.is_na("maybe_score"),
            "batch_as_missing": ap.na_if("batch", "A"),
        }
    )

    assert mutated.obs["next_score"].tolist()[:4] == [2.0, 0.5, 3.0, 2.5]
    assert mutated.obs["previous_score"].tolist() == [-1.0, 1.0, 2.0, 0.5, 3.0]
    assert mutated.obs["score_or_zero"].tolist() == [1.0, 0.0, 0.5, 0.0, 2.5]
    assert mutated.obs["score_is_missing"].tolist() == [False, True, False, True, False]
    assert mutated.obs["batch_as_missing"].iloc[2:4].tolist() == ["B", "B"]
    assert mutated.obs["batch_as_missing"].isna().tolist() == [True, True, False, False, True]


def test_lag_default_only_fills_shift_boundaries(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["maybe_score"] = [1.0, None, 0.5, 3.0, 2.5]

    mutated = adata.ap.mutate(obs={"previous": ap.lag("maybe_score", default=0)})

    assert mutated.obs["previous"].iloc[0] == 0
    assert mutated.obs["previous"].iloc[1] == 1.0
    assert pd.isna(mutated.obs["previous"].iloc[2])


def test_summarize_matrix_source_keeps_feature_when_group_name_collides(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["g0"] = [1, 1, 2, 2, 1]

    summary = adata.ap.summarize(x={"mean_g0": ap.mean("g0")}, by="g0")

    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame({"g0": [1, 2], "mean_g0": [2.0, 2.5]}),
    )


def test_rank_and_cumulative_helpers_are_axis_ordered(dense_adata: AnnData) -> None:
    mutated = dense_adata.ap.mutate(
        obs={
            "score_rank": ap.min_rank("score"),
            "score_dense_rank": ap.dense_rank("score"),
            "cum_score": ap.cum_sum("score"),
            "cum_max_score": ap.cum_max("score"),
        }
    )

    assert mutated.obs["score_rank"].tolist() == [2, 3, 1, 5, 4]
    assert mutated.obs["score_dense_rank"].tolist() == [2, 3, 1, 5, 4]
    assert mutated.obs["cum_score"].tolist() == [1.0, 3.0, 3.5, 6.5, 9.0]
    assert mutated.obs["cum_max_score"].tolist() == [1.0, 2.0, 2.0, 3.0, 3.0]
