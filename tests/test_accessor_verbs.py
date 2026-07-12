from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

import annplyr as ap


def test_import_registers_ap_accessor(dense_adata: AnnData) -> None:
    assert hasattr(dense_adata, "ap")
    assert dense_adata.ap._obj is dense_adata
    assert ap.obs_names is not None
    assert ap.var_names is not None


def test_filter_select_and_arrange_match_manual_indexing(dense_adata: AnnData) -> None:
    filtered = dense_adata.ap.filter(
        obs=ap.col("batch") == "A",
        x=ap.col("g0") > 20,
        obsm={"X_pca": ap.col("0") >= 0},
        var=ap.col("feature_type") == "rna",
        varm={"loadings": ap.col("0") > 0},
        layer="counts",
    )

    assert filtered.obs_names.tolist() == ["c4"]
    assert filtered.var_names.tolist() == ["g0", "g3"]
    np.testing.assert_array_equal(filtered.X, dense_adata[["c4"], ["g0", "g3"]].X)

    selected = dense_adata.ap.select(
        obs=ap.starts_with("cell"),
        var=ap.col("chrom"),
        x=ap.col("g0", "g3"),
    )

    assert selected.obs.columns.tolist() == ["cell_type"]
    assert selected.var.columns.tolist() == ["chrom"]
    assert selected.var_names.tolist() == ["g0", "g3"]

    arranged = dense_adata.ap.arrange(obs=ap.desc("score"), var=ap.col("length"))

    assert arranged.obs_names.tolist() == ["c3", "c4", "c1", "c0", "c2"]
    assert arranged.var_names.tolist() == ["g3", "g0", "g2", "g1"]


def test_select_everything_ignores_synthetic_name_columns_and_preserves_x_order(dense_adata: AnnData) -> None:
    selected = dense_adata.ap.select(
        obs=ap.everything(),
        var=ap.everything(),
        x=ap.col("g3", "g0"),
    )

    assert selected.obs.columns.tolist() == dense_adata.obs.columns.tolist()
    assert selected.var.columns.tolist() == dense_adata.var.columns.tolist()
    assert selected.var_names.tolist() == ["g3", "g0"]
    np.testing.assert_array_equal(selected.X, dense_adata[:, ["g3", "g0"]].X)


def test_slice_helpers_default_to_obs_axis(dense_adata: AnnData) -> None:
    assert dense_adata.ap.slice(1, 3).obs_names.tolist() == ["c1", "c3"]
    assert dense_adata.ap.slice_head(n=2).obs_names.tolist() == ["c0", "c1"]
    assert dense_adata.ap.slice_tail(n=2, axis="var").var_names.tolist() == ["g2", "g3"]
    assert dense_adata.ap.slice_min(ap.col("score"), n=2).obs_names.tolist() == ["c2", "c0"]
    assert dense_adata.ap.slice_max(ap.col("length"), n=1, axis="var").var_names.tolist() == ["g1"]

    sampled = dense_adata.ap.slice_sample(n=2, random_state=7)
    assert sampled.n_obs == 2
    assert sampled.n_vars == dense_adata.n_vars


def test_mutate_writes_obs_and_var_columns_from_read_only_sources(dense_adata: AnnData) -> None:
    original_x = dense_adata.X.copy()
    original_obsm = dense_adata.obsm["X_pca"].copy()
    original_varm = dense_adata.varm["loadings"].copy()

    mutated = dense_adata.ap.mutate(
        obs={"double_score": ap.col("score") * 2},
        var={"length_kb": ap.col("length") / 1000},
        x={"g0_counts": ap.col("g0")},
        obsm={"X_pca": {"pc0": ap.col("0")}},
        varm={"loadings": {"loading0": ap.col("0")}},
        layer="counts",
    )

    assert "double_score" not in dense_adata.obs
    assert mutated is not dense_adata
    assert mutated.obs["double_score"].tolist() == [2.0, 4.0, 1.0, 6.0, 5.0]
    assert mutated.obs["g0_counts"].tolist() == [10.0, 20.0, 0.0, 50.0, 30.0]
    assert mutated.obs["pc0"].tolist() == [0.1, -0.2, 0.3, 1.0, 0.5]
    assert mutated.var["length_kb"].tolist() == [0.1, 0.2, 0.15, 0.08]
    assert mutated.var["loading0"].tolist() == [0.5, -0.1, 0.3, 1.2]
    np.testing.assert_array_equal(mutated.X, original_x)
    np.testing.assert_array_equal(mutated.obsm["X_pca"], original_obsm)
    np.testing.assert_array_equal(mutated.varm["loadings"], original_varm)


def test_grouped_iteration_summarize_count_and_by_arguments(dense_adata: AnnData) -> None:
    grouped = dense_adata.ap.group_by(obs="batch")

    groups = list(grouped)
    assert [key for key, _ in groups] == [{"batch": "A"}, {"batch": "B"}]
    assert [group.obs_names.tolist() for _, group in groups] == [["c0", "c1", "c4"], ["c2", "c3"]]

    summary = grouped.summarize(obs={"mean_score": ap.mean("score"), "cells": ap.n()})
    expected = pd.DataFrame({"batch": ["A", "B"], "mean_score": [11 / 6, 7 / 4], "cells": [3, 2]})
    pd.testing.assert_frame_equal(summary.reset_index(drop=True), expected)

    counted = dense_adata.ap.count(by="batch")
    pd.testing.assert_frame_equal(
        counted.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "n": [3, 2]}),
    )

    one_off = dense_adata.ap.summarise(obs={"max_score": ap.max("score")}, by="batch")
    pd.testing.assert_frame_equal(
        one_off.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "max_score": [2.5, 3.0]}),
    )


def test_summarize_reads_matrix_sources_with_by_arguments(dense_adata: AnnData) -> None:
    obs_summary = dense_adata.ap.summarize(
        x={"mean_g0": ap.mean("g0")},
        obsm={"X_pca": {"mean_pc0": ap.mean("0")}},
        by="batch",
    )
    pd.testing.assert_frame_equal(
        obs_summary.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "mean_g0": [2.0, 2.5], "mean_pc0": [0.13333333333333333, 0.65]}),
    )

    var_summary = dense_adata.ap.summarize(
        varm={"loadings": {"max_loading0": ap.max("0")}},
        by="feature_type",
    )
    pd.testing.assert_frame_equal(
        var_summary.reset_index(drop=True),
        pd.DataFrame({"feature_type": ["rna", "protein"], "max_loading0": [1.2, 0.3]}),
    )


def test_pull_to_df_and_to_tidy_are_plot_friendly(dense_adata: AnnData) -> None:
    pulled = dense_adata.ap.pull(x="g3")
    assert pulled.name == "g3"
    assert pulled.tolist() == [2.0, 1.0, 0.0, 8.0, 4.0]

    wide = dense_adata.ap.to_df(obs=["batch"], x=["g0", "g3"], obsm={"X_pca": ["0"]})
    assert wide.index.tolist() == dense_adata.obs_names.tolist()
    assert wide.columns.tolist() == ["batch", "g0", "g3", "X_pca_0"]
    assert wide.loc["c3", "X_pca_0"] == 1.0

    tidy = dense_adata.ap.to_tidy(obs=["batch"], x=["g0", "g3"])
    assert tidy.columns.tolist() == ["obs_name", "feature", "value", "batch"]
    assert len(tidy) == dense_adata.n_obs * 2
    assert tidy.loc[(tidy["obs_name"] == "c4") & (tidy["feature"] == "g3"), "value"].item() == 4.0


def test_sparse_anndata_uses_same_expression_semantics(sparse_adata: AnnData) -> None:
    filtered = sparse_adata.ap.filter(x=ap.col("g3") >= 40, layer="counts")
    assert filtered.obs_names.tolist() == ["c3", "c4"]

    mutated = sparse_adata.ap.mutate(x={"g1_counts": ap.col("g1")}, layer="counts")
    assert mutated.obs["g1_counts"].tolist() == [0.0, 30.0, 10.0, 20.0, 0.0]
