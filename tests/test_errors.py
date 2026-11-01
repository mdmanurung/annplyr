from __future__ import annotations

import pytest
from anndata import AnnData

import annplyr as ap


def test_invalid_axis_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="axis"):
        dense_adata.ap.slice_head(axis="cells")


def test_summarize_rejects_mixed_axes_with_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="obs-axis or var-axis"):
        dense_adata.ap.summarize(
            {"cells": ap.n(), "genes": ap.n_distinct("var:__annplyr_var_names__")},
        )


def test_summarize_rejects_a_cross_container_aggregate(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="cannot mix columns from"):
        dense_adata.ap.summarize({"bad": ap.mean("score") + ap.mean("x:g0")})


def test_pull_requires_exactly_one_source(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="exactly one source"):
        dense_adata.ap.pull()

    with pytest.raises(ap.UnknownSourceError, match="exactly one source"):
        dense_adata.ap.pull(obs="batch", x="g0")


def test_pull_requires_single_selected_column(dense_adata: AnnData) -> None:
    with pytest.raises(ap.SelectionError, match="exactly one column"):
        dense_adata.ap.pull(obs=["batch", "score"])


@pytest.mark.parametrize("source", ["obsm", "varm", "obsp", "varp", "uns"])
@pytest.mark.parametrize("mapping", [{}, {"first": "0", "second": "0"}])
def test_pull_requires_exactly_one_mapping_key(dense_adata: AnnData, source: str, mapping: dict[str, str]) -> None:
    with pytest.raises(ap.UnknownSourceError, match="exactly one source key"):
        dense_adata.ap.pull(**{source: mapping})


def test_missing_metadata_column_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownColumnError, match="missing"):
        dense_adata.ap.select(obs="missing")


def test_missing_obsm_key_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="obsm.*missing"):
        dense_adata.ap.pull(obsm={"missing": "0"})


def test_group_by_rejects_a_qualified_key_because_axis_is_given(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="redundant"):
        dense_adata.ap.group_by("var:feature_type")


def test_group_by_rejects_an_unknown_axis(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="obs|var"):
        dense_adata.ap.group_by("batch", axis="cells")


def test_a_grouping_key_naming_the_other_axis_is_rejected(dense_adata: AnnData) -> None:
    """Resolving ``var:`` against an obs frame would group the wrong column."""
    for call in (
        lambda: dense_adata.ap.count("var:feature_type"),
        lambda: dense_adata.ap.add_count("var:feature_type"),
        lambda: dense_adata.ap.summarize({"n": ap.n()}, by="var:feature_type"),
    ):
        with pytest.raises(ap.IncompatibleAxisError, match="names the var axis"):
            call()


def test_a_var_axis_reduce_requires_a_qualified_grouping_key(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="bare name is obs metadata"):
        dense_adata.ap.summarize({"m": ap.max("var:feature_type")}, by="feature_type")


def test_pipe_keyword_collision_raises_annplyr_error(dense_adata: AnnData) -> None:
    def identity(data: AnnData) -> AnnData:
        return data

    with pytest.raises(ap.AnnplyrError, match="pipe target"):
        dense_adata.ap.pipe((identity, "data"), data=dense_adata)
