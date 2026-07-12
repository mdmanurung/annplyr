from __future__ import annotations

import pytest
from anndata import AnnData

import annplyr as ap


def test_invalid_axis_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="axis"):
        dense_adata.ap.slice_head(axis="cells")


def test_summarize_rejects_mixed_axes_with_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="obs-axis or var-axis"):
        dense_adata.ap.summarize(obs={"cells": ap.n()}, var={"genes": ap.n()})


def test_pull_requires_exactly_one_source(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="exactly one source"):
        dense_adata.ap.pull()

    with pytest.raises(ap.UnknownSourceError, match="exactly one source"):
        dense_adata.ap.pull(obs="batch", x="g0")


def test_pull_requires_single_selected_column(dense_adata: AnnData) -> None:
    with pytest.raises(ap.SelectionError, match="exactly one column"):
        dense_adata.ap.pull(obs=["batch", "score"])


def test_missing_metadata_column_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownColumnError, match="missing"):
        dense_adata.ap.select(obs="missing")


def test_missing_obsm_key_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="obsm.*missing"):
        dense_adata.ap.pull(obsm={"missing": "0"})


def test_group_by_rejects_mixed_axes_with_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="one axis"):
        dense_adata.ap.group_by(obs="batch", var="feature_type")


def test_pipe_keyword_collision_raises_annplyr_error(dense_adata: AnnData) -> None:
    def identity(data: AnnData) -> AnnData:
        return data

    with pytest.raises(ap.AnnplyrError, match="pipe target"):
        dense_adata.ap.pipe((identity, "data"), data=dense_adata)
