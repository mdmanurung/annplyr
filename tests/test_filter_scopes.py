"""Independent references for the scoped filtering conveniences."""

from __future__ import annotations

import narwhals as nw
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap
from annplyr._sources import RequestPlanner


@pytest.mark.parametrize("surface", ["root", "accessor", "chain"])
def test_metadata_and_expression_chaining(dense_adata: AnnData, surface: str) -> None:
    dense_adata.raw = dense_adata
    dense_adata.obs["g0"] = -1.0
    mask = (dense_adata.obs["score"].to_numpy() >= 1) & (dense_adata.X[:, 0] > 1)
    if surface == "root":
        result = ap.filter_exp(ap.filter_obs(dense_adata, ap.col("score") >= 1), ap.col("raw:g0") > 1)
    elif surface == "accessor":
        result = dense_adata.ap.filter_obs(ap.col("score") >= 1).ap.filter_exp(ap.col("raw:g0") > 1)
    else:
        result = ap.filter_obs(dense_adata, ap.col("score") >= 1).filter_exp(ap.col("raw:g0") > 1)
    expected = dense_adata[mask].copy()
    pd.testing.assert_frame_equal(result.obs, expected.obs)
    np.testing.assert_array_equal(result.X, expected.X)
    assert isinstance(result, ap.TidyAnnData)
    result.obs["g0"] = 100.0
    assert (dense_adata.obs["g0"] == -1).all()


@pytest.mark.parametrize("source", ["x", "raw", "counts"])
@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
def test_expression_sources_and_projection(dense_adata: AnnData, source: str, storage: str) -> None:
    dense_adata.raw = dense_adata
    if storage != "dense":
        dense_adata.X = getattr(sparse, f"{storage}_matrix")(dense_adata.X)
        dense_adata.layers["counts"] = getattr(sparse, f"{storage}_matrix")(dense_adata.layers["counts"])
    values = (
        dense_adata.layers[source] if source == "counts" else dense_adata.raw.X if source == "raw" else dense_adata.X
    )
    column = values[:, 0].toarray().ravel() if sparse.issparse(values) else values[:, 0]
    expected = dense_adata.obs_names[column > 1]
    result = dense_adata.ap.filter_exp(ap.col(f"{source}:g0") > 1, max_matrix_values=dense_adata.n_obs)
    assert result.obs_names.equals(expected)
    if storage != "dense":
        assert result.X.format == storage


@pytest.mark.parametrize(
    ("verb", "reference"),
    [
        ("filter_obs", "x:g0"),
        ("filter_obs", "var:length"),
        ("filter_exp", "score"),
        ("filter_exp", "X_pca:0"),
        ("filter_exp", "var:length"),
    ],
)
def test_scope_rejected_before_reads(
    dense_adata: AnnData, monkeypatch: pytest.MonkeyPatch, verb: str, reference: str
) -> None:
    def unread(*args: object, **kwargs: object) -> None:
        pytest.fail("an invalid scoped filter attempted a matrix read")

    monkeypatch.setattr(RequestPlanner, "execute", unread)
    with pytest.raises(ap.IncompatibleAxisError, match=verb):
        getattr(dense_adata.ap, verb)(ap.col(reference) > 0)


def test_opaque_expression_requires_explicit_matrix_source(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError):
        dense_adata.ap.filter_exp(nw.col("g0") > 1)
    result = dense_adata.ap.filter_exp(nw.col("g0") > 1, source="x")
    assert result.obs_names.equals(dense_adata.obs_names[dense_adata.X[:, 0] > 1])
    with pytest.raises(ap.AnnplyrError, match="exceeds max_matrix_values"):
        dense_adata.ap.filter_exp(nw.col("g0") > 1, source="x", max_matrix_values=dense_adata.n_obs)
    with pytest.raises(ap.IncompatibleAxisError):
        dense_adata.ap.filter_obs(ap.col("score") > 0, source="raw")


def test_raw_only_feature_empty_and_errors(dense_adata: AnnData) -> None:
    dense_adata.raw = dense_adata
    selected = dense_adata[:, ["g1"]].copy()
    assert selected.ap.filter_exp(ap.col("raw:g0") > 1).n_vars == 1
    assert selected.ap.filter_exp(ap.col("raw:g0") > 100).shape == (0, 1)
    with pytest.raises(ap.UnknownColumnError):
        selected.ap.filter_exp(ap.col("x:g0") > 0)
    with pytest.raises(ap.UnknownSourceError):
        selected.ap.filter_exp(ap.col("absent:g1") > 0)
    selected.raw = None
    with pytest.raises(ap.UnknownSourceError):
        selected.ap.filter_exp(ap.col("raw:g0") > 0)


def test_null_and_or_predicates_with_duplicate_names(dense_adata: AnnData) -> None:
    dense_adata.obs["maybe"] = pd.array([1, pd.NA, 0, 2, -1], dtype="Int64")
    dense_adata.obs_names = ["duplicate"] * dense_adata.n_obs
    result = dense_adata.ap.filter_obs((ap.col("maybe") > 0, ap.col("score") >= 1), copy=False)
    np.testing.assert_array_equal(result.X, dense_adata.X[[0, 3]])
    union = dense_adata.ap.filter_obs((ap.col("score") < 1) | (ap.col("score") > 2.5))
    np.testing.assert_array_equal(union.X, dense_adata.X[[2, 3]])


@pytest.mark.parametrize("verb", ["filter_obs", "filter_exp"])
def test_group_local_evaluation_preserves_keys(dense_adata: AnnData, verb: str) -> None:
    column = "score" if verb == "filter_obs" else "x:g0"
    values = (
        dense_adata.obs["score"]
        if verb == "filter_obs"
        else pd.Series(dense_adata.X[:, 0], index=dense_adata.obs_names)
    )
    mask = values > values.groupby(dense_adata.obs["batch"], sort=False).transform("mean")
    positions = [
        i
        for key in pd.unique(dense_adata.obs["batch"])
        for i in range(dense_adata.n_obs)
        if dense_adata.obs["batch"].iloc[i] == key and mask.iloc[i]
    ]
    grouped = getattr(dense_adata.ap.group_by("batch"), verb)(ap.col(column) > ap.mean(column))
    assert grouped.group_vars() == ["batch"]
    assert grouped.ungroup().obs_names.equals(dense_adata.obs_names[positions])


def test_no_predicate_is_the_existing_filter_noop(dense_adata: AnnData) -> None:
    assert dense_adata.ap.filter_obs() is dense_adata
    assert dense_adata.ap.filter_exp() is dense_adata
    assert dense_adata.ap.filter_exp(ap.lit(False), source="x").n_obs == 0


@pytest.mark.parametrize("predicate", [ap.lit(True), ap.mean("score") > 0, nw.lit(True)])
def test_scalar_predicates_broadcast(dense_adata: AnnData, predicate: object) -> None:
    result = dense_adata.ap.filter(predicate)
    assert result.obs_names.equals(dense_adata.obs_names)
    assert dense_adata.ap.filter_obs(ap.lit(False)).n_obs == 0


def test_pipe_preserves_callable_return(dense_adata: AnnData) -> None:
    def mutate_returning_none(adata: AnnData) -> None:
        adata.obs["flag"] = True

    assert dense_adata.ap.pipe(mutate_returning_none) is None
    assert dense_adata.obs["flag"].all()
