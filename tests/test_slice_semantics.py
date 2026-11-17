"""``slice`` is 0-based, unlike dplyr, and ``by=`` verbs refuse a conflicting ``axis=``."""

from __future__ import annotations

import pytest
from anndata import AnnData

import annplyr as ap


def test_slice_is_zero_based(dense_adata: AnnData) -> None:
    assert dense_adata.ap.slice(0).obs_names.tolist() == ["c0"]
    assert dense_adata.ap.slice(1).obs_names.tolist() == ["c1"]
    assert ap.slice(dense_adata, -1).obs_names.tolist() == ["c4"]


def test_slice_docstring_warns_about_the_base(dense_adata: AnnData) -> None:
    for doc in (
        ap.slice.__doc__,
        type(dense_adata.ap).slice.__doc__,
        dense_adata.ap.group_by("batch").slice.__doc__,
    ):
        assert doc is not None
        assert "0-based" in doc
        assert "1-based" in doc


@pytest.mark.parametrize(
    ("verb", "args"),
    [
        ("slice_head", ()),
        ("slice_tail", ()),
        ("slice_sample", ()),
        ("slice_min", ("score",)),
        ("slice_max", ("score",)),
    ],
)
@pytest.mark.parametrize("axis", ["obs", "var"])
def test_by_with_axis_raises(dense_adata: AnnData, verb: str, args: tuple, axis: str) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match=r"Drop axis="):
        getattr(dense_adata.ap, verb)(*args, by="batch", axis=axis)
    with pytest.raises(ap.AnnplyrError):
        getattr(ap, verb)(dense_adata, *args, by="batch", axis=axis)


def test_by_without_axis_still_works(dense_adata: AnnData) -> None:
    assert dense_adata.ap.slice_head(1, by="batch").n_obs == dense_adata.obs["batch"].nunique()
    assert dense_adata.ap.slice_head(2, axis="var").var_names.tolist() == dense_adata.var_names[:2].tolist()
