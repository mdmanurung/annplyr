from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

import annplyr as ap
from annplyr._features import join_features_adata


@pytest.fixture
def cells(dense_adata: AnnData) -> AnnData:
    adata = dense_adata.copy()
    adata.layers["counts"] = np.asarray(adata.X) * 10
    adata.raw = adata
    return adata


def test_long_has_every_obs_column_and_one_abundance_per_source(cells: AnnData) -> None:
    long = join_features_adata(cells, "raw:g0", "raw:g1", "counts:g0", "counts:g1", shape="long")
    assert set(cells.obs.columns) <= set(long.columns)
    assert {".feature", ".abundance_raw", ".abundance_counts"} <= set(long.columns)
    assert len(long) == cells.n_obs * 2
    row = long.set_index([".cell", ".feature"]).loc[("c1", "g1")]
    assert row[".abundance_raw"] == cells.raw[:, "g1"].X[1, 0]
    assert row[".abundance_counts"] == cells.layers["counts"][1, 1]
    assert row["batch"] == cells.obs.loc["c1", "batch"]


def test_wide_adds_features_to_obs_and_refuses_a_clash(cells: AnnData) -> None:
    wide = join_features_adata(cells, "raw:g0", shape="wide")
    np.testing.assert_array_equal(wide.obs["g0"].to_numpy(), np.asarray(cells.raw[:, "g0"].X).ravel())
    with pytest.raises(ap.DuplicateNameError):
        join_features_adata(cells, "raw:g0", "counts:g0", shape="wide")


def test_metadata_is_not_a_feature(cells: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match=r"raw\.score"):
        join_features_adata(cells, "score")


def test_long_keeps_duplicate_cell_names_apart(cells: AnnData) -> None:
    cells.obs_names = ["c0", "c0", "c2", "c3", "c4"]
    long = join_features_adata(cells, "raw:g0", "counts:g0", shape="long")
    assert len(long) == cells.n_obs
    np.testing.assert_array_equal(long[".abundance_counts"].to_numpy(), cells.layers["counts"][:, 0])
