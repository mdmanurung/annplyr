from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import annplyr as ap
from annplyr._aggregate import aggregate_cells_adata


@pytest.fixture
def cells() -> AnnData:
    rng = np.random.default_rng(0)
    counts = rng.poisson(2, size=(12, 5)).astype(np.int64)
    obs = pd.DataFrame(
        {
            "donor": pd.Categorical(["d1", "d2", "d3"] * 4),
            "cell_type": ["T"] * 6 + ["B"] * 6,
        },
        index=[f"c{i}" for i in range(12)],
    )
    adata = AnnData(X=np.log1p(counts).astype(np.float32), obs=obs)
    adata.layers["counts"] = sp.csr_matrix(counts)
    return adata


@pytest.mark.parametrize("fun", ["sum", "mean"])
def test_aggregate_matches_a_pandas_groupby(cells: AnnData, fun: str) -> None:
    pb = aggregate_cells_adata(cells, "donor", "cell_type", source="counts", fun=fun)

    dense = pd.DataFrame(cells.layers["counts"].toarray(), columns=cells.var_names)
    reference = getattr(
        dense.groupby([cells.obs["donor"].to_numpy(), cells.obs["cell_type"].to_numpy()], sort=False), fun
    )()
    got = pd.DataFrame(pb.X.toarray(), index=pd.MultiIndex.from_frame(pb.obs[["donor", "cell_type"]].astype(str)))
    np.testing.assert_allclose(got.to_numpy(), reference.loc[got.index].to_numpy())
    assert pb.obs["n_cells"].tolist() == [2] * 6
    assert sp.issparse(pb.X)
    assert pb.var_names.equals(cells.var_names)
    if fun == "sum":
        assert pb.X.dtype == np.int64


def test_aggregate_needs_a_known_source(cells: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="counts"):
        aggregate_cells_adata(cells, "donor", source="cnts")


def test_a_float32_sum_accumulates_in_float64(cells: AnnData) -> None:
    cells.X = np.full((cells.n_obs, cells.n_vars), 0.1, dtype=np.float32)
    pb = aggregate_cells_adata(cells, "donor", source="x")
    assert pb.X.dtype == np.float64
    np.testing.assert_allclose(np.asarray(pb.X)[0, 0], 4 * np.float64(np.float32(0.1)), rtol=1e-15)


def test_aggregate_keeps_columns_constant_within_each_group(cells: AnnData) -> None:
    cells.obs["condition"] = pd.Categorical(np.where(cells.obs["donor"] == "d1", "ctrl", "stim"))
    cells.obs["varies"] = np.arange(cells.n_obs)
    pb = aggregate_cells_adata(cells, "donor", "cell_type", source="counts")
    assert "varies" not in pb.obs
    assert isinstance(pb.obs["condition"].dtype, pd.CategoricalDtype)
    expected = np.where(pb.obs["donor"].astype(str) == "d1", "ctrl", "stim")
    assert pb.obs["condition"].astype(str).tolist() == expected.tolist()
