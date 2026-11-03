"""``aggregate_cells`` against Scanpy's own pseudobulk, ``sc.get.aggregate``, on PBMC3k."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from annplyr._aggregate import aggregate_cells_adata

sc = pytest.importorskip("scanpy", reason="Scanpy is installed only in the integration environment")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.filterwarnings("ignore::FutureWarning"),
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
]


@pytest.mark.parametrize("fun", ["sum", "mean"])
def test_aggregate_cells_matches_scanpy_on_raw(fun: str) -> None:
    pbmc = sc.datasets.pbmc3k_processed()
    pbmc.obs["donor"] = np.where(np.arange(pbmc.n_obs) % 2 == 0, "d1", "d2")  # synthetic donors

    ours = aggregate_cells_adata(pbmc, "donor", "louvain", source="raw", fun=fun)
    raw = pbmc.raw.to_adata()
    raw.obs = pbmc.obs
    theirs = sc.get.aggregate(raw, by=["donor", "louvain"], func=fun)

    key = ["donor", "louvain"]
    theirs_rows = theirs.obs[key].astype(str).agg("_".join, axis=1).tolist()
    position = {name: i for i, name in enumerate(theirs_rows)}
    order = [position[name] for name in ours.obs_names]
    expected = np.asarray(theirs.layers[fun])[order]
    got = ours.X.toarray() if sparse.issparse(ours.X) else np.asarray(ours.X)
    np.testing.assert_allclose(got, expected, rtol=1e-5, atol=1e-6)
    assert ours.var_names.equals(theirs.var_names)
