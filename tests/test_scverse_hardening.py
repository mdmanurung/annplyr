from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap


def test_raw_source_is_explicit_read_only_obs_axis_source(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.raw = adata.copy()
    adata.X = np.zeros_like(adata.X)

    filtered = adata.ap.filter(raw=ap.col("g0") > 2)
    assert filtered.obs_names.tolist() == ["c3", "c4"]

    mutated = adata.ap.mutate(raw={"raw_g3": ap.col("g3")})
    assert mutated.obs["raw_g3"].tolist() == [2.0, 1.0, 0.0, 8.0, 4.0]
    np.testing.assert_array_equal(adata.X, np.zeros_like(adata.X))

    summary = adata.ap.summarize(raw={"mean_raw_g0": ap.mean("g0")}, by="batch")
    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "mean_raw_g0": [2.0, 2.5]}),
    )

    wide = adata.ap.to_df(obs=["batch"], raw=["g0"])
    assert wide.columns.tolist() == ["batch", "raw_g0"]
    assert wide["raw_g0"].tolist() == [1.0, 2.0, 0.0, 5.0, 3.0]

    tidy = adata.ap.to_tidy(obs=["batch"], raw=["g3"])
    assert tidy.columns.tolist() == ["obs_name", "feature", "value", "batch"]
    assert tidy.loc[(tidy["obs_name"] == "c3") & (tidy["feature"] == "g3"), "value"].item() == 8.0


def test_as_frame_exposes_scverse_sources_with_controlled_selection(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.raw = adata.copy()
    adata.obsp["distances"] = np.arange(adata.n_obs * adata.n_obs).reshape(adata.n_obs, adata.n_obs)
    adata.varp["correlations"] = np.arange(adata.n_vars * adata.n_vars).reshape(adata.n_vars, adata.n_vars)
    adata.uns["qc"] = pd.DataFrame({"metric": ["cells", "genes"], "value": [adata.n_obs, adata.n_vars]})

    raw = adata.ap.as_frame("raw", select=["g0"])
    assert raw.index.tolist() == adata.obs_names.tolist()
    assert raw.columns.tolist() == ["g0"]

    obsp = adata.ap.as_frame("obsp", key="distances", select=["c1"])
    assert obsp.index.tolist() == adata.obs_names.tolist()
    assert obsp.columns.tolist() == ["c1"]
    assert obsp.loc["c3", "c1"] == adata.obsp["distances"][3, 1]

    varp = adata.ap.as_frame("varp", key="correlations", select=["g0"])
    assert varp.index.tolist() == adata.var_names.tolist()
    assert varp.columns.tolist() == ["g0"]
    assert varp.loc["g3", "g0"] == adata.varp["correlations"][3, 0]

    uns = adata.ap.as_frame("uns", key="qc")
    pd.testing.assert_frame_equal(uns, adata.uns["qc"])

    pulled = adata.ap.pull(obsp={"distances": "c1"})
    assert pulled.index.tolist() == adata.obs_names.tolist()
    assert pulled.name == "c1"


def test_matrix_exports_accept_public_materialization_budgets(dense_adata: AnnData) -> None:
    with pytest.raises(ap.AnnplyrError, match="materialize 20 matrix values"):
        dense_adata.ap.to_tidy(allow_all_features=True, max_matrix_values=10)

    tidy = dense_adata.ap.to_tidy(allow_all_features=True, max_matrix_values=20)
    assert len(tidy) == dense_adata.n_obs * dense_adata.n_vars

    with pytest.raises(ap.AnnplyrError, match="materialize 10 matrix values"):
        dense_adata.ap.to_df(x=["g0", "g3"], max_matrix_values=9)

    with pytest.raises(ap.AnnplyrError, match="materialize 10 matrix values"):
        dense_adata.ap.pivot_longer(x=["g0", "g3"], max_matrix_values=9)


def test_missing_raw_source_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="raw"):
        dense_adata.ap.to_df(raw=["g0"])
