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

    filtered = adata.ap.filter(ap.col("raw:g0") > 2)
    assert filtered.obs_names.tolist() == ["c3", "c4"]

    mutated = adata.ap.mutate({"raw_g3": ap.col("raw:g3")})
    assert mutated.obs["raw_g3"].tolist() == [2.0, 1.0, 0.0, 8.0, 4.0]
    np.testing.assert_array_equal(adata.X, np.zeros_like(adata.X))

    summary = adata.ap.summarize({"mean_raw_g0": ap.mean("raw:g0")}, by="batch")
    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "mean_raw_g0": [2.0, 2.5]}),
    )

    wide = adata.ap.to_wide("batch", "raw:g0")
    assert wide.columns.tolist() == ["batch", "g0"]
    assert wide["g0"].tolist() == [1.0, 2.0, 0.0, 5.0, 3.0]

    tidy = adata.ap.join_features("raw:g3")
    assert {"batch", ".cell", ".feature", ".abundance_raw"} <= set(tidy.columns)
    assert tidy.loc[(tidy[".cell"] == "c3") & (tidy[".feature"] == "g3"), ".abundance_raw"].item() == 8.0


def test_to_wide_and_pull_expose_scverse_sources(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.raw = adata.copy()
    adata.obsp["distances"] = np.arange(adata.n_obs * adata.n_obs).reshape(adata.n_obs, adata.n_obs)

    raw = adata.ap.to_wide("raw:g0")
    assert raw.index.tolist() == adata.obs_names.tolist()
    assert raw.columns.tolist() == ["g0"]

    pulled = adata.ap.pull("obsp:distances:c1")
    assert pulled.index.tolist() == adata.obs_names.tolist()
    assert pulled.name == "c1"
    assert pulled.loc["c3"] == adata.obsp["distances"][3, 1]


def test_matrix_exports_accept_public_materialization_budgets(dense_adata: AnnData) -> None:
    features = [f"x:{name}" for name in dense_adata.var_names]
    with pytest.raises(ap.AnnplyrError, match="materialize 20 matrix values"):
        dense_adata.ap.join_features(*features, max_matrix_values=10)

    tidy = dense_adata.ap.join_features(*features, max_matrix_values=20)
    assert len(tidy) == dense_adata.n_obs * dense_adata.n_vars

    with pytest.raises(ap.AnnplyrError, match="materialize 10 matrix values"):
        dense_adata.ap.to_wide("x:g0", "x:g3", max_matrix_values=9)

    with pytest.raises(ap.AnnplyrError, match="materialize 10 matrix values"):
        dense_adata.ap.pivot_longer("x:g0", "x:g3", max_matrix_values=9)


def test_missing_raw_source_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="raw"):
        dense_adata.ap.to_wide("raw:g0")
