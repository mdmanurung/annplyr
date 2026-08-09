from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse
from tests.integration.representative_fixture import make_representative_adata, write_fixture

import annplyr as ap

pytestmark = [
    pytest.mark.integration,
    pytest.mark.filterwarnings("ignore:Observation names are not unique"),
    pytest.mark.filterwarnings("ignore:Variable names are not unique"),
]

EXPECTED_GROUPED_POSITIONS = np.array([0, 5, 2, 4, 1, 7, 6, 3])


@pytest.fixture(scope="module")
def serialized_fixture(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    return write_fixture(tmp_path_factory.mktemp("representative-fixture"))


def _assert_matrix_equal(actual: Any, expected: Any, *, path: str) -> None:
    assert type(actual) is type(expected), f"{path}: {type(actual).__name__} != {type(expected).__name__}"
    assert actual.dtype == expected.dtype, f"{path}: {actual.dtype} != {expected.dtype}"
    if sparse.issparse(expected):
        np.testing.assert_array_equal(actual.toarray(), expected.toarray(), err_msg=path)
    else:
        np.testing.assert_array_equal(np.asarray(actual), np.asarray(expected), err_msg=path)


def _assert_mapping_equal(actual: Mapping[Any, Any], expected: Mapping[Any, Any], *, path: str) -> None:
    assert set(actual) == set(expected), f"{path}: keys differ"
    for key in expected:
        _assert_matrix_equal(actual[key], expected[key], path=f"{path}[{key!r}]")


def assert_adata_equal(actual: AnnData, expected: AnnData) -> None:
    """Assert values, dtypes, categories, sparse formats, order, and aligned shapes."""
    assert actual.shape == expected.shape
    assert actual.obs_names.equals(expected.obs_names), "obs axis order differs"
    assert actual.var_names.equals(expected.var_names), "var axis order differs"
    pd.testing.assert_frame_equal(actual.obs, expected.obs, check_categorical=True, check_dtype=True)
    pd.testing.assert_frame_equal(actual.var, expected.var, check_categorical=True, check_dtype=True)
    _assert_matrix_equal(actual.X, expected.X, path="X")
    _assert_mapping_equal(actual.layers, expected.layers, path="layers")
    _assert_mapping_equal(actual.obsm, expected.obsm, path="obsm")
    _assert_mapping_equal(actual.varm, expected.varm, path="varm")
    _assert_mapping_equal(actual.obsp, expected.obsp, path="obsp")
    _assert_mapping_equal(actual.varp, expected.varp, path="varp")
    assert actual.raw is not None and expected.raw is not None
    pd.testing.assert_frame_equal(actual.raw.var, expected.raw.var, check_categorical=True, check_dtype=True)
    _assert_matrix_equal(actual.raw.X, expected.raw.X, path="raw.X")
    assert actual.obsm["X_reference"].shape == (actual.n_obs, 2)
    assert actual.varm["reference_loadings"].shape == (actual.n_vars, 2)
    assert actual.obsp["connectivities"].shape == (actual.n_obs, actual.n_obs)
    assert actual.varp["feature_similarity"].shape == (actual.n_vars, actual.n_vars)
    assert actual.uns["fixture"] == expected.uns["fixture"]


def _run_grouped_pipeline(source: AnnData) -> tuple[Any, AnnData]:
    grouped = source.ap.group_by(obs="cell_type")
    result = grouped.mutate(obs={"within_type": ap.row_number()}).arrange(obs=ap.desc("quality_score"))
    return result, result.ungroup()


def _assert_independent_ownership(actual: AnnData, source: AnnData) -> None:
    assert actual is not source
    assert not np.shares_memory(actual.X.data, source.X.data)
    assert not np.shares_memory(actual.layers["log1p"].data, source.layers["log1p"].data)
    assert actual.raw is not None and source.raw is not None
    assert not np.shares_memory(actual.raw.X.data, source.raw.X.data)
    assert not np.shares_memory(actual.obsm["X_reference"], source.obsm["X_reference"])
    assert not np.shares_memory(actual.varm["reference_loadings"], source.varm["reference_loadings"])
    assert not np.shares_memory(actual.obsp["connectivities"].data, source.obsp["connectivities"].data)
    assert not np.shares_memory(actual.varp["feature_similarity"].data, source.varp["feature_similarity"].data)


def test_generated_fixture_contract_and_grouped_positional_identity() -> None:
    source = make_representative_adata()
    grouped, result = _run_grouped_pipeline(source)

    assert grouped.group_vars() == ["cell_type"]
    assert grouped.group_keys()["cell_type"].tolist() == ["T cell", "B cell", "Monocyte"]
    assert result.obs["fixture_row_id"].to_numpy().tolist() == EXPECTED_GROUPED_POSITIONS.tolist()
    assert result.obs["within_type"].tolist() == [1, 3, 2, 2, 1, 3, 2, 1]
    expected = source[EXPECTED_GROUPED_POSITIONS, :].copy()
    expected.obs["within_type"] = [1, 3, 2, 2, 1, 3, 2, 1]
    assert_adata_equal(result, expected)
    _assert_independent_ownership(result, source)

    result.obs.iloc[0, result.obs.columns.get_loc("quality_score")] = -1
    result.X.data[0] = -1
    assert source.obs["quality_score"].tolist() == [8, 6, 5, 3, 9, 7, 10, 4]
    assert source.X.toarray()[0].tolist() == [8, 0, 1, 0, 3, 0]


@pytest.mark.parametrize("storage", ["h5ad", "zarr"])
def test_serialization_round_trip_preserves_pipeline(
    storage: str,
    serialized_fixture: tuple[Path, Path],
) -> None:
    source = make_representative_adata()
    h5ad_path, zarr_path = serialized_fixture
    restored = ad.read_h5ad(h5ad_path) if storage == "h5ad" else ad.read_zarr(zarr_path)

    assert_adata_equal(restored, source)
    expected_grouped, expected = _run_grouped_pipeline(source)
    actual_grouped, actual = _run_grouped_pipeline(restored)
    assert actual_grouped.group_vars() == expected_grouped.group_vars() == ["cell_type"]
    assert_adata_equal(actual, expected)
    _assert_independent_ownership(actual, restored)


def test_concatenation_reconstructs_every_aligned_container() -> None:
    source = make_representative_adata()
    reconstructed = ad.concat(
        [source[:4, :].copy(), source[4:, :].copy()],
        axis="obs",
        join="inner",
        merge="same",
        uns_merge="same",
        pairwise=True,
        index_unique=None,
    )

    assert_adata_equal(reconstructed, source)
    _, expected = _run_grouped_pipeline(source)
    grouped, actual = _run_grouped_pipeline(reconstructed)
    assert grouped.group_vars() == ["cell_type"]
    assert_adata_equal(actual, expected)


def test_scanpy_preprocessing_accepts_annplyr_output() -> None:
    scanpy = pytest.importorskip("scanpy", reason="Scanpy is installed only in the integration environment")
    source = make_representative_adata()
    result = source.ap.filter(obs=ap.col("qc_pass").fill_null(False)).ap.select(x=["GATA3", "MS4A1", "LYZ", "NKG7"])
    source_snapshot = source.copy()

    scanpy.pp.normalize_total(result, target_sum=10_000)
    scanpy.pp.log1p(result)
    scanpy.pp.pca(result, n_comps=2, mask_var=np.ones(result.n_vars, dtype=bool), random_state=0)

    assert sparse.isspmatrix_csr(result.X)
    assert result.obsm["X_pca"].shape == (result.n_obs, 2)
    assert result.varm["PCs"].shape == (result.n_vars, 2)
    assert np.isfinite(result.X.data).all()
    assert np.isfinite(result.obsm["X_pca"]).all()
    assert_adata_equal(source, source_snapshot)
