from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap
from annplyr._errors import SelectionError, UnknownColumnError
from annplyr.data import X


def _aligned_adata(*, duplicate_names: bool = False) -> AnnData:
    obs_names = ["cell", "cell", "other", "last"] if duplicate_names else ["c0", "c1", "c2", "c3"]
    var_names = ["gene", "gene", "other"] if duplicate_names else ["g0", "g1", "g2"]
    x = np.arange(12, dtype=np.float64).reshape(4, 3)
    adata = AnnData(
        X=x,
        obs=pd.DataFrame(
            {"key": ["b", "a", "b", "a"], "value": [2, 1, 2, 3]},
            index=obs_names,
        ),
        var=pd.DataFrame({"kind": ["b", "a", "b"]}, index=var_names),
    )
    adata.layers["counts"] = x + 100
    adata.obsm["embedding"] = np.arange(8, dtype=np.float64).reshape(4, 2)
    adata.varm["loadings"] = np.arange(6, dtype=np.float64).reshape(3, 2)
    adata.obsp["neighbors"] = np.arange(16, dtype=np.float64).reshape(4, 4)
    adata.varp["similarity"] = np.arange(9, dtype=np.float64).reshape(3, 3)
    adata.raw = adata.copy()
    return adata


def test_shape_copy_is_independent_across_every_aligned_container() -> None:
    source = _aligned_adata()
    result = source.ap.slice([0, 1, 2], axis="obs")

    result.X[0, 0] = -1
    result.layers["counts"][0, 0] = -2
    result.obsm["embedding"][0, 0] = -3
    result.varm["loadings"][0, 0] = -4
    result.obsp["neighbors"][0, 0] = -5
    result.varp["similarity"][0, 0] = -6
    assert result.raw is not None
    result.raw.X[0, 0] = -7

    assert source.X[0, 0] == 0
    assert source.layers["counts"][0, 0] == 100
    assert source.obsm["embedding"][0, 0] == 0
    assert source.varm["loadings"][0, 0] == 0
    assert source.obsp["neighbors"][0, 0] == 0
    assert source.varp["similarity"][0, 0] == 0
    assert source.raw is not None and source.raw.X[0, 0] == 0


def test_subsetting_verbs_do_not_mutate_the_input() -> None:
    source = _aligned_adata()
    before = source.copy()
    result = source.ap.filter(ap.col("value") > 1)
    assert result.obs_names.tolist() == ["c0", "c2", "c3"]
    assert not result.is_view
    np.testing.assert_array_equal(source.X, before.X)
    pd.testing.assert_frame_equal(source.obs, before.obs)


def test_backed_subsetting_returns_the_view_without_loading(tmp_path) -> None:
    path = tmp_path / "backed.h5ad"
    _aligned_adata().write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        for result in (
            backed.ap.slice([3, 1], axis="obs"),
            backed.ap.filter(ap.col("value") > 1),
            backed.ap.select(X.g1),
        ):
            assert result.isbacked
            assert result.is_view
        assert backed.ap.slice([3, 1], axis="obs").obs_names.tolist() == ["c3", "c1"]
        assert backed.ap.filter(ap.col("value") > 1).obs_names.tolist() == ["c0", "c2", "c3"]
    finally:
        backed.file.close()


def test_same_shape_verbs_share_matrices_and_copy_metadata() -> None:
    adata = _aligned_adata()
    result = adata.ap.mutate({"twice": ap.col("value") * 2})
    assert result is not adata
    assert result.X is adata.X
    assert result.layers["counts"] is adata.layers["counts"]
    assert result.obsm["embedding"] is adata.obsm["embedding"]
    assert result.raw is not None and adata.raw is not None
    assert result.raw.X is adata.raw.X
    assert "twice" in result.obs
    assert "twice" not in adata.obs
    result.obs.loc["c0", "value"] = -1
    assert adata.obs.loc["c0", "value"] == 2


def test_same_shape_verbs_return_new_objects_and_validate_before_writing() -> None:
    adata = _aligned_adata()
    result = adata.ap.rename(renamed="value")
    assert result is not adata
    assert result.obs.columns.tolist() == ["key", "renamed"]
    assert adata.obs.columns.tolist() == ["key", "value"]

    snapshot = adata.obs.copy()
    with pytest.raises(UnknownColumnError):
        adata.ap.rename(ok="key", bad=ap.col("var:missing"))
    pd.testing.assert_frame_equal(adata.obs, snapshot)

    with pytest.raises(SelectionError):
        adata.ap.mutate({"would_write": ap.lit(1), "bad": ap.col("missing")})
    assert "would_write" not in adata.obs


def test_relocate_x_and_sorted_counts_leave_the_input_untouched() -> None:
    adata = _aligned_adata()
    expected_x = np.asarray(adata.X)[:, [2, 0, 1]].copy()
    expected_layer = np.asarray(adata.layers["counts"])[:, [2, 0, 1]].copy()
    result = adata.ap.relocate("x:g2")
    assert result is not adata
    assert result.var_names.tolist() == ["g2", "g0", "g1"]
    np.testing.assert_array_equal(result.X, expected_x)
    np.testing.assert_array_equal(result.layers["counts"], expected_layer)
    assert adata.var_names.tolist() == ["g0", "g1", "g2"]

    adata.obs["key"] = ["b", "a", "b", "b"]
    counted = adata.ap.add_count("key", sort=True)
    assert counted is not adata
    assert counted.obs_names.tolist() == ["c0", "c2", "c3", "c1"]
    assert counted.obs["n"].tolist() == [3, 3, 3, 1]
    assert adata.obs_names.tolist() == ["c0", "c1", "c2", "c3"]
    assert "n" not in adata.obs


@pytest.mark.filterwarnings("ignore:Observation names are not unique")
@pytest.mark.filterwarnings("ignore:Variable names are not unique")
def test_positional_verbs_and_joins_are_unambiguous_with_duplicate_axis_names() -> None:
    adata = _aligned_adata(duplicate_names=True)

    assert adata.ap.filter(ap.col("value") >= 2).X[:, 0].tolist() == [0, 6, 9]
    assert adata.ap.arrange(ap.col("value")).X[:, 0].tolist() == [3, 0, 6, 9]
    assert adata.ap.distinct("value").X[:, 0].tolist() == [0, 3, 9]
    assert adata.ap.slice(1, 0).X[:, 0].tolist() == [3, 0]
    assert adata.ap.select("value").X[:, 0].tolist() == [0, 3, 6, 9]

    other = pd.DataFrame({"key": ["a", "b"], "label": [10, 20]})
    left = adata.ap.left_join(other, by="key")
    inner = adata.ap.inner_join(other, by="key")
    right = adata.ap.right_join(other, by="key")
    full = adata.ap.full_join(other, by="key")
    semi = adata.ap.semi_join(other, by="key")
    anti = adata.ap.anti_join(pd.DataFrame({"key": ["a"]}), by="key")
    assert left.X[:, 0].tolist() == [0, 3, 6, 9]
    assert inner.X[:, 0].tolist() == [0, 3, 6, 9]
    assert right.X[:, 0].tolist() == [3, 9, 0, 6]
    assert full.X[:, 0].tolist() == [0, 3, 6, 9]
    assert semi.X[:, 0].tolist() == [0, 3, 6, 9]
    assert anti.X[:, 0].tolist() == [0, 6]

    assert adata.ap.filter(ap.col("var:kind") == "b").X[0, :].tolist() == [0, 2]
    assert adata.ap.arrange(ap.col("var:kind")).X[0, :].tolist() == [1, 0, 2]
    assert adata.ap.distinct("var:kind").X[0, :].tolist() == [0, 1]
    assert adata.ap.slice(1, 0, axis="var").X[0, :].tolist() == [1, 0]
    assert adata.ap.select(ap.col("var:kind")).X[0, :].tolist() == [0, 1, 2]

    var_other = pd.DataFrame({"kind": ["a", "b"], "label": [10, 20]})
    assert adata.ap.left_join(var_other, by="kind", axis="var").X[0, :].tolist() == [0, 1, 2]


def test_transmute_returns_a_new_object_that_shares_matrices() -> None:
    source = _aligned_adata()
    result = source.ap.transmute({"twice": ap.col("value") * 2})
    assert result is not source
    assert result.X is source.X
    assert result.obs.columns.tolist() == ["twice"]
    assert source.obs.columns.tolist() == ["key", "value"]


def test_same_shape_exported_utilities_return_new_objects() -> None:
    adata = _aligned_adata()
    renamed = ap.rename_obs_names(adata, lambda name: f"x_{name}")
    prefixed = ap.add_name_prefix(adata, "v", axis="var")
    painted = ap.store_palette(adata, "key", ["#000", "#fff"])
    for result in (renamed, prefixed, painted):
        assert result is not adata
        assert result.X is adata.X
    assert renamed.obs_names.tolist() == ["x_c0", "x_c1", "x_c2", "x_c3"]
    assert prefixed.var_names.tolist() == ["v_g0", "v_g1", "v_g2"]
    assert "key_colors" in painted.uns
    assert adata.obs_names.tolist() == ["c0", "c1", "c2", "c3"]
    assert adata.var_names.tolist() == ["g0", "g1", "g2"]
    assert "key_colors" not in adata.uns
