from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

import annplyr as ap
from annplyr._frames import source_frame


def test_left_join_enriches_obs_without_reordering_or_duplicating(dense_adata: AnnData) -> None:
    metadata = pd.DataFrame(
        {
            "batch": ["B", "A"],
            "treatment": ["drug", "vehicle"],
        }
    )

    joined = dense_adata.ap.left_join(metadata, by="batch", axis="obs")

    assert joined.obs_names.tolist() == dense_adata.obs_names.tolist()
    assert joined.obs["treatment"].tolist() == ["vehicle", "vehicle", "drug", "drug", "vehicle"]
    assert "treatment" not in dense_adata.obs


def test_join_relationship_and_unmatched_axis_records_are_checked(dense_adata: AnnData) -> None:
    duplicated = pd.DataFrame({"batch": ["A", "A"], "label": ["x", "y"]})
    with pytest.raises(ap.JoinRelationshipError, match="multiple"):
        dense_adata.ap.left_join(duplicated, by="batch", axis="obs")

    extra = pd.DataFrame({"batch": ["A", "B", "C"], "label": ["x", "y", "z"]})
    with pytest.raises(ap.JoinRelationshipError, match="add axis records"):
        dense_adata.ap.right_join(extra, by="batch", axis="obs")

    with pytest.raises(ap.JoinRelationshipError, match="unmatched"):
        dense_adata.ap.inner_join(pd.DataFrame({"batch": ["A"], "label": ["x"]}), by="batch", unmatched="error")

    with pytest.raises(ap.JoinRelationshipError, match="one-to-one"):
        dense_adata.ap.left_join(
            pd.DataFrame({"batch": ["A", "B"], "label": ["x", "y"]}),
            by="batch",
            relationship="one-to-one",
        )


def test_join_na_matches_never_does_not_enrich_missing_keys(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["donor"] = ["d1", None, "d2", None, "d1"]
    other = pd.DataFrame({"donor": ["d1", None], "site": ["blood", "missing"]})

    joined = adata.ap.left_join(other, by="donor", na_matches="never")

    assert joined.obs.loc["c0", "site"] == "blood"
    assert pd.isna(joined.obs.loc["c1", "site"])
    assert pd.isna(joined.obs.loc["c3", "site"])


def test_semi_join_and_anti_join_subset_anndata_axes(dense_adata: AnnData) -> None:
    keep_obs = pd.DataFrame({"batch": ["A"]})
    assert dense_adata.ap.semi_join(keep_obs, by="batch").obs_names.tolist() == ["c0", "c1", "c4"]
    assert dense_adata.ap.anti_join(keep_obs, by="batch").obs_names.tolist() == ["c2", "c3"]

    keep_var = pd.DataFrame({"feature_type": ["protein"]})
    assert dense_adata.ap.semi_join(keep_var, by="feature_type", axis="var").var_names.tolist() == ["g2"]


def test_pivot_longer_and_pivot_wider_round_trip_selected_values(dense_adata: AnnData) -> None:
    long = dense_adata.ap.pivot_longer(obs=["batch"], x=["g0", "g3"])

    assert long.columns.tolist() == ["obs_name", "batch", "name", "value"]
    assert len(long) == dense_adata.n_obs * 2
    assert long.loc[(long["obs_name"] == "c3") & (long["name"] == "g3"), "value"].item() == 8.0

    wide = ap.pivot_wider(long, id_cols=["obs_name", "batch"], names_from="name", values_from="value")
    assert wide.columns.tolist() == ["obs_name", "batch", "g0", "g3"]
    assert wide.loc[wide["obs_name"] == "c4", "g3"].item() == 4.0


def test_pivot_longer_rejects_cross_source_duplicate_names(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["g0"] = ["meta"] * adata.n_obs

    with pytest.raises(ap.NameRepairError, match="Duplicate"):
        adata.ap.pivot_longer(obs=["g0"], x=["g0"])


def test_to_tidy_requires_explicit_features_by_default(dense_adata: AnnData) -> None:
    with pytest.raises(ap.SelectionError, match="explicit x"):
        dense_adata.ap.to_tidy(obs=["batch"])


def test_to_tidy_rejects_reserved_name_collisions(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.var_names = ["obs_name", "feature", "value", "g3"]

    with pytest.raises(ap.NameRepairError, match="reserved"):
        adata.ap.to_tidy(obs=["batch"], x=["obs_name"])

    adata = dense_adata.copy()
    adata.obs["feature"] = ["real"] * adata.n_obs
    with pytest.raises(ap.NameRepairError, match="reserved"):
        adata.ap.to_tidy(obs=["feature"], x=["g0"])


def test_nest_by_and_unnest_return_pandas_tables(dense_adata: AnnData) -> None:
    nested = dense_adata.ap.nest_by(by="batch", obs=["cell_type", "score"])

    assert nested.columns.tolist() == ["batch", "data"]
    assert nested["batch"].tolist() == ["A", "B"]
    assert nested.loc[0, "data"].columns.tolist() == ["cell_type", "score"]

    unnested = ap.unnest(nested, "data")
    assert unnested.columns.tolist() == ["batch", "cell_type", "score"]
    assert unnested["batch"].tolist() == ["A", "A", "A", "B", "B"]


def test_source_registry_supports_pairwise_and_controlled_uns_sources(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obsp["connectivities"] = np.eye(adata.n_obs)
    adata.varp["correlations"] = np.eye(adata.n_vars)
    adata.uns["qc"] = pd.DataFrame({"metric": ["cells"], "value": [adata.n_obs]})

    obsp = source_frame(adata, "obsp", key="connectivities")
    varp = source_frame(adata, "varp", key="correlations")
    uns = source_frame(adata, "uns", key="qc")

    assert obsp.index.tolist() == adata.obs_names.tolist()
    assert obsp.columns.tolist() == adata.obs_names.tolist()
    assert varp.index.tolist() == adata.var_names.tolist()
    assert uns["metric"].tolist() == ["cells"]

    adata.uns["scalar_mapping"] = {"a": 1, "b": 2}
    with pytest.raises(ap.UnknownSourceError, match="DataFrame"):
        source_frame(adata, "uns", key="scalar_mapping")


def test_sparse_frames_preserve_sparse_dtypes_for_matrix_sources(sparse_adata: AnnData) -> None:
    frame = source_frame(sparse_adata, "x")

    assert all(isinstance(dtype, pd.SparseDtype) for dtype in frame.dtypes)
    filtered = sparse_adata.ap.filter(x=ap.col("g3") >= 40, layer="counts")
    assert filtered.obs_names.tolist() == ["c3", "c4"]


def test_backed_mutating_verbs_raise_typed_error(tmp_path, dense_adata: AnnData) -> None:
    path = tmp_path / "adata.h5ad"
    dense_adata.write_h5ad(path)
    backed = read_h5ad(path, backed="r")

    try:
        counted = backed.ap.count(by="batch")
        pd.testing.assert_frame_equal(
            counted.reset_index(drop=True),
            pd.DataFrame({"batch": ["A", "B"], "n": [3, 2]}),
        )

        with pytest.raises(ap.AnnplyrError, match="backed"):
            backed.ap.mutate(obs={"double_score": ap.col("score") * 2})
    finally:
        backed.file.close()


def test_pairwise_containers_remain_aligned_after_axis_reordering(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obsp["connectivities"] = np.arange(adata.n_obs * adata.n_obs).reshape(adata.n_obs, adata.n_obs)
    adata.varp["correlations"] = np.arange(adata.n_vars * adata.n_vars).reshape(adata.n_vars, adata.n_vars)

    arranged = adata.ap.arrange(obs=ap.desc("score"), var="length")

    obs_order = [adata.obs_names.get_loc(name) for name in arranged.obs_names]
    var_order = [adata.var_names.get_loc(name) for name in arranged.var_names]
    np.testing.assert_array_equal(
        arranged.obsp["connectivities"], adata.obsp["connectivities"][np.ix_(obs_order, obs_order)]
    )
    np.testing.assert_array_equal(
        arranged.varp["correlations"], adata.varp["correlations"][np.ix_(var_order, var_order)]
    )


def test_expression_verbs_do_not_call_anndata_to_df(monkeypatch, dense_adata: AnnData) -> None:
    def fail_to_df(*args, **kwargs):
        raise AssertionError("AnnData.to_df should not be used by internal expression verbs")

    monkeypatch.setattr(dense_adata, "to_df", fail_to_df)

    filtered = dense_adata.ap.filter(x=ap.col("g0") > 1)
    assert filtered.obs_names.tolist() == ["c1", "c3", "c4"]
    mutated = dense_adata.ap.mutate(x={"g3_counts": ap.col("g3")})
    assert mutated.obs["g3_counts"].tolist() == [2.0, 1.0, 0.0, 8.0, 4.0]


def test_metadata_assignment_preserves_common_extension_dtypes(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["cluster"] = pd.Categorical(["a", "a", "b", "b", "a"], categories=["a", "b"], ordered=True)
    adata.obs["nullable_n"] = pd.Series([1, None, 2, None, 3], index=adata.obs_names, dtype="Int64")
    adata.obs["flag"] = pd.Series([True, None, False, True, None], index=adata.obs_names, dtype="boolean")
    adata.obs["label"] = pd.Series(["x", "y", None, "z", "x"], index=adata.obs_names, dtype="string")

    mutated = adata.ap.mutate(
        obs={
            "cluster_copy": ap.col("cluster"),
            "nullable_copy": ap.col("nullable_n"),
            "flag_copy": ap.col("flag"),
            "label_copy": ap.col("label"),
        }
    )

    assert isinstance(mutated.obs["cluster_copy"].dtype, pd.CategoricalDtype)
    assert str(mutated.obs["nullable_copy"].dtype) == "Int64"
    assert str(mutated.obs["flag_copy"].dtype) == "boolean"
    assert str(mutated.obs["label_copy"].dtype) == "string"
