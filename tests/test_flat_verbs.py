"""The flat verb surfaces: one payload, containers named per expression."""

from __future__ import annotations

from typing import Any

import anndata
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap
from annplyr._errors import IncompatibleAxisError, UnknownColumnError, UnknownSourceError


@pytest.fixture
def multi() -> AnnData:
    """An object with every container the qualifiers can name."""
    rng = np.random.default_rng(0)
    adata = AnnData(
        X=rng.normal(size=(8, 4)).astype("float32"),
        obs=pd.DataFrame(
            {
                "n_genes": [100, 200, 300, 400, 500, 600, 700, 800],
                "group": ["a", "a", "b", "b", "a", "b", "a", "b"],
            },
            index=[f"c{i}" for i in range(8)],
        ),
        var=pd.DataFrame(index=["MS4A1", "CD79A", "NKG7", "LYZ"]),
    )
    adata.layers["counts"] = (adata.X > 0).astype("float32")
    adata.raw = adata
    adata.obsm["X_pca"] = rng.normal(size=(8, 3))
    adata.varm["loadings"] = rng.normal(size=(4, 2))
    return adata


class TestFilter:
    def test_bare_name_reads_obs(self, multi: AnnData) -> None:
        result = multi.ap.filter(ap.col("n_genes") >= 500)

        assert result.obs_names.tolist() == ["c4", "c5", "c6", "c7"]

    def test_several_predicates_are_anded(self, multi: AnnData) -> None:
        result = multi.ap.filter(ap.col("n_genes") >= 400, ap.col("group") == "a")

        assert result.obs_names.tolist() == ["c4", "c6"]

    def test_a_tuple_is_accepted_as_sugar(self, multi: AnnData) -> None:
        assert multi.ap.filter((ap.col("n_genes") >= 400,)).obs_names.tolist() == ["c3", "c4", "c5", "c6", "c7"]

    def test_x_and_a_layer_are_readable_in_one_call(self, multi: AnnData) -> None:
        # the old layer= switch made this impossible
        result = multi.ap.filter(ap.col("x:MS4A1") > 0, ap.col("counts:LYZ") > 0)

        assert result.n_obs < multi.n_obs

    def test_obs_raw_and_obsm_together(self, multi: AnnData) -> None:
        result = multi.ap.filter(ap.col("n_genes") >= 400, ap.col("raw:LYZ") < 0.5, ap.col("X_pca:0") >= 0)

        assert set(result.obs_names).issubset(set(multi.obs_names))

    def test_a_matrix_column_needs_a_qualifier(self, multi: AnnData) -> None:
        with pytest.raises(UnknownColumnError):
            multi.ap.filter(ap.col("MS4A1") > 0)

    def test_var_qualifier_subsets_features(self, multi: AnnData) -> None:
        result = multi.ap.filter(ap.col("var:__annplyr_var_names__").str.contains("MS4A1"))

        assert result.var_names.tolist() == ["MS4A1"]

    def test_cross_container_predicate_is_rejected(self, multi: AnnData) -> None:
        with pytest.raises(IncompatibleAxisError, match="cannot mix columns from obs, raw"):
            multi.ap.filter(ap.col("n_genes") > 100 + ap.col("raw:LYZ"))

    def test_a_mistyped_qualifier_is_reported(self, multi: AnnData) -> None:
        with pytest.raises(UnknownSourceError, match="not a known container"):
            multi.ap.filter(ap.col("nosuch:LYZ") > 0)

    def test_source_points_an_opaque_expression_at_a_matrix(self, multi: AnnData) -> None:
        import narwhals as nw

        positives = int((multi.X[:, 0] > 0).sum())

        result = multi.ap.filter(nw.col("MS4A1") > 0, source="x")

        assert result.n_obs == positives

    def test_budget_is_enforced_before_any_read(self, multi: AnnData) -> None:
        from annplyr._errors import AnnplyrError

        with pytest.raises(AnnplyrError, match="exceeds max_matrix_values=7"):
            multi.ap.filter(ap.col("raw:LYZ") > 0, max_matrix_values=7)


class TestArrange:
    def test_bare_name_sorts_obs(self, multi: AnnData) -> None:
        result = multi.ap.arrange(ap.col("n_genes"))

        assert result.obs_names.tolist() == [f"c{i}" for i in range(8)]

    def test_desc_reverses(self, multi: AnnData) -> None:
        result = multi.ap.arrange(ap.desc(ap.col("n_genes")))

        assert result.obs_names.tolist() == [f"c{i}" for i in reversed(range(8))]

    def test_keys_apply_left_to_right(self, multi: AnnData) -> None:
        result = multi.ap.arrange(ap.col("group"), ap.desc(ap.col("n_genes")))

        assert result.obs_names.tolist() == ["c6", "c4", "c1", "c0", "c7", "c5", "c3", "c2"]

    def test_a_matrix_sort_key(self, multi: AnnData) -> None:
        result = multi.ap.arrange(ap.col("raw:LYZ"))

        assert sorted(result.obs_names) == sorted(multi.obs_names)

    def test_a_var_sort_key(self, multi: AnnData) -> None:
        result = multi.ap.arrange(ap.col("var:__annplyr_var_names__"))

        assert result.var_names.tolist() == sorted(multi.var_names)


class TestMutate:
    def test_bare_name_reads_obs_and_writes_obs(self, multi: AnnData) -> None:
        result = multi.ap.mutate({"doubled": ap.col("n_genes") * 2})

        assert result.obs["doubled"].tolist() == [200, 400, 600, 800, 1000, 1200, 1400, 1600]

    def test_into_var_writes_the_feature_axis(self, multi: AnnData) -> None:
        result = multi.ap.mutate(
            {"length": ap.col("var:__annplyr_var_names__").str.len_chars()},
            into="var",
        )

        assert "length" in result.var.columns
        assert "length" not in result.obs.columns

    def test_a_matrix_value_reads_that_container(self, multi: AnnData) -> None:
        result = multi.ap.mutate({"raw_lyz": ap.col("raw:LYZ")})

        assert result.obs["raw_lyz"].tolist() == multi.raw.X[:, 3].tolist()

    def test_assignments_are_sequential(self, multi: AnnData) -> None:
        result = multi.ap.mutate(
            {
                "value": ap.col("raw:LYZ"),
                "doubled": ap.col("value") * 2,
            }
        )

        assert result.obs["doubled"].tolist() == [v * 2 for v in result.obs["value"]]

    def test_a_later_assignment_sees_an_earlier_overwrite(self, multi: AnnData) -> None:
        result = multi.ap.mutate(
            {
                "n_genes": ap.col("n_genes") + 1,
                "after": ap.col("n_genes") * 2,
            }
        )

        assert result.obs["after"].tolist() == [(v + 1) * 2 for v in range(100, 900, 100)]

    def test_a_literal_recycles_across_rows(self, multi: AnnData) -> None:
        result = multi.ap.mutate({"flag": ap.lit(True)})

        assert result.obs["flag"].tolist() == [True] * 8

    def test_independent_assignments_share_one_evaluation(self, multi: AnnData) -> None:
        import annplyr._frames as frames

        calls = 0
        original = frames.nw.from_native

        def counted(*args: Any, **kwargs: Any) -> Any:
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        frames.nw.from_native = counted
        try:
            multi.ap.mutate({"a": ap.col("n_genes") + 1, "b": ap.col("group")})
        finally:
            frames.nw.from_native = original

        assert calls == 1

    def test_containers_may_be_mixed_in_one_call(self, multi: AnnData) -> None:
        result = multi.ap.mutate(
            {
                "from_x": ap.col("x:MS4A1"),
                "from_counts": ap.col("counts:LYZ"),
            }
        )

        assert result.obs["from_x"].tolist() == multi.X[:, 0].tolist()

    def test_extension_dtypes_survive(self) -> None:
        adata = AnnData(
            X=np.zeros((3, 2), dtype="float32"),
            obs=pd.DataFrame(
                {
                    "nullable": pd.array([1, 2, pd.NA], dtype="Int64"),
                    "flag": pd.array([True, pd.NA, False], dtype="boolean"),
                    "label": pd.Categorical(["u", "v", "u"]),
                },
                index=["a", "b", "c"],
            ),
            var=pd.DataFrame(index=["g0", "g1"]),
        )

        result = adata.ap.mutate(
            {
                "nullable_copy": ap.col("nullable"),
                "flag_copy": ap.col("flag"),
                "label_copy": ap.col("label"),
            }
        )

        assert str(result.obs["nullable_copy"].dtype) == "Int64"
        assert str(result.obs["flag_copy"].dtype) == "boolean"
        assert isinstance(result.obs["label_copy"].dtype, pd.CategoricalDtype)
        assert result.obs["nullable_copy"].isna().tolist() == [False, False, True]

    def test_inplace_still_returns_the_same_object(self, multi: AnnData) -> None:
        result = multi.ap.mutate({"flag": ap.lit(1)}, inplace=True)

        assert result is multi

    def test_across_expands_within_the_default_frame(self) -> None:
        adata = AnnData(
            X=np.zeros((2, 2), dtype="float32"),
            obs=pd.DataFrame({"score": [1.0, 2.0], "other": [3.0, 4.0]}),
            var=pd.DataFrame(index=["g0", "g1"]),
        )

        result = adata.ap.mutate(
            {"across": ap.across(ap.starts_with("s"), lambda column: ap.col(column) * 2, names="{col}_x2")}
        )

        assert result.obs["score_x2"].tolist() == [2.0, 4.0]

    def test_cross_container_assignment_is_rejected(self, multi: AnnData) -> None:
        with pytest.raises(IncompatibleAxisError, match="cannot mix columns from obs, raw"):
            multi.ap.mutate({"bad": ap.col("n_genes") + ap.col("raw:LYZ")})

    def test_a_backed_object_is_rejected(self, multi: AnnData, tmp_path: Any) -> None:
        from annplyr._errors import AnnplyrError

        path = tmp_path / "b.h5ad"
        multi.write_h5ad(path)
        backed = anndata.read_h5ad(path, backed="r")
        try:
            with pytest.raises(AnnplyrError, match="backed mode"):
                backed.ap.mutate({"flag": ap.lit(1)})
        finally:
            backed.file.close()

    def test_transmute_keeps_only_the_named_columns(self, multi: AnnData) -> None:
        result = multi.ap.transmute({"only": ap.col("n_genes")})

        assert result.obs.columns.tolist() == ["only"]

    def test_grouped_mutate_retains_its_keys(self, multi: AnnData) -> None:
        result = multi.ap.group_by("group").transmute({"rank": ap.row_number()}).ungroup()

        assert set(result.obs.columns) == {"group", "rank"}


class TestSummarize:
    def test_bare_names_reduce_obs(self, multi: AnnData) -> None:
        summary = multi.ap.summarize({"cells": ap.n()}, by="group")

        assert summary["cells"].tolist() == [4, 4]

    def test_obs_and_raw_mix_in_one_dict(self, multi: AnnData) -> None:
        summary = multi.ap.summarize(
            {"cells": ap.n(), "mean_lyz": ap.mean("raw:LYZ")},
            by="group",
        )

        assert list(summary.columns) == ["group", "cells", "mean_lyz"]
        assert summary["mean_lyz"].notna().all()

    def test_several_containers_at_once(self, multi: AnnData) -> None:
        summary = multi.ap.summarize(
            {
                "n": ap.n(),
                "from_x": ap.mean("x:MS4A1"),
                "from_counts": ap.mean("counts:LYZ"),
                "from_pca": ap.mean("X_pca:0"),
            },
            by="group",
        )

        assert list(summary.columns) == ["group", "n", "from_x", "from_counts", "from_pca"]

    def test_without_by_reduces_to_one_row(self, multi: AnnData) -> None:
        assert len(multi.ap.summarize({"cells": ap.n()})) == 1

    def test_a_var_qualifier_reduces_features(self, multi: AnnData) -> None:
        # grouping the feature axis needs a real var column, since the frame has
        # no virtual column to select
        multi.var["kind"] = ["receptor", "receptor", "signaling", "signaling"]
        summary = multi.ap.summarize(
            {"features": ap.n_distinct("var:__annplyr_var_names__")},
            by="var:kind",
        )

        assert list(summary.columns) == ["kind", "features"]
        assert summary["features"].tolist() == [2, 2]

    def test_mixed_axes_are_rejected(self, multi: AnnData) -> None:
        with pytest.raises(IncompatibleAxisError, match="obs-axis or var-axis"):
            multi.ap.summarize({"a": ap.n(), "b": ap.n_distinct("var:__annplyr_var_names__")})

    def test_a_cross_container_aggregate_is_rejected(self, multi: AnnData) -> None:
        with pytest.raises(IncompatibleAxisError, match="cannot mix columns from obs, raw"):
            multi.ap.summarize({"bad": ap.mean("n_genes") + ap.mean("raw:LYZ")})

    def test_grouped_summarize_uses_the_stored_keys(self, multi: AnnData) -> None:
        grouped = multi.ap.group_by("group").summarize({"cells": ap.n()})

        assert list(grouped.columns) == ["group", "cells"]
        assert grouped["cells"].tolist() == [4, 4]

    def test_budget_is_enforced_before_any_read(self, multi: AnnData) -> None:
        from annplyr._errors import AnnplyrError

        with pytest.raises(AnnplyrError, match="exceeds max_matrix_values=7"):
            multi.ap.summarize({"m": ap.mean("raw:LYZ")}, max_matrix_values=7)

    def test_sparse_dtypes_still_promote(self) -> None:
        from scipy import sparse

        adata = AnnData(
            X=sparse.csr_matrix(np.array([[1, 1, 1, 1], [2, 2, 2, 2]], dtype=np.float32)),
            obs=pd.DataFrame({"label": ["a", "b"]}),
            var=pd.DataFrame(index=["g0", "g1", "g2", "g3"]),
        )

        summary = adata.ap.summarize({"total": ap.sum("x:g0") + ap.sum("x:g1")}, by="label")

        assert summary["total"].tolist() == [2.0, 4.0]


class TestSharedContracts:
    @pytest.mark.parametrize(
        ("call", "expected"),
        [
            (lambda a: a.ap.filter(ap.col("n_genes") >= 500), (4, 4)),
            (lambda a: a.ap.arrange(ap.col("n_genes")), (8, 4)),
            (lambda a: a.ap.mutate({"x": ap.lit(1)}), (8, 4)),
            (lambda a: a.ap.transmute({"x": ap.lit(1)}), (8, 4)),
        ],
    )
    def test_shape_verbs_still_return_anndata(self, multi: AnnData, call: Any, expected: tuple) -> None:
        from anndata import AnnData as AnnDataType

        result = call(multi)

        assert isinstance(result, AnnDataType)
        assert result.shape == expected

    def test_a_flat_verb_does_not_mutate_its_input(self, multi: AnnData) -> None:
        before_obs = multi.obs.copy()
        before_var = multi.var.copy()

        multi.ap.filter(ap.col("n_genes") >= 500)
        multi.ap.mutate({"flag": ap.lit(1)})
        multi.ap.arrange(ap.col("n_genes"))
        multi.ap.summarize({"n": ap.n()})

        pd.testing.assert_frame_equal(multi.obs, before_obs)
        pd.testing.assert_frame_equal(multi.var, before_var)

    def test_a_scalar_value_recycles_across_the_axis(self, multi: AnnData) -> None:
        import narwhals as nw

        result = multi.ap.mutate({"mean": nw.col("n_genes").mean()})

        assert result.obs["mean"].nunique() == 1
