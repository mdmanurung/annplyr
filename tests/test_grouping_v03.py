from __future__ import annotations

import inspect

import narwhals as nw
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap
from annplyr._errors import AnnplyrError, JoinRelationshipError, SelectionError, UnknownColumnError
from annplyr._grouped import GroupedAnnData
from annplyr._groups import GroupPlan


def _grouped_fixture(*, duplicate_names: bool = False) -> AnnData:
    names = ["x", "x", "y", "x", "z", "x"] if duplicate_names else [f"c{i}" for i in range(6)]
    obs = pd.DataFrame(
        {
            "group": pd.Categorical(
                ["b", "a", "b", None, "a", "b"],
                categories=["a", "b", "unused"],
                ordered=True,
            ),
            "value": [3, 2, 1, 5, 4, 3],
            "id": [10, 11, 12, 13, 14, 15],
            "nullable": pd.array([1, pd.NA, 3, 4, pd.NA, 6], dtype="Int64"),
            "flag": pd.array([True, pd.NA, False, True, False, pd.NA], dtype="boolean"),
            "text": pd.array(["a", "b", pd.NA, "d", "e", "f"], dtype="string"),
        },
        index=names,
    )
    var_names = ["g", "g", "h"] if duplicate_names else ["g0", "g1", "g2"]
    var = pd.DataFrame({"kind": ["rna", "protein", "rna"]}, index=var_names)
    return AnnData(X=np.arange(18, dtype=np.float64).reshape(6, 3), obs=obs, var=var)


def test_group_spec_resolves_once_and_rejects_computed_virtual_and_empty_selection() -> None:
    adata = _grouped_fixture()
    grouped = adata.ap.group_by(obs=ap.starts_with("gro"))
    assert grouped.group_vars() == ["group"]
    assert grouped.spec.columns == ("group",)

    with pytest.raises(SelectionError, match="computed"):
        adata.ap.group_by(obs=ap.col("group"))
    with pytest.raises(SelectionError):
        adata.ap.group_by(obs=nw.col("group"))
    with pytest.raises(SelectionError):
        adata.ap.group_by(obs=ap.obs_names)
    with pytest.raises(SelectionError, match="at least one"):
        adata.ap.group_by(obs=ap.any_of(["missing"]))
    with pytest.raises(SelectionError):
        adata.ap.group_by(obs=ap.col("value") + 1)

    del adata.obs["group"]
    with pytest.raises(UnknownColumnError, match="Stored grouping"):
        grouped.group_keys()


def test_group_plan_is_first_seen_na_inclusive_observed_and_positional() -> None:
    grouped = _grouped_fixture().ap.group_by(obs="group")
    keys = grouped.group_keys()
    assert keys["group"].astype(object).where(keys["group"].notna(), "NA").tolist() == ["b", "a", "NA"]
    assert isinstance(keys["group"].dtype, pd.CategoricalDtype)
    assert keys["group"].cat.categories.tolist() == ["a", "b", "unused"]
    assert grouped.group_data()[".rows"].tolist() == [[0, 2, 5], [1, 4], [3]]
    assert "unused" not in keys["group"].dropna().tolist()

    plan = GroupPlan.build(grouped.ungroup(), grouped.spec)
    assert plan.group_ids.tolist() == [0, 1, 0, 2, 1, 0]
    assert plan.first_positions.tolist() == [0, 1, 3]


@pytest.mark.filterwarnings("ignore:Observation names are not unique")
@pytest.mark.filterwarnings("ignore:Variable names are not unique")
def test_grouped_shape_verbs_use_positions_and_persist_with_duplicate_axis_names() -> None:
    grouped = _grouped_fixture(duplicate_names=True).ap.group_by(obs="group")
    operations = {
        "filter": grouped.filter(obs=ap.row_number() == 1),
        "arrange": grouped.arrange(obs="value"),
        "distinct": grouped.distinct(obs="value"),
        "slice": grouped.slice(1),
        "head": grouped.slice_head(1),
        "tail": grouped.slice_tail(1),
        "min": grouped.slice_min("value", 1),
        "max": grouped.slice_max("value", 1),
        "sample": grouped.slice_sample(1, random_state=7),
    }
    assert all(isinstance(result, GroupedAnnData) for result in operations.values())
    assert operations["filter"].ungroup().X[:, 0].tolist() == [0, 3, 9]
    assert operations["arrange"].ungroup().X[:, 0].tolist() == [6, 0, 15, 3, 12, 9]
    assert operations["slice"].ungroup().X[:, 0].tolist() == [6, 12]
    assert operations["head"].ungroup().X[:, 0].tolist() == [0, 3, 9]
    assert operations["tail"].ungroup().X[:, 0].tolist() == [15, 12, 9]
    assert operations["min"].ungroup().X[:, 0].tolist() == [6, 3, 9]
    assert operations["max"].ungroup().X[:, 0].tolist() == [0, 12, 9]
    assert all(result.group_vars() == ["group"] for result in operations.values())


def test_grouped_select_transmute_rename_and_relocate_retain_or_update_keys() -> None:
    grouped = _grouped_fixture().ap.group_by(obs=["group", "id"])
    selected = grouped.select(obs="value")
    assert selected.ungroup().obs.columns.tolist() == ["group", "id", "value"]

    transmuted = grouped.transmute(obs={"double": ap.col("value") * 2})
    assert transmuted.ungroup().obs.columns.tolist() == ["group", "id", "double"]
    assert transmuted.group_vars() == ["group", "id"]

    renamed = grouped.rename(obs={"batch": "group"})
    assert renamed.group_vars() == ["batch", "id"]
    renamed_with = renamed.rename_with(str.upper, obs=ap.all_of(["batch"]))
    assert renamed_with.group_vars() == ["BATCH", "id"]
    relocated = renamed_with.relocate(obs="id", after="value")
    assert relocated.group_vars() == ["BATCH", "id"]


def test_grouped_inplace_identity_and_key_mutation_use_old_plan_then_regroup() -> None:
    grouped = _grouped_fixture().ap.group_by(obs="group")
    underlying = grouped.ungroup()
    result = grouped.mutate(
        obs={
            "within": ap.row_number(),
            "group": ap.if_else(ap.row_number() == 1, "first", "rest"),
        },
        inplace=True,
    )
    assert result is grouped
    assert result.ungroup() is underlying
    assert underlying.obs["within"].tolist() == [1, 1, 2, 1, 2, 3]
    assert result.group_keys()["group"].tolist() == ["first", "rest"]

    renamed = grouped.rename(obs={"partition": "group"}, inplace=True)
    assert renamed is grouped
    assert grouped.group_vars() == ["partition"]
    relocated = grouped.relocate(obs="partition", after="text", inplace=True)
    assert relocated is grouped
    counted = grouped.add_count(sort=True, inplace=True)
    assert counted is grouped
    tallied = grouped.add_tally(inplace=True)
    assert tallied is grouped


def test_grouped_mutate_scatter_preserves_extension_and_categorical_dtypes() -> None:
    grouped = _grouped_fixture().ap.group_by(obs="group")
    result = grouped.mutate(
        obs={
            "nullable_copy": ap.col("nullable"),
            "flag_copy": ap.col("flag"),
            "text_copy": ap.col("text"),
            "category_copy": ap.col("group"),
        }
    ).ungroup()
    assert str(result.obs["nullable_copy"].dtype) == "Int64"
    assert str(result.obs["flag_copy"].dtype) == "boolean"
    assert str(result.obs["text_copy"].dtype).startswith("string")
    assert isinstance(result.obs["category_copy"].dtype, pd.CategoricalDtype)


def test_grouped_summary_count_and_tally_share_plan_order_and_key_dtypes() -> None:
    grouped = _grouped_fixture().ap.group_by(obs="group")
    summary = grouped.summarize(obs={"mean": ap.mean("value")})
    count = grouped.count()
    tally = grouped.tally(wt="value")
    assert summary["group"].astype(object).where(summary["group"].notna(), "NA").tolist() == ["b", "a", "NA"]
    assert count["n"].tolist() == [3, 2, 1]
    assert tally["n"].tolist() == [7, 6, 5]
    assert isinstance(summary["group"].dtype, pd.CategoricalDtype)
    assert isinstance(count["group"].dtype, pd.CategoricalDtype)


def test_all_six_grouped_joins_execute_globally_and_update_suffixed_key() -> None:
    adata = _grouped_fixture()
    grouped = adata.ap.group_by(obs="group")
    right = pd.DataFrame(
        {
            "id": [12, 10, 15, 11, 14, 13],
            "group": ["r12", "r10", "r15", "r11", "r14", "r13"],
            "label": np.arange(6),
        }
    )
    left = grouped.left_join(right, by="id", suffixes=("_left", "_right"))
    inner = grouped.inner_join(right, by="id")
    ordered = grouped.right_join(right, by="id")
    full = grouped.full_join(right, by="id")
    semi = grouped.semi_join(right, by="id")
    anti = grouped.anti_join(right.iloc[:2], by="id")
    assert left.group_vars() == ["group_left"]
    assert ordered.ungroup().obs["id"].tolist() == [12, 10, 15, 11, 14, 13]
    assert inner.ungroup().n_obs == full.ungroup().n_obs == semi.ungroup().n_obs == 6
    assert anti.ungroup().obs["id"].tolist() == [11, 13, 14, 15]
    assert all(isinstance(value, GroupedAnnData) for value in [left, inner, ordered, full, semi, anti])

    right_only = pd.concat([right, pd.DataFrame({"id": [99], "group": ["new"], "label": [99]})], ignore_index=True)
    with pytest.raises(JoinRelationshipError, match="add axis records"):
        grouped.right_join(right_only, by="id", unmatched="drop")
    with pytest.raises(JoinRelationshipError, match="add axis records"):
        grouped.full_join(right_only, by="id", unmatched="drop")


def test_grouped_var_axis_and_non_group_axis_join_remain_positional() -> None:
    adata = _grouped_fixture(duplicate_names=True)
    adata.var["gene_id"] = [0, 1, 2]
    grouped_var = adata.ap.group_by(var="kind")
    sliced = grouped_var.slice_head(1)
    assert sliced.ungroup().X[0, :].tolist() == [0, 1]
    assert sliced.group_vars() == ["kind"]

    grouped_obs = adata.ap.group_by(obs="group")
    right = pd.DataFrame({"gene_id": [2, 0, 1], "label": ["c", "a", "b"]})
    joined = grouped_obs.right_join(right, by="gene_id", axis="var")
    assert joined.ungroup().X[0, :].tolist() == [2, 0, 1]
    assert joined.group_vars() == ["group"]


def test_grouped_shape_verbs_have_explicit_copy_true_signatures() -> None:
    methods = [
        "filter",
        "select",
        "arrange",
        "distinct",
        "slice",
        "slice_head",
        "slice_tail",
        "slice_min",
        "slice_max",
        "slice_sample",
        "left_join",
        "inner_join",
        "right_join",
        "full_join",
        "semi_join",
        "anti_join",
    ]
    for method in methods:
        signature = inspect.signature(getattr(GroupedAnnData, method))
        assert signature.parameters["copy"].default is True
        assert not any(parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in signature.parameters.values())

    assert inspect.signature(GroupedAnnData.left_join).parameters["unmatched"].default == "drop"
    assert inspect.signature(GroupedAnnData.right_join).parameters["unmatched"].default == "error"


def test_grouped_pipe_passes_wrapper_and_supports_keyword_injection() -> None:
    grouped = _grouped_fixture().ap.group_by(obs="group")
    assert grouped.pipe(lambda value: value) is grouped
    assert grouped.pipe((lambda *, data: data, "data")) is grouped
    with pytest.raises(AnnplyrError, match="both the pipe target"):
        grouped.pipe((lambda *, data: data, "data"), data=grouped)


def test_grouped_matrix_methods_expose_budget_and_reject_before_read(monkeypatch) -> None:
    import annplyr._verbs as verbs

    methods = ["filter", "arrange", "distinct", "mutate", "transmute", "summarize", "summarise"]
    for method in methods:
        assert "max_matrix_values" in inspect.signature(getattr(GroupedAnnData, method)).parameters

    class Sentinel:
        reads = 0
        shape = (6, 3)
        names = ("a", "b", "c")
        dtypes = (np.dtype("float64"),) * 3
        schema = pd.DataFrame({name: pd.Series([], dtype="float64") for name in names})

        def read(self, rows, columns):
            self.reads += 1
            return pd.DataFrame(np.zeros((len(rows), len(columns))))

    sentinel = Sentinel()
    monkeypatch.setattr(verbs, "source_adapter", lambda *args, **kwargs: sentinel)
    grouped = _grouped_fixture().ap.group_by(obs="group")
    with pytest.raises(AnnplyrError, match="materialize 6 matrix values"):
        grouped.filter(x=ap.col("a") > 0, max_matrix_values=5)
    assert sentinel.reads == 0


def test_high_cardinality_plan_has_one_nonempty_group_per_first_seen_key() -> None:
    size = 50_000
    adata = AnnData(
        X=np.zeros((size, 1), dtype=np.float32),
        obs=pd.DataFrame({"key": pd.array(np.arange(size), dtype="Int64")}),
    )
    grouped = adata.ap.group_by(obs="key")
    plan = grouped._plan()
    assert len(plan.positions) == size
    assert plan.first_positions.tolist() == list(range(size))
    assert all(len(positions) == 1 for positions in plan.positions)
