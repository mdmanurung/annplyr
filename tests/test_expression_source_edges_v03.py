from __future__ import annotations

from collections.abc import Callable

import narwhals as nw
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap
from annplyr._errors import (
    AnnplyrError,
    DuplicateNameError,
    SelectionError,
    SizeMismatchError,
    UnknownColumnError,
    UnknownSourceError,
)
from annplyr._expr import AnnplyrExpr, expression_dependencies, to_narwhals
from annplyr._frames import (
    _aligned_assignment_series,
    evaluate_assignments,
    evaluate_filter,
    evaluate_select,
    matrix_frame,
    obs_frame,
    source_frame,
    var_frame,
    with_row_number,
)
from annplyr._sources import (
    RequestPlanner,
    _name_positions,
    _positions,
    adapter_for,
    request_dependencies,
    resolve_columns,
    source_adapter,
)


def _evaluate(frame: pd.DataFrame, expression: AnnplyrExpr) -> pd.Series:
    return evaluate_assignments(frame, {"result": expression})["result"]


def test_expression_operator_proxy_matches_python_arithmetic_and_boolean_semantics() -> None:
    frame = pd.DataFrame({"x": [2, 3], "y": [5, 7], "flag": [True, False]})
    expressions = {
        "radd": 10 + ap.col("x"),
        "rsub": 10 - ap.col("x"),
        "rmul": 10 * ap.col("x"),
        "rdiv": 12 / ap.col("x"),
        "rfloor": 13 // ap.col("x"),
        "rmod": 13 % ap.col("x"),
        "rpow": 2 ** ap.col("x"),
        "and": True & ap.col("flag"),
        "or": False | ap.col("flag"),
        "eq": ap.col("x") == 2,
        "ne": ap.col("x") != 2,
        "lt": ap.col("x") < ap.col("y"),
        "le": ap.col("x") <= 2,
        "gt": ap.col("y") > ap.col("x"),
        "ge": ap.col("x") >= 3,
        "invert": ~ap.col("flag"),
        "neg": -ap.col("x"),
    }
    result = evaluate_assignments(frame, expressions)
    assert result.to_dict("list") == {
        "radd": [12, 13],
        "rsub": [8, 7],
        "rmul": [20, 30],
        "rdiv": [6.0, 4.0],
        "rfloor": [6, 4],
        "rmod": [1, 1],
        "rpow": [4, 8],
        "and": [True, False],
        "or": [True, False],
        "eq": [True, False],
        "ne": [False, True],
        "lt": [True, True],
        "le": [True, False],
        "gt": [True, True],
        "ge": [False, True],
        "invert": [False, True],
        "neg": [-2, -3],
    }

    expression = ap.col("x")
    assert +expression is expression
    with pytest.raises(TypeError, match="Python boolean"):
        bool(expression)


def test_expression_metadata_is_conservative_for_nested_and_opaque_inputs() -> None:
    expression = ap.col("x").clip(lower_bound=ap.col("lower"), upper_bound=10)
    assert expression.dependencies == frozenset({"x", "lower"})
    assert expression.cardinality == "row"
    assert to_narwhals({"x": [expression]}) == {"x": [expression.to_narwhals()]}

    opaque = ap.coalesce(ap.col("x"), nw.col("fallback"))
    assert opaque.dependencies is None
    assert opaque.cardinality == "unknown"
    assert expression_dependencies(ap.desc(expression)) == frozenset({"x", "lower"})
    assert expression_dependencies(object()) == frozenset()


def test_shallow_evaluation_frames_do_not_expose_whole_column_writes() -> None:
    adata = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"value": [1, 2]}, index=["a", "b"]),
        var=pd.DataFrame({"kind": ["rna"]}, index=["g"]),
    )
    obs = obs_frame(adata)
    var = var_frame(adata)
    work = with_row_number(obs)
    obs["value"] = [10, 20]
    var["kind"] = ["protein"]
    work["value"] = [30, 40]
    assert adata.obs["value"].tolist() == [1, 2]
    assert adata.var["kind"].tolist() == ["rna"]


def test_public_expression_helpers_have_value_parity_on_rows_and_aggregates() -> None:
    frame = pd.DataFrame(
        {
            "x": [3, 1, 1, 4],
            "nullable": [1.0, np.nan, 3.0, np.nan],
            "flag": [False, True, False, True],
        }
    )

    aggregate_cases = [
        (ap.n(), 4),
        (ap.n_distinct("x"), 3),
        (ap.mean("x"), 2.25),
        (ap.median("x"), 2.0),
        (ap.sum("x"), 9),
        (ap.min("x"), 1),
        (ap.max("x"), 4),
        (ap.first("x"), 3),
        (ap.last("x"), 4),
        (ap.nth("x", 1), 1),
        (ap.nth("x", -1), 4),
        (ap.nth("x", 99, default=-1), -1),
    ]
    for expression, expected in aggregate_cases:
        result = _evaluate(frame, expression)
        assert len(result) == len(frame)
        assert result.iloc[0] == expected
    assert _evaluate(frame, ap.sd("x")).iloc[0] == pytest.approx(pd.Series(frame["x"]).std())

    assert _evaluate(frame, ap.lead("x", default=-1)).tolist() == [1, 1, 4, -1]
    assert _evaluate(frame, ap.lag("x", default=-1)).tolist() == [-1, 3, 1, 1]
    assert _evaluate(frame, ap.lead("x", -1, default=-1)).tolist() == [-1, 3, 1, 1]
    assert _evaluate(frame, ap.lag("x", -1, default=-1)).tolist() == [1, 1, 4, -1]
    assert _evaluate(frame, ap.lead("x", 0)).tolist() == frame["x"].tolist()
    assert _evaluate(frame, ap.lag("x", 0)).tolist() == frame["x"].tolist()
    assert pd.isna(_evaluate(frame, ap.lead("x")).iloc[-1])
    assert pd.isna(_evaluate(frame, ap.lag("x")).iloc[0])

    row_cases = {
        "coalesce": (ap.coalesce("nullable", "x"), [1.0, 1.0, 3.0, 4.0]),
        "replace_na": (ap.replace_na("nullable", 0), [1.0, 0.0, 3.0, 0.0]),
        "is_na": (ap.is_na("nullable"), [False, True, False, True]),
        "min_rank": (ap.min_rank("x"), [3, 1, 1, 4]),
        "max_rank": (ap.max_rank("x"), [3, 2, 2, 4]),
        "average_rank": (ap.average_rank("x"), [3.0, 1.5, 1.5, 4.0]),
        "dense_rank": (ap.dense_rank("x"), [2, 1, 1, 3]),
        "cum_sum": (ap.cum_sum("x"), [3, 4, 5, 9]),
        "cum_min": (ap.cum_min("x"), [3, 1, 1, 1]),
        "cum_max": (ap.cum_max("x"), [3, 3, 3, 4]),
        "cum_prod": (ap.cum_prod("x"), [3, 3, 3, 12]),
        "cumany": (ap.cumany("flag"), [False, True, True, True]),
        "cumall": (ap.cumall("flag"), [False, False, False, False]),
        "between": (ap.between("x", 1, 3), [True, True, True, False]),
        "if_else": (ap.if_else(ap.col("x") == 1, "one", "other"), ["other", "one", "one", "other"]),
        "case_when": (
            ap.case_when((ap.col("x") == 1, "one"), (ap.col("x") == 3, "three"), default="other"),
            ["three", "one", "one", "other"],
        ),
        "case_match": (
            ap.case_match("x", ((1, 2), "small"), (3, "three"), default="other"),
            ["three", "small", "small", "other"],
        ),
        "recode": (ap.recode("x", {1: 10}), [3, 10, 10, 4]),
        "near": (ap.near("x", ap.col("x") + 1e-10), [True, True, True, True]),
    }
    for name, (expression, expected) in row_cases.items():
        assert _evaluate(frame, expression).tolist() == expected, name

    assert _evaluate(frame, ap.percent_rank("x")).tolist() == pytest.approx([2 / 3, 0, 0, 1])
    assert _evaluate(frame, ap.cume_dist("x")).tolist() == pytest.approx([0.75, 0.5, 0.5, 1])
    assert _evaluate(frame, ap.cummean("x")).tolist() == pytest.approx([3, 2, 5 / 3, 2.25])
    assert len(_evaluate(frame, ap.ntile("x", 3))) == len(frame)
    with pytest.raises(UnknownColumnError, match="positive integer"):
        ap.ntile("x", 0)


def test_selector_algebra_preserves_public_order_and_hides_virtual_columns() -> None:
    frame = pd.DataFrame({"alpha1": [1], "alpha2": [2], "beta": [3], "text": ["x"]})
    frame["__annplyr_row_number__"] = 1
    frame.attrs["annplyr_virtual_columns"] = {"__annplyr_row_number__"}

    def columns(selector: object) -> list[str]:
        return evaluate_select(frame, selector).columns.tolist()

    assert columns(None) == ["alpha1", "alpha2", "beta", "text"]
    assert columns([]) == []
    assert columns(ap.starts_with("alpha") | ap.ends_with("ta")) == ["alpha1", "alpha2", "beta"]
    assert columns(ap.contains(r"alpha[12]", regex=True) & ap.matches(r"2$")) == ["alpha2"]
    assert columns(~ap.all_of(["alpha2", "text"])) == ["alpha1", "beta"]
    assert columns(~~ap.all_of("alpha1")) == ["alpha1"]
    assert columns(ap.everything() & ap.any_of(["beta", "missing"])) == ["beta"]
    assert columns(ap.last_col()) == ["text"]
    assert columns(ap.last_col(1) | ap.pick("alpha1")) == ["beta", "alpha1"]
    assert columns(ap.num_range("alpha", [1, 2], width=0)) == ["alpha1", "alpha2"]

    with pytest.raises(UnknownColumnError, match="Unknown column"):
        columns(ap.all_of("missing"))
    with pytest.raises(UnknownColumnError, match="outside"):
        columns(ap.last_col(99))
    with pytest.raises(SelectionError, match="Selection failed"):
        columns(42)
    with pytest.raises(SelectionError, match="must return one bool"):
        columns(ap.where(lambda series: series.dtype))


def test_across_and_if_any_all_cover_empty_multi_function_and_duplicate_cases() -> None:
    frame = pd.DataFrame({"a": [1, -1], "b": [2, -2], "text": ["x", "y"]})
    selected = evaluate_assignments(frame, ap.across(ap.starts_with("a")))
    assert selected["a"].tolist() == [1, -1]

    mapped = evaluate_assignments(
        frame,
        {"ignored": ap.across(["a", "b"], {"double": lambda name: ap.col(name) * 2}, names="{col}_{fn}")},
    )
    assert mapped.columns.tolist() == ["a_double", "b_double"]
    assert mapped["b_double"].tolist() == [4, -4]

    sequenced = evaluate_assignments(
        frame,
        ap.across("a", [lambda name: ap.col(name) + 1, lambda name: ap.col(name) - 1]),
    )
    assert sequenced.columns.tolist() == ["a_fn1", "a_fn2"]

    assert evaluate_filter(frame, ap.if_any(["a", "b"], lambda name: ap.col(name) > 0)).tolist() == [0]
    assert evaluate_filter(frame, ap.if_all(["a", "b"], lambda name: ap.col(name) > 0)).tolist() == [0]
    assert evaluate_filter(frame, ap.if_any(ap.any_of("missing"), lambda name: ap.col(name) > 0)).tolist() == []
    assert evaluate_filter(frame, ap.if_all(ap.any_of("missing"), lambda name: ap.col(name) > 0)).tolist() == [0, 1]

    with pytest.raises(DuplicateNameError, match="duplicate output"):
        evaluate_assignments(frame, ap.across(["a", "b"], names="same"))
    with pytest.raises(UnknownColumnError, match="across fns"):
        evaluate_assignments(frame, ap.across("a", fns=[object()]))


@pytest.mark.parametrize(
    ("indexer", "message"),
    [
        ([[0, 1]], "one-dimensional"),
        ([True], "wrong length"),
        ([0.5], "integer positions"),
        ([3], "out of bounds"),
    ],
)
def test_source_positions_reject_invalid_projection_before_read(indexer: object, message: str) -> None:
    with pytest.raises(SelectionError, match=message):
        _positions(indexer, 3, dimension="row")
    assert _positions([-1, 0], 3, dimension="row").tolist() == [2, 0]


def test_dense_sparse_and_backed_adapters_preserve_projection_order_and_schema() -> None:
    vector = adapter_for(np.array([4, 5, 6]), names=["value"])
    assert vector.shape == (3, 1)
    assert vector.read(np.array([2, 0]), np.array([0])).values.tolist() == [[6], [4]]

    empty = adapter_for(sparse.csr_matrix((2, 0)), names=[])
    assert empty.read(np.array([0, 1]), np.array([], dtype=np.intp)).shape == (0, 0)

    class BackedLike:
        shape = (3, 4)
        dtype = np.dtype("int64")

        def __init__(self) -> None:
            self.values = np.arange(12).reshape(self.shape)

        def __getitem__(self, key: object) -> object:
            return self.values[key]

    backed = adapter_for(BackedLike(), names=list("abcd"))
    full_rows = backed.read(np.array([0, 1, 2]), np.array([3, 1]))
    full_columns = backed.read(np.array([2, 0]), np.array([0, 1, 2, 3]))
    duplicate_order = backed.read(np.array([2, 0, 2]), np.array([3, 1, 3]))
    assert full_rows.values.tolist() == [[3, 1], [7, 5], [11, 9]]
    assert full_columns.values.tolist() == [[8, 9, 10, 11], [0, 1, 2, 3]]
    assert duplicate_order.values.tolist() == [[11, 9, 11], [3, 1, 3], [11, 9, 11]]

    with pytest.raises(UnknownSourceError, match="shape"):
        adapter_for(object())
    with pytest.raises(SelectionError, match="width"):
        adapter_for(np.zeros((2, 2)), names=["only_one"])


def test_source_registry_supports_every_aligned_container_and_typed_failures() -> None:
    adata = AnnData(
        X=np.arange(6).reshape(3, 2),
        obs=pd.DataFrame({"group": ["a", "b", "a"]}, index=list("abc")),
        var=pd.DataFrame({"kind": ["x", "y"]}, index=list("uv")),
    )
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["dense"] = np.arange(6).reshape(3, 2)
    adata.obsm["table"] = pd.DataFrame({"one": [1, 2, 3]}, index=adata.obs_names)
    adata.varm["dense"] = np.arange(4).reshape(2, 2)
    adata.obsp["graph"] = sparse.eye(3, format="csr")
    adata.varp["graph"] = sparse.eye(2, format="csc")
    adata.uns["series"] = pd.Series([1, 2], name="value")
    adata.uns["mapping"] = {"left": [1], "right": [2]}

    cases = [
        ("x", None, None, (3, 2)),
        ("x", None, "counts", (3, 2)),
        ("raw", None, None, (3, 2)),
        ("obs", None, None, (3, 1)),
        ("var", None, None, (2, 1)),
        ("obsm", "dense", None, (3, 2)),
        ("obsm", "table", None, (3, 1)),
        ("varm", "dense", None, (2, 2)),
        ("obsp", "graph", None, (3, 3)),
        ("varp", "graph", None, (2, 2)),
        ("uns", "series", None, (2, 1)),
        ("uns", "mapping", None, (1, 2)),
    ]
    for source, key, layer, expected_shape in cases:
        assert source_adapter(adata, source, key=key, layer=layer).shape == expected_shape

    failures: list[tuple[Callable[[], object], str]] = [
        (lambda: source_adapter(adata, "x", layer="missing"), "Unknown layer"),
        (lambda: source_adapter(adata, "raw", layer="counts"), "does not support layer"),
        (lambda: source_adapter(adata, "obsm"), "requires a key"),
        (lambda: source_adapter(adata, "obsm", key="missing"), "Unknown obsm key"),
        (lambda: source_adapter(adata, "uns"), "requires a key"),
        (lambda: source_adapter(adata, "uns", key="missing"), "Unknown uns key"),
        (lambda: source_adapter(adata, "missing"), "Unknown AnnData source"),
    ]
    adata.uns["invalid"] = 42
    failures.append((lambda: source_adapter(adata, "uns", key="invalid"), "cannot be represented"))
    for call, message in failures:
        with pytest.raises(UnknownSourceError, match=message):
            call()

    without_raw = AnnData(X=np.ones((1, 1)))
    with pytest.raises(UnknownSourceError, match="no raw"):
        source_adapter(without_raw, "raw")


def test_request_resolution_and_planner_charge_are_schema_first_and_cumulative() -> None:
    frame = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6]})
    adapter = adapter_for(frame)
    schema = adapter.schema
    assert resolve_columns(adapter, ap.desc(ap.col("b"))).tolist() == [1]
    assert resolve_columns(adapter, {"sum": ap.col("a") + ap.col("c")}).tolist() == [0, 2]
    assert resolve_columns(adapter, nw.col("b") + 1).tolist() == [0, 1, 2]
    assert request_dependencies([ap.col("a"), ap.col("b")], schema) == frozenset({"a", "b"})
    assert request_dependencies([ap.col("a"), nw.col("b")], schema) is None
    with pytest.raises(UnknownColumnError, match="Unknown source column"):
        resolve_columns(adapter, ap.col("missing"))
    assert _name_positions(["a", "a", "b"], ["a", "a"]).tolist() == [0, 1]
    with pytest.raises(UnknownColumnError, match="over-selected"):
        _name_positions(["a"], ["a", "a"])

    planner = RequestPlanner(max_matrix_values=4)
    charged = planner.add(adapter, request=["c", "a"], mode="selection", row_positions=[1, 0], context="charged")
    uncharged = planner.add(adapter, request="b", row_positions=slice(None), context="schema", charge=False)
    assert planner.projected_cells == 4
    result = planner.execute()
    assert result[charged].values.tolist() == [[6, 2], [5, 1]]
    assert result[uncharged]["b"].tolist() == [3, 4]

    rejected = RequestPlanner(max_matrix_values=3)
    rejected.add(adapter, request=["a", "b"], mode="selection", context="matrix export")
    with pytest.raises(AnnplyrError, match="matrix export.*4 matrix values"):
        rejected.execute()


def test_frame_registry_preserves_shapes_and_reports_typed_source_errors() -> None:
    adata = AnnData(
        X=np.arange(6).reshape(3, 2),
        obs=pd.DataFrame({"group": ["a", "b", "a"]}, index=list("abc")),
        var=pd.DataFrame({"kind": ["x", "y"]}, index=list("uv")),
    )
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["embedding"] = np.arange(6).reshape(3, 2)
    adata.varm["loading"] = np.arange(4).reshape(2, 2)
    adata.obsp["graph"] = sparse.eye(3, format="csr")
    adata.varp["graph"] = sparse.eye(2, format="csc")
    adata.uns["table"] = pd.DataFrame({"value": [1, 2]})
    adata.uns["series"] = pd.Series([1, 2])
    adata.uns["mapping"] = {"value": [1, 2]}

    cases = [
        ("obs", None, None, (3, 2)),
        ("var", None, None, (2, 2)),
        ("x", None, "counts", (3, 2)),
        ("raw", None, None, (3, 2)),
        ("obsm", "embedding", None, (3, 2)),
        ("varm", "loading", None, (2, 2)),
        ("obsp", "graph", None, (3, 3)),
        ("varp", "graph", None, (2, 2)),
        ("uns", "table", None, (2, 1)),
        ("uns", "series", None, (2, 1)),
        ("uns", "mapping", None, (2, 1)),
    ]
    for source, key, layer, expected_shape in cases:
        assert source_frame(adata, source, key=key, layer=layer).shape == expected_shape

    assert matrix_frame(pd.DataFrame({"x": [1, 2, 3]}), adata.obs_names).index.tolist() == list("abc")
    assert matrix_frame(np.array([1, 2, 3]), adata.obs_names).shape == (3, 1)
    assert matrix_frame(sparse.csr_matrix((3, 0)), adata.obs_names).shape == (3, 0)

    with pytest.raises(UnknownSourceError, match="does not support layer"):
        source_frame(adata, "raw", layer="counts")
    for source in ["obsm", "varm", "obsp", "varp", "uns"]:
        with pytest.raises(UnknownSourceError, match="requires a key"):
            source_frame(adata, source)
    with pytest.raises(UnknownSourceError, match="Unknown AnnData source"):
        source_frame(adata, "missing")


def test_assignment_size_and_selection_failures_are_typed() -> None:
    frame = pd.DataFrame({"a": [1, 2, 3]})
    with pytest.raises(SizeMismatchError, match="returned 2 rows"):
        _aligned_assignment_series(pd.Series([1, 2]), frame.index, name="bad")
    with pytest.raises(UnknownColumnError, match="Unknown column in assignment"):
        evaluate_assignments(frame, {"bad": ap.col("missing")})
    with pytest.raises(UnknownColumnError, match="Unknown column in filter"):
        evaluate_filter(frame, ap.col("missing") > 0)
    with pytest.raises(UnknownColumnError, match="Unknown column in selector"):
        evaluate_select(frame, nw.col("missing"))
