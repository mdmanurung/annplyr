from __future__ import annotations

import inspect
from dataclasses import dataclass

import anndata as ad
import narwhals as nw
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap
from annplyr._errors import AnnplyrError, SelectionError
from annplyr._expr import AnnplyrExpr, expression_dependencies
from annplyr._frames import evaluate_select
from annplyr._sources import RequestPlanner, adapter_for, resolve_columns, source_adapter


def test_annplyr_expr_tracks_exact_metadata_through_operators_methods_and_namespaces() -> None:
    expression = (ap.col("left") + ap.col("right") * 2).cast(nw.Float64).alias("total")
    assert isinstance(expression, AnnplyrExpr)
    assert expression.dependencies == frozenset({"left", "right"})
    assert expression.output_width == 1
    assert expression.cardinality == "row"
    assert isinstance(expression.to_narwhals(), nw.Expr)

    string = ap.col("label").str.to_uppercase().alias("upper")
    assert string.dependencies == frozenset({"label"})
    assert string.cardinality == "row"
    aggregate = ap.mean(expression)
    assert aggregate.dependencies == frozenset({"left", "right"})
    assert aggregate.cardinality == "scalar"


def test_raw_narwhals_is_accepted_as_an_opaque_expression() -> None:
    raw = nw.col("x") + 1
    assert expression_dependencies(raw) is None
    frame = pd.DataFrame({"x": [1, 2]})
    selected = evaluate_select(frame, raw.alias("out"))
    assert selected["out"].tolist() == [2, 3]


def test_expression_helpers_return_wrappers_and_virtuals_have_exact_dependencies() -> None:
    helpers = [
        ap.lit(1),
        ap.mean("x"),
        ap.n(),
        ap.coalesce(ap.col("x"), 0),
        ap.if_else(ap.col("x") > 0, 1, 0),
        ap.row_number(),
        ap.obs_names,
        ap.var_names,
    ]
    assert all(isinstance(expression, AnnplyrExpr) for expression in helpers)
    assert ap.obs_names.dependencies == frozenset({"__annplyr_obs_names__"})
    assert ap.var_names.dependencies == frozenset({"__annplyr_var_names__"})


def test_every_public_matrix_reading_accessor_exposes_the_cumulative_budget() -> None:
    methods = [
        "filter",
        "arrange",
        "distinct",
        "mutate",
        "transmute",
        "summarize",
        "summarise",
        "pull",
        "to_df",
        "to_tidy",
        "pivot_longer",
        "as_frame",
    ]
    accessor = type(AnnData(X=np.zeros((1, 1))).ap)
    for method in methods:
        assert "max_matrix_values" in inspect.signature(getattr(accessor, method)).parameters


def test_where_is_dtype_only_and_receives_zero_length_typed_series() -> None:
    seen: list[tuple[int, object]] = []

    def numeric(series: pd.Series) -> bool:
        seen.append((len(series), series.dtype))
        return pd.api.types.is_numeric_dtype(series.dtype)

    frame = pd.DataFrame({"number": [1, 2], "text": pd.Series(["a", "b"], dtype="string")})
    selected = evaluate_select(frame, ap.where(numeric))
    assert selected.columns.tolist() == ["number"]
    assert seen == [(0, frame["number"].dtype), (0, frame["text"].dtype)]

    with pytest.raises(SelectionError, match="value-dependent"):
        evaluate_select(frame, ap.where(lambda series: series.iloc[0] > 0))


def test_mutate_remains_sequential_and_preserves_extension_and_scalar_dtypes() -> None:
    adata = AnnData(
        X=np.ones((3, 1)),
        obs=pd.DataFrame(
            {
                "value": pd.array([1, 2, pd.NA], dtype="Int64"),
                "flag": pd.array([True, pd.NA, False], dtype="boolean"),
            },
            index=["a", "b", "c"],
        ),
    )
    result = adata.ap.mutate(
        obs={
            "value": ap.col("value").fill_null(0) + 1,
            "after_overwrite": ap.col("value") * 2,
            "scalar": ap.lit("x"),
            "flag_copy": ap.col("flag"),
        }
    )
    assert result.obs["value"].tolist() == [2, 3, 1]
    assert result.obs["after_overwrite"].tolist() == [4, 6, 2]
    assert result.obs["scalar"].tolist() == ["x", "x", "x"]
    assert str(result.obs["value"].dtype) == "Int64"
    assert str(result.obs["flag_copy"].dtype) == "boolean"


def test_assignment_batching_combines_only_proven_independent_expressions(monkeypatch) -> None:
    import annplyr._frames as frames

    calls = 0
    original = frames.nw.from_native

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(frames.nw, "from_native", counted)
    adata = AnnData(X=np.ones((3, 1)), obs=pd.DataFrame({"a": [1, 2, 3], "b": [3, 2, 1]}))
    result = adata.ap.mutate(
        obs={
            "left": ap.col("a") + 1,
            "right": ap.col("b") + 1,
            "dependent": ap.col("left") * 2,
        }
    )
    assert result.obs["dependent"].tolist() == [4, 6, 8]
    assert calls == 2


def test_noop_metadata_assignment_preserves_existing_series_object() -> None:
    adata = AnnData(X=np.ones((3, 1)), obs=pd.DataFrame({"value": pd.array([1, 2, pd.NA], dtype="Int64")}))
    result = adata.ap.mutate(obs={"value": ap.col("value")})
    assert str(result.obs["value"].dtype) == "Int64"
    assert result.obs["value"].tolist() == [1, 2, pd.NA]


@pytest.mark.parametrize(
    "value",
    [
        np.arange(20).reshape(4, 5),
        sparse.csr_matrix(np.arange(20).reshape(4, 5)),
        sparse.csc_matrix(np.arange(20).reshape(4, 5)),
        sparse.csr_array(np.arange(20).reshape(4, 5)),
        sparse.csc_array(np.arange(20).reshape(4, 5)),
        pd.DataFrame(np.arange(20).reshape(4, 5), columns=list("abcde")),
    ],
)
def test_adapters_project_reordered_integer_rows_and_columns(value) -> None:
    adapter = adapter_for(value, names=None if isinstance(value, pd.DataFrame) else list("abcde"))
    planner = RequestPlanner(max_matrix_values=6)
    planner.add(adapter, request=["e", "b"], mode="selection", row_positions=[3, 1, 3], context="test")
    frame = planner.execute()[0]
    assert frame.columns.tolist() == ["e", "b"]
    assert frame.astype(int).values.tolist() == [[19, 16], [9, 6], [19, 16]]
    if sparse.issparse(value):
        assert all(isinstance(dtype, pd.SparseDtype) for dtype in frame.dtypes)


@pytest.mark.parametrize("storage", ["csr", "csc"])
def test_anndata_view_adapters_and_shape_outputs_preserve_sparse_storage(storage: str) -> None:
    values = np.arange(20, dtype=np.float32).reshape(4, 5)
    matrix = getattr(sparse, f"{storage}_matrix")(values)
    adata = AnnData(X=matrix, var=pd.DataFrame(index=list("abcde")))
    view = adata[[3, 1], [4, 1]]
    adapter = adapter_for(view.X, names=["e", "b"])
    planner = RequestPlanner(max_matrix_values=4)
    planner.add(adapter, request=["b", "e"], mode="selection", context="view")
    assert planner.execute()[0].sparse.to_dense().values.tolist() == [[16, 19], [6, 9]]

    result = adata.ap.filter(x=ap.col("a") >= 0)
    expected_type = sparse.csr_matrix if storage == "csr" else sparse.csc_matrix
    assert isinstance(result.X, expected_type)


@dataclass
class _ReadSentinel:
    reads: int = 0
    schema_reads: int = 0
    last_rows: list[int] | None = None
    last_columns: list[int] | None = None

    @property
    def shape(self) -> tuple[int, int]:
        return 5, 4

    @property
    def names(self) -> tuple[str, ...]:
        return "a", "b", "c", "d"

    @property
    def dtypes(self) -> tuple[object, ...]:
        return (np.dtype("float64"),) * 4

    @property
    def schema(self) -> pd.DataFrame:
        self.schema_reads += 1
        return pd.DataFrame({name: pd.Series([], dtype="float64") for name in self.names})

    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        self.reads += 1
        self.last_rows = row_positions.tolist()
        self.last_columns = column_positions.tolist()
        return pd.DataFrame(
            np.zeros((len(row_positions), len(column_positions))),
            columns=[self.names[position] for position in column_positions],
        )


@pytest.mark.parametrize("budget", [-1, 4])
def test_request_planner_rejects_before_any_read_including_later_cumulative_failure(budget: int) -> None:
    first = _ReadSentinel()
    later = _ReadSentinel()
    planner = RequestPlanner(max_matrix_values=budget)
    planner.add(first, request="a", row_positions=[0, 1, 2], context="first")
    planner.add(later, request="b", row_positions=[0, 1, 2], context="later")
    with pytest.raises(AnnplyrError):
        planner.execute()
    assert first.reads == 0
    assert later.reads == 0


def test_plain_name_selection_resolves_without_materializing_source_schema() -> None:
    adapter = _ReadSentinel()
    assert resolve_columns(adapter, ["d", "a"], mode="selection").tolist() == [3, 0]
    assert adapter.schema_reads == 0


def test_public_planner_rejects_later_cumulative_request_before_first_read(monkeypatch) -> None:
    import annplyr._verbs as verbs

    first = _ReadSentinel()
    later = _ReadSentinel()
    sentinels = iter([first, later])
    monkeypatch.setattr(verbs, "source_adapter", lambda *args, **kwargs: next(sentinels))
    adata = AnnData(X=np.zeros((5, 4)))
    with pytest.raises(AnnplyrError, match="materialize 10 matrix values"):
        adata.ap.to_df(x="a", raw="b", max_matrix_values=9)
    assert first.reads == 0
    assert later.reads == 0


def test_public_projection_reads_only_requested_columns_in_requested_order(monkeypatch) -> None:
    import annplyr._verbs as verbs

    sentinel = _ReadSentinel()
    monkeypatch.setattr(verbs, "source_adapter", lambda *args, **kwargs: sentinel)
    adata = AnnData(X=np.zeros((5, 4)), obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]))
    result = adata.ap.to_tidy(x=["d", "b"], max_matrix_values=10)
    assert sentinel.reads == 1
    assert sentinel.last_rows == [0, 1, 2, 3, 4]
    assert sentinel.last_columns == [3, 1]
    assert result["feature"].drop_duplicates().tolist() == ["d", "b"]


def test_public_opaque_expression_charges_full_source_and_where_reads_only_schema(monkeypatch) -> None:
    import annplyr._verbs as verbs

    opaque = _ReadSentinel()
    monkeypatch.setattr(verbs, "source_adapter", lambda *args, **kwargs: opaque)
    adata = AnnData(X=np.zeros((5, 4)))
    with pytest.raises(AnnplyrError, match="materialize 20 matrix values"):
        adata.ap.filter(x=nw.col("a") > 0, max_matrix_values=5)
    assert opaque.reads == 0

    schema_only = _ReadSentinel()
    monkeypatch.setattr(verbs, "source_adapter", lambda *args, **kwargs: schema_only)
    with pytest.raises(AnnplyrError, match="materialize 20 matrix values"):
        adata.ap.as_frame(
            "x", select=ap.where(lambda series: pd.api.types.is_numeric_dtype(series.dtype)), max_matrix_values=0
        )
    assert schema_only.reads == 0


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
def test_backed_adapters_project_without_materializing_full_source(tmp_path, storage: str) -> None:
    values = np.arange(30, dtype=np.float32).reshape(6, 5)
    matrix = values if storage == "dense" else getattr(sparse, f"{storage}_matrix")(values)
    path = tmp_path / f"{storage}.h5ad"
    AnnData(X=matrix, var=pd.DataFrame(index=list("abcde"))).write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        adapter = source_adapter(backed, "x")
        planner = RequestPlanner(max_matrix_values=4)
        planner.add(adapter, request=["e", "b"], mode="selection", row_positions=[5, 1], context="backed")
        result = planner.execute()[0]
        assert (
            result.sparse.to_dense().values.tolist() == [[29, 26], [9, 6]]
            if storage != "dense"
            else result.values.tolist() == [[29, 26], [9, 6]]
        )
    finally:
        backed.file.close()
