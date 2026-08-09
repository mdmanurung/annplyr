from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap
from annplyr._errors import AnnplyrError
from annplyr._groups import GroupPlan
from annplyr._sources import RequestPlanner, adapter_for


def _reduction_fixture(matrix: object) -> AnnData:
    obs = pd.DataFrame(
        {
            "group": pd.Categorical(["b", "a", None, "b", "a", None, "b", "a", "b"]),
            "nullable": pd.array([1, pd.NA, 2, 1, pd.NA, 2, 1, pd.NA, 1], dtype="Int64"),
        },
        index=["cell", "cell", "other", "cell", "last", "cell", "other", "cell", "cell"],
    )
    return AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=["g0", "g1", "g2"]))


def _assignments() -> dict[str, object]:
    return {
        "rows": ap.n(),
        "g2_sum": ap.sum("g2"),
        "g0_mean": ap.mean("g0"),
        "g0_sd": ap.sd("g0"),
        "g0_min": ap.min("g0"),
        "g0_max": ap.max("g0"),
        "g0_first": ap.first("g0"),
        "g0_last": ap.last("g0"),
        "g0_median": ap.median("g0"),
        "g0_unique": ap.n_distinct("g0"),
        "derived_mean": ap.mean(ap.col("g0") + ap.col("g1")),
        "all_nonnegative": (ap.col("g2") >= 0).all(),
        "any_positive": (ap.col("g2") > 0).any(),
    }


@pytest.fixture
def reduction_values() -> np.ndarray:
    values = np.array(
        [
            [1, 2, 0],
            [np.nan, 4, 2],
            [3, 6, 0],
            [4, 8, 3],
            [5, 10, 0],
            [6, 12, 4],
            [7, 14, 0],
            [8, 16, 5],
            [9, 18, 0],
        ],
        dtype=np.float32,
    )
    return values


@pytest.mark.parametrize("by", [None, "group", "nullable"])
@pytest.mark.parametrize("storage", ["dense", "csr", "csc", "backed_dense", "backed_csr", "backed_csc"])
@pytest.mark.filterwarnings("ignore:Observation names are not unique")
def test_chunked_reductions_exactly_match_one_chunk_for_every_storage(
    tmp_path,
    monkeypatch,
    reduction_values: np.ndarray,
    storage: str,
    by: str | None,
) -> None:
    import annplyr._verbs as verbs

    monkeypatch.setattr(verbs, "DEFAULT_REDUCTION_CHUNK_VALUES", reduction_values.size + 1)
    expected = _reduction_fixture(reduction_values).ap.summarize(x=_assignments(), by=by)

    memory_storage = storage.removeprefix("backed_")
    matrix = (
        reduction_values if memory_storage == "dense" else getattr(sparse, f"{memory_storage}_matrix")(reduction_values)
    )
    candidate = _reduction_fixture(matrix)
    backed = None
    if storage.startswith("backed_"):
        path = tmp_path / f"{storage}.h5ad"
        candidate.write_h5ad(path)
        backed = ad.read_h5ad(path, backed="r")
        candidate = backed

    try:
        monkeypatch.setattr(verbs, "DEFAULT_REDUCTION_CHUNK_VALUES", 4)
        actual = candidate.ap.summarize(x=_assignments(), by=by)
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
    finally:
        if backed is not None:
            backed.file.close()


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
def test_var_axis_chunking_preserves_nullable_groups_and_assignment_order(
    monkeypatch,
    storage: str,
) -> None:
    import annplyr._verbs as verbs

    values = np.arange(18, dtype=np.float32).reshape(6, 3)
    loadings = values if storage == "dense" else getattr(sparse, f"{storage}_matrix")(values)
    adata = AnnData(
        X=np.zeros((2, 6), dtype=np.float32),
        var=pd.DataFrame(
            {
                "kind": pd.Categorical(["rna", "protein", None, "rna", "protein", None]),
            },
            index=[f"g{position}" for position in range(6)],
        ),
    )
    adata.varm["loadings"] = loadings
    monkeypatch.setattr(verbs, "DEFAULT_REDUCTION_CHUNK_VALUES", 2)
    result = adata.ap.summarize(
        varm={"loadings": {"second": ap.mean("1"), "first": ap.sum("0")}},
        by="kind",
    )
    expected = pd.DataFrame(
        {
            "kind": pd.Categorical(["rna", "protein", None], categories=["protein", "rna"]),
            "second": np.array([5.5, 8.5, 11.5], dtype=np.float32),
            "first": np.array([9.0, 15.0, 21.0], dtype=np.float32),
        }
    )
    pd.testing.assert_frame_equal(result, expected, check_exact=True)


def test_chunk_plan_is_explicit_deterministic_and_preserves_requested_positions() -> None:
    values = np.arange(21).reshape(7, 3)
    planner = RequestPlanner(max_matrix_values=10)
    token = planner.add(
        adapter_for(values, names=["a", "b", "c"]),
        request=["c", "a"],
        mode="selection",
        row_positions=[6, 1, 6, 2, 0],
        context="chunk plan",
    )
    chunks = planner.chunk_plan(token, target_values=5)
    assert chunks.rows_per_chunk == 2
    assert chunks.chunk_count == 3
    assert chunks.columns_per_chunk == 1
    assert chunks.column_chunk_count == 2
    assert [rows.tolist() for rows in chunks.row_chunks()] == [[6, 1], [6, 2], [0]]
    assert [frame.values.tolist() for _, frame in chunks.read_chunks()] == [
        [[20, 18], [5, 3]],
        [[20, 18], [8, 6]],
        [[2, 0]],
    ]
    assert [columns.tolist() for columns, _ in chunks.read_column_chunks()] == [[2], [0]]


@dataclass
class _CountingAdapter:
    values: np.ndarray
    reads: int = 0

    @property
    def shape(self) -> tuple[int, int]:
        return self.values.shape

    @property
    def names(self) -> tuple[str, ...]:
        return "g0", "g1"

    @property
    def dtypes(self) -> tuple[object, ...]:
        return (self.values.dtype, self.values.dtype)

    @property
    def schema(self) -> pd.DataFrame:
        return pd.DataFrame(
            {name: pd.Series([], dtype=dtype) for name, dtype in zip(self.names, self.dtypes, strict=True)}
        )

    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        self.reads += 1
        return pd.DataFrame(
            self.values[np.ix_(row_positions, column_positions)],
            columns=[self.names[position] for position in column_positions],
        )


def test_shared_reductions_read_each_chunk_once_and_budget_before_first_read(monkeypatch) -> None:
    import annplyr._verbs as verbs

    adapter = _CountingAdapter(np.arange(14, dtype=np.float32).reshape(7, 2))
    monkeypatch.setattr(verbs, "source_adapter", lambda *args, **kwargs: adapter)
    monkeypatch.setattr(verbs, "DEFAULT_REDUCTION_CHUNK_VALUES", 4)
    adata = AnnData(
        X=np.zeros((7, 2), dtype=np.float32),
        obs=pd.DataFrame({"group": pd.Categorical(["a", "a", "b", "b", "a", "b", "a"])}),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    result = adata.ap.group_by(obs="group").summarize(
        x={"sum": ap.sum("g0"), "mean": ap.mean("g0"), "other": ap.max("g1")}
    )
    assert result.columns.tolist() == ["group", "sum", "mean", "other"]
    assert adapter.reads == 4

    adapter.reads = 0
    with pytest.raises(AnnplyrError, match="exceeds max_matrix_values=13"):
        adata.ap.summarize(x={"sum": ap.sum("g0"), "other": ap.sum("g1")}, max_matrix_values=13)
    assert adapter.reads == 0


def test_grouped_summary_reuses_the_callers_positional_plan(monkeypatch, reduction_values: np.ndarray) -> None:
    calls = 0
    original = GroupPlan.build.__func__

    def counted(cls, *args, **kwargs):
        nonlocal calls
        calls += 1
        return original(cls, *args, **kwargs)

    monkeypatch.setattr(GroupPlan, "build", classmethod(counted))
    result = (
        _reduction_fixture(reduction_values)
        .ap.group_by(obs="group")
        .summarize(x={"mean": ap.mean("g0"), "sum": ap.sum("g1")})
    )
    assert result.shape == (3, 3)
    assert calls == 1


@pytest.mark.parametrize("dtype", [np.float32, np.int32])
def test_empty_chunked_reductions_have_stable_values_and_dtypes(dtype: object) -> None:
    adata = AnnData(
        X=np.empty((0, 2), dtype=dtype),
        obs=pd.DataFrame({"group": pd.Categorical([])}),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    result = adata.ap.summarize(x={"n": ap.n(), "sum": ap.sum("g0"), "mean": ap.mean("g0")})
    assert result.loc[0, "n"] == 0
    assert result.loc[0, "sum"] == 0
    assert pd.isna(result.loc[0, "mean"])
    assert str(result["n"].dtype) == "int64"
    grouped = adata.ap.summarize(x={"n": ap.n(), "sum": ap.sum("g0")}, by="group")
    assert grouped.empty
    assert str(grouped["group"].dtype) == "category"
