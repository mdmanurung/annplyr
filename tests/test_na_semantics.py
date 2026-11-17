"""NA in reductions and cumulative functions follows dplyr and base R.

A reduction returns NA when its input holds one, unless ``na_rm=True``. A running
``cum_*`` is NA from the first NA onward. ``n_distinct`` counts NA as a value
unless ``na_rm=True``. Every execution path must agree: the chunked accumulator
(row chunks, dense and sparse column chunks, float32 pairwise sums), the plain
Narwhals fallback, and the obs metadata frame.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap

REDUCTIONS = {
    "sum": ap.sum,
    "mean": ap.mean,
    "median": ap.median,
    "sd": ap.sd,
    "min": ap.min,
    "max": ap.max,
    "n_distinct": ap.n_distinct,
}


def _oracle(operation: str, values: np.ndarray, na_rm: bool) -> float:
    """dplyr's answer for one group, computed in float64."""
    values = values.astype(np.float64)
    missing = np.isnan(values)
    if operation == "n_distinct":
        return float(len(np.unique(values[~missing])) + (0 if na_rm or not missing.any() else 1))
    if na_rm:
        values = values[~missing]
    elif missing.any():
        return np.nan
    if operation == "sum":
        return float(values.sum())
    if not len(values):
        return np.nan
    if operation == "sd":
        return float(values.std(ddof=1)) if len(values) > 1 else np.nan
    return float({"mean": np.mean, "median": np.median, "min": np.min, "max": np.max}[operation](values))


# group a holds the NaN in g0; group c is all NaN in g0; g1 never has one.
GROUPS = np.array(["a", "b", "a", "b", "c", "a", "b", "c", "b"])
G0 = np.array([1, 2, 3, 4, np.nan, np.nan, 7, np.nan, 7], dtype=np.float64)
G1 = np.arange(9, dtype=np.float64) * 2


def _adata(dtype: type, storage: str) -> AnnData:
    values = np.column_stack([G0, G1]).astype(dtype)
    matrix = values if storage == "dense" else getattr(sparse, f"{storage}_matrix")(values)
    obs = pd.DataFrame({"group": pd.Categorical(GROUPS)}, index=[f"c{i}" for i in range(9)])
    return AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=["g0", "g1"]))


def _by_group(result: pd.DataFrame, column: str, by: str | None) -> dict[str, float]:
    if by is None:
        return {"all": float(result[column].iloc[0])}
    return dict(zip(result["group"].astype(str), result[column].astype(float), strict=True))


def _expected(operation: str, values: np.ndarray, by: str | None, na_rm: bool) -> dict[str, float]:
    if by is None:
        return {"all": _oracle(operation, values, na_rm)}
    return {key: _oracle(operation, values[GROUPS == key], na_rm) for key in sorted(set(GROUPS))}


@pytest.mark.parametrize("na_rm", [False, True])
@pytest.mark.parametrize("by", [None, "group"])
@pytest.mark.parametrize("chunks", ["column_chunks", "row_chunks"])
@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_reductions_over_a_matrix_follow_na_rm_on_every_chunk_path(
    monkeypatch, dtype: type, storage: str, chunks: str, by: str | None, na_rm: bool
) -> None:
    import annplyr._verbs as verbs

    if chunks == "row_chunks":
        monkeypatch.setattr(verbs, "DEFAULT_REDUCTION_CHUNK_VALUES", 4)
    adata = _adata(dtype, storage)
    for column, values in {"g0": G0, "g1": G1}.items():
        stats = {name: fun(f"x:{column}", na_rm=na_rm) for name, fun in REDUCTIONS.items()}
        result = adata.ap.summarize(stats, by=by)
        for name in REDUCTIONS:
            actual = _by_group(result, name, by)
            expected = _expected(name, values.astype(dtype), by, na_rm)
            assert actual.keys() == expected.keys()
            np.testing.assert_allclose(
                list(actual.values()), list(expected.values()), rtol=1e-5, err_msg=f"{name} {column}"
            )


@pytest.mark.filterwarnings("ignore:Found complex group-by expression")
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("by", [None, "group"])
def test_derived_inputs_and_the_narwhals_fallback_agree_with_the_chunked_paths(dtype: type, by: str | None) -> None:
    adata = _adata(dtype, "dense")
    for name in ["sum", "mean", "sd"]:
        fun = REDUCTIONS[name]
        derived = adata.ap.summarize({"v": fun(ap.col("x:g0") * 1, na_rm=False)}, by=by)  # generic accumulator
        direct = adata.ap.summarize({"v": fun("x:g0", na_rm=False)}, by=by)  # column chunks
        composite = adata.ap.summarize({"v": fun("x:g0", na_rm=False) + 0}, by=by)  # eager Narwhals expression
        assert pd.isna(derived["v"]).tolist() == pd.isna(direct["v"]).tolist()
        assert pd.isna(composite["v"]).tolist() == pd.isna(direct["v"]).tolist()
        np.testing.assert_allclose(derived["v"].astype(float), direct["v"].astype(float), rtol=1e-5)


@pytest.mark.parametrize("na_rm", [False, True])
def test_x_raw_layers_and_grouped_summaries_agree(na_rm: bool) -> None:
    adata = _adata(np.float32, "csr")
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()
    results = [
        adata.ap.summarize(
            {"m": ap.mean(f"{source}:g0", na_rm=na_rm), "s": ap.sum(f"{source}:g0", na_rm=na_rm)}, by="group"
        )
        for source in ("x", "raw", "counts")
    ]
    grouped = adata.ap.group_by("group").summarize(
        {"m": ap.mean("x:g0", na_rm=na_rm), "s": ap.sum("x:g0", na_rm=na_rm)}
    )
    for other in [*results[1:], grouped]:
        pd.testing.assert_frame_equal(other, results[0])
    assert results[0]["m"].isna().tolist() == ([True, False, True] if not na_rm else [False, False, True])


def test_var_axis_reductions_follow_na_rm() -> None:
    loadings = np.array([[1, 5], [np.nan, 6], [3, 7], [4, np.nan]], dtype=np.float32)
    adata = AnnData(X=np.zeros((2, 4), dtype=np.float32), var=pd.DataFrame({"kind": list("xxyy")}, index=list("abcd")))
    adata.varm["loadings"] = loadings
    keep = adata.ap.summarize({"s": ap.sum("loadings:0"), "m": ap.max("loadings:1")}, by="var:kind")
    skip = adata.ap.summarize(
        {"s": ap.sum("loadings:0", na_rm=True), "m": ap.max("loadings:1", na_rm=True)}, by="var:kind"
    )
    assert keep["s"].isna().tolist() == [True, False] and keep["m"].isna().tolist() == [False, True]
    assert skip["s"].tolist() == [1.0, 7.0] and skip["m"].tolist() == [6.0, 7.0]


def test_default_is_na_propagating_and_na_rm_true_skips() -> None:
    adata = _adata(np.float64, "dense")
    result = adata.ap.summarize({"keep": ap.mean("x:g0"), "skip": ap.mean("x:g0", na_rm=True)})
    assert np.isnan(result.loc[0, "keep"])
    assert result.loc[0, "skip"] == pytest.approx(np.nanmean(G0))


def test_float32_pairwise_sum_propagates_nan_over_many_rows() -> None:
    rng = np.random.default_rng(0)
    values = rng.random(5000, dtype=np.float32).reshape(-1, 1)
    values[4321, 0] = np.nan
    adata = AnnData(X=values, obs=pd.DataFrame({"group": pd.Categorical(["a", "b"] * 2500)}))
    adata.var_names = ["g0"]
    keep = adata.ap.summarize({"s": ap.sum("x:g0"), "m": ap.mean("x:g0")}, by="group")
    skip = adata.ap.summarize({"s": ap.sum("x:g0", na_rm=True), "m": ap.mean("x:g0", na_rm=True)}, by="group")
    hit = "b" if 4321 % 2 else "a"
    assert keep.set_index("group").loc[hit].isna().all()
    assert keep.set_index("group").drop(index=hit).notna().all().all()
    assert skip.notna().all().all()


@pytest.mark.parametrize("dtype", ["Int64", "Float64", "float32", "float64"])
@pytest.mark.parametrize("by", [None, "g"])
def test_obs_metadata_reductions_follow_na_rm(dtype: str, by: str | None) -> None:
    values = pd.array([1, None, 3, 4, 5, 5], dtype=dtype)
    obs = pd.DataFrame({"v": values, "g": ["a", "a", "b", "b", "b", "b"]}, index=[f"c{i}" for i in range(6)])
    adata = AnnData(X=np.zeros((6, 1), dtype=np.float32), obs=obs)
    raw = np.array([1, np.nan, 3, 4, 5, 5])
    groups = np.array(["a", "a", "b", "b", "b", "b"])
    for na_rm in (False, True):
        stats = {name: fun("v", na_rm=na_rm) for name, fun in REDUCTIONS.items()}
        result = adata.ap.summarize(stats, by=by)
        for name in REDUCTIONS:
            expected = (
                [_oracle(name, raw, na_rm)]
                if by is None
                else [_oracle(name, raw[groups == key], na_rm) for key in ("a", "b")]
            )
            np.testing.assert_allclose(result[name].astype(float).to_numpy(), expected, err_msg=f"{name} {na_rm}")


def test_n_distinct_counts_na_as_a_value_like_dplyr() -> None:
    obs = pd.DataFrame({"v": ["x", None, "x", None, "y"], "g": ["a", "a", "a", "b", "b"]})
    adata = AnnData(X=np.zeros((5, 1), dtype=np.float32), obs=obs.set_index(pd.Index([f"c{i}" for i in range(5)])))
    result = adata.ap.summarize({"keep": ap.n_distinct("v"), "drop": ap.n_distinct("v", na_rm=True)}, by="g")
    assert result["keep"].tolist() == [2, 2]
    assert result["drop"].tolist() == [1, 1]
    assert adata.ap.summarize({"keep": ap.n_distinct("v")}).loc[0, "keep"] == 3


def test_first_last_nth_return_the_element_even_when_it_is_na() -> None:
    obs = pd.DataFrame({"v": [np.nan, 2.0, np.nan]}, index=["a", "b", "c"])
    adata = AnnData(X=np.zeros((3, 1), dtype=np.float32), obs=obs)
    result = adata.ap.summarize({"f": ap.first("v"), "l": ap.last("v"), "n": ap.nth("v", 1)})
    assert pd.isna(result.loc[0, "f"]) and pd.isna(result.loc[0, "l"]) and result.loc[0, "n"] == 2.0


def test_n_counts_every_row_including_na() -> None:
    adata = _adata(np.float64, "dense")
    assert adata.ap.summarize({"n": ap.n()}).loc[0, "n"] == 9


def test_integer_reductions_keep_their_dtype_when_nothing_is_missing() -> None:
    obs = pd.DataFrame({"v": [1, 2, 3, 4]}, index=list("abcd"))
    adata = AnnData(X=np.zeros((4, 1), dtype=np.float32), obs=obs)
    result = adata.ap.summarize({"s": ap.sum("v"), "lo": ap.min("v"), "hi": ap.max("v"), "n": ap.n_distinct("v")})
    assert [str(result[c].dtype) for c in ["s", "lo", "hi", "n"]] == ["int64"] * 4
    assert result.loc[0, "s"] == 10


def test_reductions_in_mutate_and_filter_propagate_na_per_group() -> None:
    obs = pd.DataFrame({"v": [1.0, np.nan, 3.0, 4.0], "g": ["a", "a", "b", "b"]}, index=list("abcd"))
    adata = AnnData(X=np.zeros((4, 1), dtype=np.float32), obs=obs)
    out = adata.ap.group_by("g").mutate(m=ap.mean("v"), m_rm=ap.mean("v", na_rm=True)).ungroup()
    assert out.obs["m"].isna().tolist() == [True, True, False, False]
    assert out.obs["m_rm"].tolist() == [1.0, 1.0, 3.5, 3.5]
    kept = adata.ap.group_by("g").filter(ap.col("v") > ap.mean("v")).ungroup()
    assert kept.obs_names.tolist() == ["d"]  # group a compares with NA, so it drops out


def test_across_uses_the_dplyr_default_and_a_lambda_sets_na_rm() -> None:
    obs = pd.DataFrame({"u": [1.0, np.nan, 3.0], "w": [1.0, 2.0, 3.0]}, index=list("abc"))
    adata = AnnData(X=np.zeros((3, 1), dtype=np.float32), obs=obs)
    default = adata.ap.summarize({"m": ap.across(["u", "w"], ap.mean)})
    assert np.isnan(default.loc[0, "u"]) and default.loc[0, "w"] == 2.0
    skipped = adata.ap.summarize({"m": ap.across(["u", "w"], lambda c: ap.mean(c, na_rm=True))})
    assert skipped.loc[0, "u"] == 2.0


CUMULATIVE = {
    "cum_sum": (ap.cum_sum, lambda x: np.cumsum(x)),
    "cum_prod": (ap.cum_prod, lambda x: np.cumprod(x)),
    "cum_min": (ap.cum_min, lambda x: np.minimum.accumulate(x)),
    "cum_max": (ap.cum_max, lambda x: np.maximum.accumulate(x)),
    "cummean": (ap.cummean, lambda x: np.cumsum(x) / np.arange(1, len(x) + 1)),
}


@pytest.mark.parametrize("dtype", ["float64", "float32", "Float64", "Int64"])
@pytest.mark.parametrize("name", list(CUMULATIVE))
def test_cumulative_functions_are_na_from_the_first_na_onward(name: str, dtype: str) -> None:
    fun, reference = CUMULATIVE[name]
    values = [2, 1, np.nan, 3, 1, np.nan] if dtype != "Int64" else pd.array([2, 1, None, 3, 1, None], dtype="Int64")
    obs = pd.DataFrame({"v": pd.Series(values, dtype=dtype, index=list("abcdef"))})
    adata = AnnData(X=np.zeros((6, 1), dtype=np.float32), obs=obs)
    out = adata.ap.mutate(r=fun("v")).obs["r"].astype(float).to_numpy()
    clean = reference(np.array([2.0, 1.0]))
    np.testing.assert_allclose(out[:2], clean, rtol=1e-6)
    assert np.isnan(out[2:]).all()


def test_cumulative_functions_restart_per_group_and_keep_integers() -> None:
    obs = pd.DataFrame(
        {"v": [1.0, np.nan, 2.0, 3.0], "i": [1, 2, 3, 4], "g": ["a", "a", "b", "b"]},
        index=list("abcd"),
    )
    adata = AnnData(X=np.zeros((4, 1), dtype=np.float32), obs=obs)
    out = adata.ap.group_by("g").mutate(s=ap.cum_sum("v"), t=ap.cum_sum("i")).ungroup().obs
    np.testing.assert_array_equal(out["s"].to_numpy(), [1.0, np.nan, 2.0, 5.0])
    assert out["t"].tolist() == [1, 3, 3, 7]
    assert str(out["t"].dtype) == "int64"


# aggregate_cells is not a dplyr verb: it sums or averages a matrix over cells
# with matrix products, so a NaN reaches its group and nothing drops it.
@pytest.mark.parametrize("storage", ["dense", "csr"])
@pytest.mark.parametrize("fun", ["sum", "mean"])
def test_aggregate_cells_propagates_nan_into_its_group_only(storage: str, fun: str) -> None:
    adata = _adata(np.float64, storage)
    pseudobulk = ap.aggregate_cells(adata, "group", source="x", fun=fun)
    values = pseudobulk.X.toarray() if sparse.issparse(pseudobulk.X) else pseudobulk.X
    assert np.isnan(values[:, 0]).tolist() == [True, False, True]  # a and c hold NaN in g0, b does not
    assert not np.isnan(values[:, 1]).any()
    assert pseudobulk.obs["n_cells"].tolist() == [3, 4, 2]  # NaN cells still count


def test_obs_derived_input_and_the_flat_function_follow_na_rm() -> None:
    obs = pd.DataFrame(
        {"u": [1.0, np.nan, 3.0, 4.0], "w": [1.0, 1.0, 1.0, 1.0], "g": ["a", "a", "b", "b"]}, index=list("abcd")
    )
    adata = AnnData(X=np.zeros((4, 1), dtype=np.float32), obs=obs)
    derived = adata.ap.summarize({"m": ap.mean(ap.col("u") + ap.col("w"))}, by="g")
    assert derived["m"].isna().tolist() == [True, False]
    assert ap.summarize(adata, {"m": ap.mean("u")})["m"].isna().tolist() == [True]
    assert ap.summarize(adata, {"m": ap.mean("u", na_rm=True)}).loc[0, "m"] == pytest.approx(8 / 3)
