"""annplyr against independent pandas / NumPy / SciPy / Scanpy references on real PBMC3k.

Processed PBMC3k has a scaled dense ``X`` (1,838 highly variable genes) and a
log-normalized sparse ``.raw`` (13,714 genes). Neither is a count matrix, so the
checks are about identical arithmetic and alignment, not biology. Every
expected value comes from pandas, NumPy, SciPy or ``scanpy.get`` rather than
from annplyr. Tolerances: ``RTOL_EXACT`` for selection and copying (values must
not change), ``RTOL_F32`` for float32 element-wise arithmetic, ``RTOL_AGG`` for
aggregates, whose accumulation order may differ from the reference.

Needs network (or ``./data/pbmc3k_processed.h5ad``) on first use.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap

sc = pytest.importorskip("scanpy", reason="Scanpy is installed only in the integration environment")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.filterwarnings("ignore::FutureWarning"),
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
]

RTOL_EXACT = 0.0
RTOL_F32 = 1e-6  # one float32 rounding step (eps ~ 1.2e-7) with margin
RTOL_AGG = 1e-5  # float32 inputs summed in a different order
GENES = ["CST3", "NKG7", "MS4A1"]  # present in both X (HVG) and raw
RAW_ONLY = "LYZ"  # in raw, not among the 1,838 highly variable genes
KEY = "louvain"


# --------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def _pbmc() -> AnnData:
    return sc.datasets.pbmc3k_processed()


@pytest.fixture
def pbmc(_pbmc: AnnData) -> AnnData:
    return _pbmc.copy()


def gene(adata: AnnData, name: str, *, raw: bool) -> pd.Series:
    """Independent reference read of one gene: scanpy's own getter, float64."""
    return sc.get.obs_df(adata, keys=[name], use_raw=raw)[name].astype("float64")


def by_group(adata: AnnData) -> Any:
    return adata.obs.groupby(KEY, sort=False, observed=True)


def assert_same_cells(actual: AnnData, expected_names: pd.Index) -> None:
    assert actual.obs_names.tolist() == expected_names.tolist()


# ------------------------------------------------------------------ 1. values per verb
@pytest.mark.parametrize(
    ("verb", "predicates", "mask"),
    [
        pytest.param(
            "filter",
            [ap.col("n_genes") >= 1000, ap.col(f"raw:{RAW_ONLY}") > 1.0],
            lambda a: (a.obs["n_genes"] >= 1000) & (gene(a, RAW_ONLY, raw=True) > 1.0),
            id="filter-obs-and-raw",
        ),
        pytest.param(
            "filter",
            [ap.col("x:CST3") > 0],
            lambda a: gene(a, "CST3", raw=False) > 0,
            id="filter-x",
        ),
    ],
)
def test_filter_family_selects_the_same_cells_as_a_boolean_mask(
    pbmc: AnnData, verb: str, predicates: list[Any], mask: Callable[[AnnData], pd.Series]
) -> None:
    keep = mask(pbmc).to_numpy()
    result = getattr(pbmc.ap, verb)(*predicates)

    assert 0 < result.n_obs < pbmc.n_obs
    assert_same_cells(result, pbmc.obs_names[keep])
    np.testing.assert_allclose(result.X, pbmc.X[keep], rtol=RTOL_EXACT)
    np.testing.assert_allclose(result.raw.X.toarray(), pbmc.raw.X[keep].toarray(), rtol=RTOL_EXACT)
    assert result.n_vars == pbmc.n_vars  # filter_* subset observations only


def test_grouped_filter_uses_group_local_statistics(pbmc: AnnData) -> None:
    keep = (pbmc.obs["n_counts"] >= pbmc.obs.groupby(KEY, observed=True)["n_counts"].transform("mean")).to_numpy()

    result = pbmc.ap.group_by(KEY).filter(ap.col("n_counts") >= ap.mean("n_counts")).ungroup()

    # api-contract-v0.3.md (`grouped.filter`): "group plan then positional survivors",
    # i.e. first-seen group order, original order inside each group.
    rank = pbmc.obs[KEY].map({k: i for i, k in enumerate(pbmc.obs[KEY].unique())}).to_numpy()[keep]
    assert_same_cells(result, pbmc.obs_names[keep][np.argsort(rank, kind="stable")])


def test_select_keeps_named_columns_and_genes_in_requested_order(pbmc: AnnData) -> None:
    result = pbmc.ap.select("louvain", "n_genes", "x:NKG7", "x:CST3")

    assert result.obs.columns.tolist() == ["louvain", "n_genes"]
    assert result.var_names.tolist() == ["NKG7", "CST3"]
    pd.testing.assert_frame_equal(result.obs, pbmc.obs[["louvain", "n_genes"]])
    np.testing.assert_allclose(result.X, pbmc[:, ["NKG7", "CST3"]].X, rtol=RTOL_EXACT)
    assert result.raw.n_vars == pbmc.raw.n_vars  # select(x=) does not touch .raw


def test_arrange_matches_stable_pandas_sort(pbmc: AnnData) -> None:
    # n_genes has many ties, so this also pins the tie rule (stable, original order).
    keys = pd.DataFrame({"k": pbmc.obs[KEY], "g": pbmc.obs["n_genes"]})
    expected = keys.sort_values(["k", "g"], ascending=[True, False], kind="stable")
    assert_same_cells(pbmc.ap.arrange(KEY, ap.desc("n_genes")), expected.index)

    raw_cst3 = gene(pbmc, "CST3", raw=True)
    by_raw = pbmc.ap.arrange(ap.desc("raw:CST3"))
    assert_same_cells(by_raw, raw_cst3.sort_values(ascending=False, kind="stable").index)


def test_arrange_orders_a_categorical_by_category_order(pbmc: AnnData) -> None:
    expected = pbmc.obs.sort_values(KEY, kind="stable")
    assert_same_cells(pbmc.ap.arrange(KEY), expected.index)


def test_mutate_and_transmute_match_vectorised_arithmetic(pbmc: AnnData) -> None:
    lyz = gene(pbmc, RAW_ONLY, raw=True)
    out = pbmc.ap.mutate({"kcounts": ap.col("n_counts") / 1000, "lyz2": ap.col(f"raw:{RAW_ONLY}") * 2})

    np.testing.assert_allclose(out.obs["kcounts"], pbmc.obs["n_counts"].astype("float64") / 1000, rtol=RTOL_F32)
    np.testing.assert_allclose(out.obs["lyz2"], lyz * 2, rtol=RTOL_F32)
    pd.testing.assert_frame_equal(out.obs[pbmc.obs.columns], pbmc.obs)  # existing columns untouched
    assert "kcounts" not in pbmc.obs  # input unchanged

    only = pbmc.ap.transmute({"lyz2": ap.col(f"raw:{RAW_ONLY}") * 2})
    assert only.obs.columns.tolist() == ["lyz2"]
    np.testing.assert_allclose(only.obs["lyz2"], lyz * 2, rtol=RTOL_F32)

    flag = pbmc.ap.mutate({"high": ap.col("x:CST3") > 0}, into="obs")
    assert flag.obs["high"].tolist() == (gene(pbmc, "CST3", raw=False) > 0).tolist()


@pytest.mark.parametrize("source", ["raw", "x"])
def test_summarize_ungrouped_and_grouped_match_pandas(pbmc: AnnData, source: str) -> None:
    values = gene(pbmc, "CST3", raw=source == "raw")
    stats = {
        "n": ap.n(),
        "mean": ap.mean(f"{source}:CST3"),
        "sd": ap.sd(f"{source}:CST3"),
        "median": ap.median(f"{source}:CST3"),
        "total": ap.sum(f"{source}:CST3"),
        "lo": ap.min(f"{source}:CST3"),
        "hi": ap.max(f"{source}:CST3"),
        "distinct": ap.n_distinct("n_genes"),
    }

    whole = pbmc.ap.summarize(stats)
    expected = {
        "n": len(values), "mean": values.mean(), "sd": values.std(ddof=1), "median": values.median(),
        "total": values.sum(), "lo": values.min(), "hi": values.max(), "distinct": pbmc.obs["n_genes"].nunique(),
    }  # fmt: skip
    for name, value in expected.items():
        np.testing.assert_allclose(whole[name].iloc[0], value, rtol=RTOL_AGG, atol=1e-5, err_msg=name)

    grouped = pbmc.ap.group_by(KEY).summarize(stats)
    frame = pd.DataFrame({"v": values, KEY: pbmc.obs[KEY], "g": pbmc.obs["n_genes"]})
    ref = frame.groupby(KEY, sort=False, observed=True)
    assert grouped[KEY].tolist() == list(ref.size().index)  # first-seen order
    for name, series in {
        "n": ref.size(), "mean": ref["v"].mean(), "sd": ref["v"].std(ddof=1), "median": ref["v"].median(),
        "total": ref["v"].sum(), "lo": ref["v"].min(), "hi": ref["v"].max(), "distinct": ref["g"].nunique(),
    }.items():  # fmt: skip
        np.testing.assert_allclose(grouped[name], series.to_numpy(), rtol=RTOL_AGG, atol=1e-5, err_msg=name)

    # same numbers through the `by=` spelling
    by = pbmc.ap.summarize(stats, by=KEY)
    pd.testing.assert_frame_equal(by, grouped)


def test_count_and_add_count_match_value_counts(pbmc: AnnData) -> None:
    ref = pbmc.obs[KEY].value_counts(sort=False).reindex(pbmc.obs[KEY].unique())

    counted = pbmc.ap.count(KEY)
    assert counted[KEY].tolist() == ref.index.tolist()
    assert counted["n"].tolist() == ref.tolist()
    assert pbmc.ap.count(KEY, sort=True)["n"].tolist() == sorted(ref.tolist(), reverse=True)
    assert pbmc.ap.count()["n"].iloc[0] == pbmc.n_obs

    weighted = pbmc.ap.count(KEY, wt=ap.col("n_counts"))
    np.testing.assert_allclose(weighted["n"], pbmc.obs.groupby(KEY, sort=False, observed=True)["n_counts"].sum())

    added = pbmc.ap.add_count(KEY)
    assert added.obs["n"].tolist() == pbmc.obs[KEY].map(ref).tolist()
    assert "n" not in pbmc.obs


def test_distinct_and_slices_match_pandas(pbmc: AnnData) -> None:
    first = pbmc.obs[~pbmc.obs[KEY].duplicated()]
    assert_same_cells(pbmc.ap.distinct(KEY, keep_all=True), first.index)

    assert_same_cells(pbmc.ap.slice_head(n=7), pbmc.obs_names[:7])
    assert_same_cells(pbmc.ap.slice_tail(n=7), pbmc.obs_names[-7:])
    assert_same_cells(pbmc.ap.slice_min(ap.col("n_counts"), n=5), pbmc.obs["n_counts"].nsmallest(5, keep="first").index)
    assert_same_cells(pbmc.ap.slice_max(ap.col("n_counts"), n=5), pbmc.obs["n_counts"].nlargest(5, keep="first").index)

    # group-local top cell: a global nlargest would pick one group only
    top = pbmc.ap.group_by(KEY).slice_max(ap.col("n_counts"), n=1).ungroup()
    expected = by_group(pbmc)["n_counts"].idxmax()
    assert sorted(top.obs_names) == sorted(expected)


@pytest.mark.parametrize("how", ["left", "inner", "semi", "anti"])
def test_joins_on_obs_match_pandas_merge(pbmc: AnnData, how: str) -> None:
    lineage = pd.DataFrame(
        {
            KEY: ["CD4 T cells", "CD8 T cells", "NK cells", "B cells"],  # myeloid groups deliberately absent
            "lineage": ["T", "T", "NK", "B"],
        }
    )
    in_table = pbmc.obs[KEY].isin(lineage[KEY]).to_numpy()
    mapped = pbmc.obs[KEY].astype(str).map(lineage.set_index(KEY)["lineage"])

    result = getattr(pbmc.ap, f"{how}_join")(lineage, by=KEY)
    if how == "left":
        assert_same_cells(result, pbmc.obs_names)  # no reorder, no duplication
        assert result.obs["lineage"].isna().tolist() == (~in_table).tolist()  # unmatched -> NA
        assert result.obs["lineage"].dropna().tolist() == mapped[in_table].tolist()
    elif how in ("inner", "semi"):
        assert_same_cells(result, pbmc.obs_names[in_table])
        assert ("lineage" in result.obs) == (how == "inner")
    else:
        assert_same_cells(result, pbmc.obs_names[~in_table])
    np.testing.assert_allclose(result.X, pbmc[result.obs_names].X, rtol=RTOL_EXACT)
    # key stays categorical with the original category order; semi/anti go through anndata
    # subsetting, which prunes unused categories, so require an order-preserving subset.
    cats = result.obs[KEY].cat.categories.tolist()
    assert cats == [c for c in pbmc.obs[KEY].cat.categories if c in set(cats)]
    assert set(result.obs[KEY]) <= set(cats)


def test_right_join_with_every_key_equals_left_join(pbmc: AnnData) -> None:
    full = pd.DataFrame({KEY: pbmc.obs[KEY].cat.categories, "idx": range(8)})
    right = pbmc.ap.right_join(full, by=KEY)
    expected = pbmc.obs[KEY].map(dict(zip(full[KEY], full["idx"], strict=True)))
    # right_join orders by the right table (joins.md), then by left position
    order = np.argsort(expected.astype(int).to_numpy(), kind="stable")
    assert_same_cells(right, pbmc.obs_names[order])
    assert right.obs["idx"].tolist() == expected.to_numpy()[order].tolist()


def test_join_on_var_axis_annotates_features_by_name(pbmc: AnnData) -> None:
    panel = pd.DataFrame({"gene": ["MS4A1", "NKG7", RAW_ONLY, "CD79A"], "marker": ["B", "NK", "Mono", "B"]})
    labelled = pbmc.ap.mutate({"gene": ap.var_names}, into="var")
    result = labelled.ap.left_join(panel, by="gene", axis="var")

    expected = pbmc.var_names.to_series().map(panel.set_index("gene")["marker"])
    assert result.var["marker"].tolist() == expected.tolist()
    assert RAW_ONLY not in result.var_names  # LYZ is unmatched in var, not invented
    assert result.var["marker"].notna().sum() == 3


# -------------------------------------------------- 2. expression sources x layout x dtype
@pytest.fixture(scope="module")
def _slice600(_pbmc: AnnData) -> AnnData:
    return _pbmc[:600].copy()  # 600 x 13,714 raw stays small enough to densify in float64


def with_layout(adata: AnnData, *, layout: str, dtype: str) -> AnnData:
    out = adata.copy()
    raw = out.raw.to_adata()
    make = (
        (lambda m: sparse.csr_matrix(m, dtype=dtype))
        if layout == "csr"
        else (lambda m: np.asarray(m.toarray() if sparse.issparse(m) else m, dtype=dtype))
    )
    out.X = make(out.X)
    raw.X = make(raw.X)
    out.raw = raw
    return out


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("layout", ["dense", "csr"])
def test_x_and_raw_sources_are_read_independently(_slice600: AnnData, layout: str, dtype: str) -> None:
    adata = with_layout(_slice600, layout=layout, dtype=dtype)
    # references come from the untouched float32 original, so layout/dtype cannot leak in
    ref_x = gene(_slice600, "CST3", raw=False)
    ref_raw = gene(_slice600, "CST3", raw=True)
    assert not np.allclose(ref_x, ref_raw)  # the two sources genuinely differ

    rtol = RTOL_F32 if dtype == "float32" else RTOL_AGG  # float64 copies of float32 data stay exact
    out = adata.ap.mutate({"x": ap.col("x:CST3"), "raw": ap.col("raw:CST3")})
    np.testing.assert_allclose(out.obs["x"], ref_x, rtol=rtol)
    np.testing.assert_allclose(out.obs["raw"], ref_raw, rtol=rtol)

    summary = adata.ap.group_by(KEY).summarize({"x": ap.mean("x:CST3"), "raw": ap.mean("raw:CST3")})
    frame = pd.DataFrame({KEY: _slice600.obs[KEY], "x": ref_x, "raw": ref_raw})
    expected = frame.groupby(KEY, sort=False, observed=True).mean()
    np.testing.assert_allclose(summary["x"], expected["x"], rtol=RTOL_AGG)
    np.testing.assert_allclose(summary["raw"], expected["raw"], rtol=RTOL_AGG)
    assert not np.allclose(summary["x"], summary["raw"])  # raw: did not silently read X

    kept = adata.ap.filter(ap.col("raw:CST3") > 1.0)
    assert_same_cells(kept, _slice600.obs_names[(ref_raw > 1.0).to_numpy()])
    assert sparse.issparse(kept.X) == (layout == "csr")  # format preserved
    assert kept.X.dtype == dtype


def test_raw_only_gene_is_unreachable_through_x(pbmc: AnnData) -> None:
    assert RAW_ONLY not in pbmc.var_names
    with pytest.raises(ap.UnknownColumnError):
        pbmc.ap.mutate({"v": ap.col(f"x:{RAW_ONLY}")})
    got = pbmc.ap.mutate({"v": ap.col(f"raw:{RAW_ONLY}")}).obs["v"]
    np.testing.assert_allclose(got, gene(pbmc, RAW_ONLY, raw=True), rtol=RTOL_F32)


# ---------------------------------------------------------------------- 3. alignment
def assert_aligned(result: AnnData, original: AnnData) -> None:
    """Every container of ``result`` is the original's, reindexed by name."""
    rows = original.obs_names.get_indexer(result.obs_names)
    cols = original.var_names.get_indexer(result.var_names)
    assert (rows >= 0).all()
    assert (cols >= 0).all()
    assert result.obs_names.is_unique
    np.testing.assert_allclose(result.X, original.X[rows][:, cols], rtol=RTOL_EXACT)
    assert result.raw is not None
    assert result.raw.obs_names.tolist() == result.obs_names.tolist()
    raw_cols = original.raw.var_names.get_indexer(result.raw.var_names)
    np.testing.assert_allclose(result.raw.X.toarray(), original.raw.X[rows][:, raw_cols].toarray(), rtol=RTOL_EXACT)
    for key in original.obsm:
        np.testing.assert_allclose(result.obsm[key], original.obsm[key][rows], rtol=RTOL_EXACT, err_msg=key)
    for key in original.varm:
        np.testing.assert_allclose(result.varm[key], original.varm[key][cols], rtol=RTOL_EXACT, err_msg=key)
    for key in original.obsp:
        want = original.obsp[key][rows][:, rows]
        np.testing.assert_allclose(result.obsp[key].toarray(), want.toarray(), rtol=RTOL_EXACT, err_msg=key)
    for key in original.layers:
        if key is not None:  # anndata 0.13 lists X itself under a None key
            np.testing.assert_allclose(result.layers[key], original.layers[key][rows][:, cols], rtol=RTOL_EXACT)
    pd.testing.assert_frame_equal(
        result.obs[original.obs.columns.intersection(result.obs.columns)],
        original.obs.iloc[rows][original.obs.columns.intersection(result.obs.columns)],
    )
    assert set(result.uns) == set(original.uns)
    np.testing.assert_array_equal(result.uns[f"{KEY}_colors"], original.uns[f"{KEY}_colors"])


ALIGNMENT_PIPELINES: dict[str, Callable[[AnnData], AnnData]] = {
    "filter": lambda a: a.ap.filter(ap.col("n_genes") > 800, ap.col("raw:CST3") > 0),
    "arrange": lambda a: a.ap.arrange(ap.desc("n_counts"), "n_genes"),
    "slice": lambda a: a.ap.slice_max(ap.col("n_counts"), n=50),
    "select-genes": lambda a: a.ap.select("x:NKG7", "x:CST3", "x:MS4A1"),
    "filter-then-select": lambda a: a.ap.filter(ap.col("x:NKG7") > 0).select(
        *[f"x:{g}" for g in a.var_names if g.startswith("C")]
    ),
    "mutate": lambda a: a.ap.mutate({"k": ap.col("n_counts") / 1000}),
    "left_join": lambda a: a.ap.left_join(pd.DataFrame({KEY: ["B cells"], "tag": ["b"]}), by=KEY),
    "grouped-filter": lambda a: a.ap.group_by(KEY).filter(ap.row_number() <= 3).ungroup(),
    "distinct": lambda a: a.ap.distinct(KEY, keep_all=True),
}


@pytest.mark.parametrize("name", list(ALIGNMENT_PIPELINES))
def test_containers_stay_aligned_after_verbs(pbmc: AnnData, name: str) -> None:
    result = ALIGNMENT_PIPELINES[name](pbmc)
    assert result.obsm["X_pca"].shape[0] == result.n_obs
    assert result.varm["PCs"].shape[0] == result.n_vars
    assert result.obsp["connectivities"].shape == (result.n_obs, result.n_obs)
    assert_aligned(result, pbmc)


@pytest.mark.parametrize("name", ["filter", "slice", "mutate", "select-genes", "left_join", "distinct"])
def test_shuffled_input_gives_the_same_cells_and_values_once_aligned_by_name(pbmc: AnnData, name: str) -> None:
    perm = np.random.default_rng(0).permutation(pbmc.n_obs)
    shuffled = pbmc[perm].copy()
    assert shuffled.obs_names.tolist() != pbmc.obs_names.tolist()

    plain = ALIGNMENT_PIPELINES[name](pbmc)
    mixed = ALIGNMENT_PIPELINES[name](shuffled)

    assert_aligned(mixed, pbmc)  # still the right values under the right names
    if name in ("filter", "mutate", "select-genes", "left_join"):  # order-free: same cell set
        assert sorted(mixed.obs_names) == sorted(plain.obs_names)
        common = plain.obs.columns
        pd.testing.assert_frame_equal(mixed.obs.loc[plain.obs_names, common], plain.obs[common])
    if name == "slice":  # top-n by n_counts: ties could swap members, so compare the values selected
        np.testing.assert_array_equal(
            np.sort(mixed.obs["n_counts"].to_numpy()), np.sort(plain.obs["n_counts"].to_numpy())
        )
    if name == "distinct":  # first-seen representative depends on order by design; one cell per group
        assert sorted(mixed.obs[KEY]) == sorted(plain.obs[KEY])


def test_shuffled_summaries_agree_after_sorting_by_group(pbmc: AnnData) -> None:
    perm = np.random.default_rng(1).permutation(pbmc.n_obs)
    shuffled = pbmc[perm].copy()
    stats = {"n": ap.n(), "m": ap.mean("raw:CST3"), "s": ap.sd("x:NKG7"), "q": ap.median("n_counts")}

    plain = pbmc.ap.group_by(KEY).summarize(stats).set_index(KEY).sort_index()
    mixed = shuffled.ap.group_by(KEY).summarize(stats).set_index(KEY).sort_index()

    np.testing.assert_allclose(mixed.to_numpy(float), plain.to_numpy(float), rtol=RTOL_AGG)
    # and the first-seen group order follows the shuffled data, as documented in user_guide/grouping.md
    assert shuffled.ap.count(KEY)[KEY].tolist() == list(shuffled.obs[KEY].unique())


# ------------------------------------------------------------- 4. silent-failure paths
def test_many_to_many_and_cardinality_violations_raise(pbmc: AnnData) -> None:
    twice = pd.DataFrame({KEY: ["B cells", "B cells"], "tag": ["a", "b"]})
    for relationship in ("many-to-one", "many-to-many", "one-to-one"):
        with pytest.raises(ap.JoinRelationshipError):
            pbmc.ap.left_join(twice, by=KEY, relationship=relationship)  # type: ignore[arg-type]
    with pytest.raises(ap.JoinRelationshipError):  # duplicate left keys under one-to-one
        pbmc.ap.left_join(pd.DataFrame({KEY: ["B cells"], "tag": ["a"]}), by=KEY, relationship="one-to-one")
    # joins.md: a right join may not add a cell the object has no data for
    with pytest.raises(ap.JoinRelationshipError):
        pbmc.ap.right_join(pd.DataFrame({KEY: ["Basophils"]}), by=KEY)
    with pytest.raises(ap.JoinRelationshipError):
        pbmc.ap.inner_join(pd.DataFrame({KEY: ["B cells"], "tag": ["a"]}), by=KEY, unmatched="error")
    assert pbmc.n_obs == 2638  # input untouched by the failed joins


def test_unmatched_join_key_is_na_never_a_default(pbmc: AnnData) -> None:
    table = pd.DataFrame({KEY: ["B cells"], "n_tag": [1.0]})
    obs = pbmc.ap.left_join(table, by=KEY).obs
    assert obs["n_tag"].notna().sum() == (pbmc.obs[KEY] == "B cells").sum()
    assert obs.loc[pbmc.obs[KEY] != "B cells", "n_tag"].isna().all()  # joins.md: unmatched -> NA


def test_unused_categories_and_emptied_groups_follow_the_documented_grouping_contract(pbmc: AnnData) -> None:
    # docs/user_guide/grouping.md: groups are first-seen, include NA, keep categorical
    # dtype and order, and omit unobserved categories.
    pbmc.obs[KEY] = pbmc.obs[KEY].cat.add_categories("Unused")
    for table in (pbmc.ap.count(KEY), pbmc.ap.group_by(KEY).summarize({"n": ap.n()})):
        assert "Unused" not in table[KEY].tolist()
        assert table["n"].tolist() == pbmc.obs[KEY].value_counts(sort=False).reindex(pbmc.obs[KEY].unique()).tolist()
        assert isinstance(table[KEY].dtype, pd.CategoricalDtype)

    emptied = pbmc.ap.filter(ap.col(KEY) != "NK cells")  # anndata prunes unused categories on subset
    assert "NK cells" not in emptied.obs[KEY].cat.categories
    assert "NK cells" not in emptied.ap.count(KEY)[KEY].tolist()
    assert "NK cells" not in emptied.ap.group_by(KEY).summarize({"m": ap.mean("raw:CST3")})[KEY].tolist()

    none = pbmc.ap.filter(ap.col("n_genes") > 10**9)
    assert none.n_obs == 0
    assert len(none.ap.group_by(KEY).summarize({"n": ap.n()})) == 0
    whole = none.ap.summarize({"n": ap.n(), "m": ap.mean("n_counts")})
    assert whole["n"].iloc[0] == 0
    assert np.isnan(whole["m"].iloc[0])  # mean of nothing is NA, not 0


def test_all_na_group_gives_na_not_a_number(pbmc: AnnData) -> None:
    pbmc.obs["score"] = pbmc.obs["n_counts"].astype("float64")
    pbmc.obs.loc[pbmc.obs[KEY] == "NK cells", "score"] = np.nan
    out = pbmc.ap.group_by(KEY).summarize({"m": ap.mean("score"), "n": ap.n()}).set_index(KEY)
    assert np.isnan(out.loc["NK cells", "m"])
    assert out.loc["NK cells", "n"] == (pbmc.obs[KEY] == "NK cells").sum()  # n() still counts rows
    ref = pbmc.obs.groupby(KEY, observed=True)["score"].mean()
    np.testing.assert_allclose(out["m"].reindex(ref.index), ref, rtol=RTOL_AGG)


def test_unknown_names_raise_named_errors(pbmc: AnnData) -> None:
    with pytest.raises(ap.UnknownColumnError):
        pbmc.ap.filter(ap.col("no_such_column") > 1)
    with pytest.raises(ap.UnknownColumnError):
        pbmc.ap.filter(ap.col("raw:NOT_A_GENE") > 1)
    with pytest.raises(ap.UnknownColumnError):
        pbmc.ap.summarize({"m": ap.mean("x:NOT_A_GENE")})
    with pytest.raises(ap.UnknownColumnError):
        pbmc.ap.arrange("no_such_column")
    with pytest.raises(ap.UnknownColumnError):
        pbmc.ap.group_by("no_such_column")
    with pytest.raises(ap.UnknownSourceError):
        pbmc.ap.filter(ap.col("not_a_container:CST3") > 1)
    with pytest.raises(ap.UnknownSourceError):
        pbmc.ap.filter(ap.col("X_nope:0") > 1)


@pytest.mark.parametrize(
    "call",
    [
        pytest.param(lambda a: a.ap.slice_max(ap.col("raw:CST3"), n=3), id="slice_max"),
        pytest.param(lambda a: a.ap.slice_min(ap.col("x:CST3"), n=3), id="slice_min"),
    ],
)
def test_slice_by_matrix_reference_fails_with_an_annplyr_error(pbmc: AnnData, call: Callable[[AnnData], Any]) -> None:
    # api-contract-v0.3.md: slice_min/slice_max take a metadata expression, no matrix reads,
    # errors are "selection, column". They must not leak a narwhals exception.
    with pytest.raises(ap.AnnplyrError):
        call(pbmc)


def test_duplicate_names_raise_duplicate_name_error(pbmc: AnnData) -> None:
    with pytest.raises(ap.DuplicateNameError):
        pbmc.ap.select("n_genes", "n_genes")
    with pytest.raises(ap.DuplicateNameError):
        pbmc.ap.rename(n_counts="n_genes")
    with pytest.raises(ap.DuplicateNameError):
        pbmc.ap.mutate({"n_genes": ap.lit(0), "again": ap.across(["n_genes"], names="n_genes")})


def test_wrong_length_assignment_is_rejected(pbmc: AnnData) -> None:
    with pytest.raises(ap.AnnplyrError):
        pbmc.ap.mutate({"short": pd.Series([1, 2, 3])})
    assert "short" not in pbmc.obs


@pytest.mark.parametrize("verb", ["mutate", "add_count", "transmute"])
def test_working_on_a_view_never_touches_the_parent(_pbmc: AnnData, verb: str) -> None:
    parent = _pbmc.copy()
    before_obs = parent.obs.copy()
    before_x = parent.X.copy()
    view = parent[:200]
    assert view.is_view

    if verb == "mutate":
        view.ap.mutate({"z": ap.lit(1)})
    elif verb == "add_count":
        view.ap.add_count(KEY)
    else:
        view.ap.transmute({"z": ap.lit(1)})

    pd.testing.assert_frame_equal(parent.obs, before_obs)
    np.testing.assert_array_equal(parent.X, before_x)
    assert parent.n_obs == _pbmc.n_obs
    assert "z" not in parent.obs
    assert "n" not in parent.obs


def test_subset_results_do_not_alias_the_input(pbmc: AnnData) -> None:
    out = pbmc.ap.filter(ap.col("n_genes") > 0)
    out.X[0, 0] = 12345.0
    out.raw.X.data[0] = 12345.0
    out.obs.iloc[0, out.obs.columns.get_loc("n_genes")] = -1
    assert pbmc.X[0, 0] != 12345.0
    assert pbmc.raw.X.data[0] != 12345.0
    assert pbmc.obs["n_genes"].iloc[0] != -1


# --------------------------------------------------------------------- 5. round trips
@pytest.fixture(scope="module")
def stores(_pbmc: AnnData, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("pbmc3k")
    _pbmc.write_h5ad(root / "pbmc3k.h5ad")
    _pbmc.write_zarr(root / "pbmc3k.zarr")
    return {"h5ad": root / "pbmc3k.h5ad", "zarr": root / "pbmc3k.zarr"}


def verb_bundle(adata: AnnData) -> tuple[AnnData, pd.DataFrame, pd.DataFrame]:
    """Filter, arrange, mutate, join, and summarize in one call chain."""
    table = pd.DataFrame({KEY: ["B cells", "NK cells"], "tag": ["lymph", "lymph"]})
    result = (
        adata.ap.filter(ap.col("n_genes") > 700, ap.col("raw:CST3") > 0)
        .arrange(ap.desc("n_counts"))
        .mutate({"k": ap.col("n_counts") / 1000, "cst3": ap.col("raw:CST3")})
        .left_join(table, by=KEY)
    )
    summary = adata.ap.group_by(KEY).summarize({"n": ap.n(), "m": ap.mean("raw:CST3"), "s": ap.sd("x:NKG7")})
    counts = adata.ap.count(KEY, sort=True)
    return result, summary, counts


def assert_results_identical(a: tuple[Any, ...], b: tuple[Any, ...]) -> None:
    assert a[0].obs_names.tolist() == b[0].obs_names.tolist()
    assert a[0].var_names.tolist() == b[0].var_names.tolist()
    pd.testing.assert_frame_equal(a[0].obs, b[0].obs, check_exact=True)
    np.testing.assert_array_equal(a[0].X, b[0].X)
    np.testing.assert_array_equal(a[0].raw.X.toarray(), b[0].raw.X.toarray())
    np.testing.assert_array_equal(a[0].obsm["X_pca"], b[0].obsm["X_pca"])
    pd.testing.assert_frame_equal(a[1], b[1], check_exact=True)
    pd.testing.assert_frame_equal(a[2], b[2], check_exact=True)


@pytest.mark.parametrize("fmt", ["h5ad", "zarr"])
def test_verbs_give_identical_results_after_a_disk_round_trip(
    _pbmc: AnnData, stores: dict[str, Path], fmt: str
) -> None:
    reader = ad.read_h5ad if fmt == "h5ad" else ad.read_zarr
    reloaded = reader(stores[fmt])

    assert reloaded.obs_names.tolist() == _pbmc.obs_names.tolist()
    pd.testing.assert_frame_equal(reloaded.obs, _pbmc.obs, check_exact=True)
    np.testing.assert_array_equal(reloaded.X, _pbmc.X)
    assert_results_identical(verb_bundle(reloaded), verb_bundle(_pbmc))


def test_verb_output_survives_write_and_read(pbmc: AnnData, tmp_path: Path) -> None:
    result, _, _ = verb_bundle(pbmc)
    for fmt, writer, reader in (("h5ad", result.write_h5ad, ad.read_h5ad), ("zarr", result.write_zarr, ad.read_zarr)):
        path = tmp_path / f"result.{fmt}"
        writer(path)
        back = reader(path)
        assert_aligned(back, pbmc)
        pd.testing.assert_frame_equal(back.obs, result.obs, check_exact=True)


# ------------------------------------------------------------------- 6. determinism
def test_verbs_are_deterministic_across_repeated_runs(pbmc: AnnData) -> None:
    assert_results_identical(verb_bundle(pbmc), verb_bundle(pbmc))
    snapshot = pbmc.obs.copy()
    verb_bundle(pbmc)
    pd.testing.assert_frame_equal(pbmc.obs, snapshot)  # running verbs does not mutate the input


def test_slice_sample_is_reproducible_only_with_a_seed(pbmc: AnnData) -> None:
    first = pbmc.ap.slice_sample(n=100, random_state=7)
    again = pbmc.ap.slice_sample(n=100, random_state=7)
    other = pbmc.ap.slice_sample(n=100, random_state=8)

    assert first.obs_names.tolist() == again.obs_names.tolist()
    assert first.obs_names.tolist() != other.obs_names.tolist()
    grouped_a = pbmc.ap.group_by(KEY).slice_sample(n=5, random_state=3).ungroup()
    grouped_b = pbmc.ap.group_by(KEY).slice_sample(n=5, random_state=3).ungroup()
    assert grouped_a.obs_names.tolist() == grouped_b.obs_names.tolist()
