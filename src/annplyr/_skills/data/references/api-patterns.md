# API Patterns

## Accessor Verbs

Wrap with `adata = ap.tidy(adata)` and call `adata.<verb>(...)`. `adata.ap.<verb>(...)` and `ap.<verb>(adata, ...)` are the same call.

For a checked downstream boundary, use one no-op cast; do not import private
accessor or grouped modules:

```python
from typing import cast

from annplyr.typing import AnnDataWithAnnplyr

typed = cast(AnnDataWithAnnplyr, adata)
filtered = typed.filter(obs.qc_pass)
```

Names are positional and carry their container (`obs.x`, `var.x`, `X.x`, `raw.x`, `counts.x`, `X_umap[0]`). Every verb returns a new object; there is no `copy=` or `inplace=`.

**Subsetting and column ops:**
- `filter(*predicates, max_matrix_values=...)` — subset observations or features, one predicate per argument; each names its own container. A `var` predicate subsets features.
- `select(*names)` — keep named metadata columns and selected features: `select(obs.louvain, var.n_cells, X.MS4A1)`. A minus sign drops an `obs` column: `select(-obs.louvain)`. Selectors such as `ap.starts_with("n_")` pick `obs` columns; index a container to pick features: `select(X[ap.starts_with("CD")])`.
- `rename(new=old)` — rename columns by mapping new names to existing columns: `rename(cell_type=obs.louvain)`. A dict handles new names that are not identifiers.
- `rename_with(str.lower, *names)` — rename selected columns by applying a function to names.
- `relocate(*names, before=..., after=...)` — reorder metadata columns. Raises `SelectionError` when the `before`/`after` anchor is among the columns being moved.
- `distinct(*names, keep_all=..., max_matrix_values=...)` — deduplicate by selected columns.
- `arrange(*keys, max_matrix_values=...)` — reorder observations or features, one sort key per argument; wrap a key in `ap.desc` to reverse it.
- `slice(idx)` / `slice_head(n=...)` / `slice_tail(n=...)` / `slice_min(order_by, n=...)` / `slice_max(order_by, n=...)` / `slice_sample(n=..., random_state=...)` — position or rank-based slicing; pass `by=obs.group` for a per-group slice in one call. A `var.` sort key ranks features.

**Mutation and summaries:**
- `mutate(name=value, ..., into="obs")` — write one metadata column per keyword, in order; a later value may read an earlier result. `into="var"` writes the feature axis. Matrix values are read-only. A dict is accepted for names that are not identifiers.
- `transmute(name=value, ..., into="obs")` — like `mutate` but drops all columns not produced by the call.
- `group_by(*by)` — return a persistent `GroupedAnnData` wrapper; call `ungroup()` for AnnData.
- `summarize(name=agg, ..., by=...)` — aggregate by grouping key, one aggregate per keyword; also callable as `summarise`. A grouped object takes no `by=`. Rows come in key order, missing keys last.
- `count(*by, wt=..., sort=False, name="n")` — count rows, optionally weighted; returns a DataFrame with groups in key order, missing keys last (`sort=True` orders by count, ties keep key order). Use `summarize(n=ap.n())` for no key.
- `add_count(*by, wt=..., name="n")` — append a count column to AnnData and return AnnData.

**Joins:**
- `left_join(right, axis=..., by=..., relationship=...)` / `inner_join` / `right_join` / `full_join` / `semi_join` / `anti_join` — metadata joins on `obs` or `var`; reject duplicated or right-only cells/features. `by` stays a string naming a column in both tables.

**Extraction and pipeline:**
- `pull(name)` — extract exactly one column from any container, such as `pull(obs.louvain)`, `pull(raw.MS4A1)` or `pull(PCs[0])`.
- `to_wide(*names, allow_all_features=False, max_matrix_values=...)` — wide DataFrame, one row per observation: `to_wide(obs.louvain, raw.MS4A1, X_umap)`. A whole `obs` or `obsm` key is allowed; a whole `X`, `raw` or layer needs `allow_all_features=True`. Renamed from `to_df`, which `AnnData` already defines differently.
- `join_features(*features, shape="long", max_matrix_values=...)` — `shape="long"` gives a DataFrame of every `obs` column plus `.cell`, `.feature` and one `.abundance_<container>` column per container; `shape="wide"` returns the object with each feature added to `obs`.
- `pivot_longer(*names, names_to=..., values_to=..., allow_all_features=False, max_matrix_values=...)` — pivot selected metadata/matrix columns to long form.
- `aggregate_cells(*by, source=..., fun="sum")` — pseudobulk AnnData with one row per observed key combination and `n_cells` in `obs`. `source` is required because `X` is often scaled.
- `nest(*by, name="data")` — DataFrame with one row per group and that group's object in `data`; `ap.unnest(table)` stacks them back.
- `glimpse()` — print the `obs` head plus layer, `obsm` and `raw` names.
- `pipe(fn, *args, **kwargs)` — pass the AnnData through an arbitrary function and return the result.

## Expressions And Selectors

Use `obs.name` for metadata and `X.name`, `raw.name`, `counts.name` for matrix features. Use selectors for schema operations on `obs`:

```python
adata.select(ap.starts_with("qc_"))
adata.select(ap.matches("^qc_"))

adata.select(ap.everything())
adata.to_wide(ap.pick(ap.starts_with("qc_")))  # pick wraps a selector into a column list
```

Pick genes with a selector inside the container: `select(X[ap.any_of(panel)])` keeps the panel genes that exist.

Virtual axis names are available as `ap.obs_names` and `ap.var_names` in expressions.

Scoped helpers apply a predicate or function across selected columns:
- `ap.across(selector, fn, names="{col}_suffix")` — expand in `mutate`, `summarize`, or `transmute`.
- `ap.if_any(selector, fn)` / `ap.if_all(selector, fn)` — row-wise predicate for `filter` or `mutate`.

Aggregation helpers: `ap.n()`, `ap.mean(...)`, `ap.median(...)`, `ap.n_distinct(...)`, `ap.sum(...)`, `ap.min(...)`, `ap.max(...)`, `ap.sd(...)`.

Rank and cumulative helpers: `ap.row_number()`, `ap.min_rank(...)`, `ap.max_rank(...)`, `ap.dense_rank(...)`, `ap.average_rank(...)`, `ap.percent_rank(...)`, `ap.cume_dist(...)`, `ap.ntile(col, n)`, `ap.cummean(...)`, `ap.cumany(...)`, `ap.cumall(...)`, `ap.cum_max(...)`, `ap.cum_min(...)`, `ap.cum_prod(...)`.

Logic and recode helpers: `ap.near(col, value, tolerance=...)`, `ap.recode(col, mapping, default=...)`, `ap.case_match(col, (values, result), ...)`, `ap.coalesce(...)`, `ap.is_na(...)`, `ap.na_if(...)`, `ap.replace_na(...)`.

Offset and window helpers: `ap.lead(col, n=1)`, `ap.lag(col, n=1)`, `ap.first(col)`, `ap.last(col)`, `ap.nth(col, n)`.

## Joins

Joins enrich or subset `obs` or `var` metadata while preserving AnnData axes.

```python
joined = adata.left_join(sample_table, axis="obs", by="sample_id")
subset = adata.inner_join(allowed_genes, axis="var", by="gene_id")
```

Use relationship checks when cardinality matters. If a join would duplicate or add cells/features, expect `JoinRelationshipError` rather than silent axis corruption.

## Tidyr-Style Extraction

The pandas-only helpers live in `ap.tidyr`: `pivot_wider`, `nest`, `unnest`, `chop`, `unchop`, `pack`, `unpack`, `separate`, `separate_rows`, `extract`, `unite`, `drop_na`, `fill`, `hoist`, `unnest_longer` and `unnest_wider`. They take a DataFrame, not an AnnData. Keep them separate from AnnData-preserving verbs. Note that `ap.nest` and `ap.unnest` at the package root work on AnnData objects.

Prefer:

```python
plot_df = adata.join_features(X.MS4A1, X.CD79A, shape="long")
raw_plot_df = adata.join_features(raw.MS4A1, shape="long")
report = ap.tidyr.pivot_wider(plot_df, id_cols=[".cell"], names_from=".feature", values_from=".abundance_X")
```

Avoid unbounded matrix-to-long exports unless the user explicitly requests whole-matrix materialization. Use `max_matrix_values=` when the workflow needs a hard budget.

## scverse Utilities

Standalone helpers for single-cell–specific metadata patterns:

- `ap.feature_present(adata, features, feature_column=..., case_check=True)` — report requested feature presence.
- `ap.sample_meta(adata, "sample_id", include=..., exclude=...)` — extract one row of sample-constant metadata per sample. Attach a sample table with `left_join`, and summarize per sample with `summarize(..., by=obs.sample_id)`.
- `ap.rename_obs_names(adata, mapping_or_fn)` / `ap.rename_var_names(...)` — safe axis label edits; each returns a new object.
- `ap.get_palette(adata, "cell_type")` / `ap.store_palette(adata, "cell_type", palette)` — Scanpy-compatible palette access and storage.
- `ap.name_duplicates(adata, axis="obs")` / `ap.add_name_prefix(adata, prefix, axis="obs")` / `ap.replace_name_suffix(adata, current, new, axis="obs")` — name utilities.
