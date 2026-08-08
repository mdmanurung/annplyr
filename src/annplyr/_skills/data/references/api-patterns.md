# API Patterns

## Accessor Verbs

Use `adata.ap.<verb>(...)` when the result should remain aligned AnnData.

**Subsetting and column ops:**
- `filter(obs=..., var=..., x=..., raw=..., layer=..., copy=True, max_matrix_values=...)` — subset observations or features.
- `select(obs=..., var=..., x=..., copy=True)` — keep named metadata columns and selected features.
- `rename(obs={new_name: old_name}, var={...})` — rename metadata columns by mapping new names to existing names.
- `rename_with(str.lower, obs=...)` — rename selected columns by applying a function to names.
- `relocate(obs=..., before=..., after=...)` — reorder metadata columns. Raises `SelectionError` when the `before`/`after` anchor is among the columns being moved.
- `distinct(obs=..., var=..., keep_all=..., copy=True, max_matrix_values=...)` — deduplicate by selected columns.
- `arrange(obs=..., var=..., x=..., raw=..., copy=True, max_matrix_values=...)` — reorder observations or features.
- `slice(idx, copy=True)` / `slice_head(n=...)` / `slice_tail(n=...)` / `slice_min(...)` / `slice_max(...)` / `slice_sample(n=..., random_state=...)` — position or rank-based slicing; default `axis="obs"`, pass `axis="var"` for features.

**Mutation and summaries:**
- `mutate(obs={...}, var={...}, x={...}, raw={...}, obsm={...}, varm={...}, layer=..., inplace=False)` — write metadata columns; matrix sources are read-only.
- `transmute(obs={...}, ..., max_matrix_values=...)` — like `mutate` but drops all columns not produced by the call; always independent and has no ownership flag.
- `group_by(obs=..., var=...)` — return a persistent `GroupedAnnData` wrapper; call `ungroup()` for AnnData.
- `summarize(obs={...}, ..., by=..., max_matrix_values=...)` — aggregate by grouping key; also callable as `summarise`.
- `count(by=..., wt=..., sort=False, axis="obs", name="n")` — count rows, optionally weighted; returns a DataFrame.
- `tally(wt=..., sort=False, axis="obs", name="n")` — like `count` without a grouping argument.
- `add_count(by=..., wt=..., name="n")` — append a count column to AnnData and return AnnData.
- `add_tally(wt=..., name="n")` — like `add_count` without a grouping argument.

**Joins:**
- `left_join(right, axis=..., by=..., relationship=..., copy=True)` / `inner_join` / `right_join` / `full_join` / `semi_join` / `anti_join` — metadata joins on `obs` or `var`; reject duplicated or right-only cells/features.

**Extraction and pipeline:**
- `pull(obs=...)` / `pull(var=...)` / `pull(x=...)` / `pull(raw=...)` — extract exactly one requested source column.
- `to_df(obs=..., x=..., raw=..., obsm={...}, obsp={...}, max_matrix_values=...)` — wide DataFrame, one row per observation.
- `to_tidy(obs=..., x=..., raw=..., allow_all_features=False, max_matrix_values=...)` — long DataFrame with observation, feature, and value columns.
- `pivot_longer(obs=..., x=..., raw=..., names_to=..., values_to=..., max_matrix_values=...)` — pivot selected metadata/matrix columns to long form.
- `as_frame(source, key=..., select=..., max_matrix_values=...)` — controlled access to any container (`obs`, `var`, `x`, `raw`, `obsm`, `varm`, `obsp`, `varp`, or tabular `uns`).
- `nest_by(by=..., obs=..., var=..., axis="obs", name="data")` — nest axis subsets into AnnData objects grouped by metadata.
- `pipe(fn, *args, **kwargs)` — pass the AnnData through an arbitrary function and return the result.

## Expressions And Selectors

Use `ap.col("name")` for metadata or matrix columns in the current source. Use selectors for schema operations:

```python
adata.ap.select(obs=ap.starts_with("qc_"))
adata.ap.select(var=ap.matches("^MT-"))
adata.ap.select(obs=ap.everything())
adata.ap.to_df(obs=ap.pick(ap.starts_with("qc_")))  # pick wraps a selector into a column list
```

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
joined = adata.ap.left_join(sample_table, axis="obs", by="sample_id")
subset = adata.ap.inner_join(allowed_genes, axis="var", by="gene_id")
```

Use relationship checks when cardinality matters. If a join would duplicate or add cells/features, expect `JoinRelationshipError` rather than silent axis corruption.

## Tidyr-Style Extraction

Use `pivot_longer`, `pivot_wider`, `nest`, `unnest`, `separate`, `unite`, `drop_na`, and `fill` for pandas-style tables returned by annplyr extraction helpers. Keep these separate from AnnData-preserving verbs.

Prefer:

```python
plot_df = adata.ap.to_tidy(obs=["cluster"], x=["MS4A1", "CD79A"])
raw_plot_df = adata.ap.to_tidy(obs=["cluster"], raw=["MS4A1"])
```

Avoid unbounded matrix-to-long exports unless the user explicitly requests whole-matrix materialization. Use `max_matrix_values=` when the workflow needs a hard budget.

## scverse Utilities

Standalone helpers for single-cell–specific metadata patterns:

- `ap.feature_present(adata, features, feature_column=..., case_check=True)` — report requested feature presence.
- `ap.add_sample_meta(adata, meta, sample="sample_id", by=..., inplace=False)` — enrich `obs` with per-sample metadata.
- `ap.sample_meta(adata, "sample_id", include=..., exclude=...)` / `ap.sample_summary(adata, "sample_id", ...)` — extract or summarize sample metadata.
- `ap.rename_obs_names(adata, mapping_or_fn, inplace=False)` / `ap.rename_var_names(...)` — safe axis label edits.
- `ap.get_palette(adata, "cell_type")` / `ap.store_palette(adata, "cell_type", palette, inplace=False)` — Scanpy-compatible palette access and storage.
- `ap.name_duplicates(adata, axis="obs")` / `ap.add_name_prefix(adata, prefix, axis="obs", inplace=False)` / `ap.replace_name_suffix(adata, current, new, axis="obs", inplace=False)` — name utilities.
