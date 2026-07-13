# API Patterns

## Accessor Verbs

Use `adata.ap.<verb>(...)` when the result should remain aligned AnnData.

- `filter(obs=..., var=..., x=..., layer=...)` subsets observations and features.
- `select(obs=..., var=..., x=...)` keeps metadata columns and selected features.
- `arrange(obs=..., var=...)` reorders observations or features.
- `slice*` helpers default to `axis="obs"`; pass `axis="var"` for features.
- `mutate(obs={...}, var={...}, x={...}, obsm={...}, varm={...})` writes metadata columns only.
- `group_by(obs=...)` returns a grouped wrapper for iteration, filtering, mutation, summaries, and counts.

## Expressions And Selectors

Use `ap.col("name")` for metadata or matrix columns in the current source. Use selectors for schema operations:

```python
adata.ap.select(obs=ap.starts_with("qc_"))
adata.ap.select(var=ap.matches("^MT-"))
adata.ap.select(obs=ap.everything())
```

Aggregation helpers such as `ap.n()`, `ap.mean("score")`, `ap.median(...)`, and `ap.n_distinct(...)` are intended for summaries.

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
```

Avoid unbounded matrix-to-long exports unless the user explicitly requests whole-matrix materialization.

