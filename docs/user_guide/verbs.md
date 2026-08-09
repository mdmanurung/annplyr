# Core Verbs

The core verbs mirror tidy dataframe workflows while respecting AnnData axes.
Filtering, ordering, distinct selection, and slicing use integer positions and
return independent AnnData objects by default. Pass `copy=False` only when
either a view or a materialized result is acceptable.

## Filter

```python
adata.ap.filter(
    obs=ap.col("sample") == "s1",
    x=ap.col("MS4A1") > 0,
)
```

## Select

`select()` keeps metadata columns and can select features through `x=`.

```python
adata.ap.select(
    obs=ap.starts_with("qc_"),
    x=["MS4A1", "CD79A"],
)
```

## Rename And Relocate

`rename()` maps new names to existing names. `rename_with()` applies a function
to selected names.

```python
adata.ap.rename(obs={"sample_id": "sample"}, var={"gene_symbol": "symbol"})
adata.ap.rename_with(str.lower, obs=ap.starts_with("QC_"))
adata.ap.relocate(obs=["sample_id"], before="cell_type")
```

## Arrange And Slice

```python
ordered = adata.ap.arrange(obs=ap.desc("n_counts"))
top_genes = adata.ap.slice_head(n=20, axis="var")
```

## Mutate And Transmute

`mutate()` and `transmute()` write metadata columns only. They may read from
matrix-like sources but do not modify those matrices.

```python
adata.ap.mutate(
    obs={"high_counts": ap.col("n_counts") > 10_000},
    raw={"raw_MS4A1": ap.col("MS4A1")},
)
```

Mutation expressions are sequential: a later assignment may read a column
created earlier in the same call. annplyr batches only assignments proven to be
independent. Use `inplace=True` for a same-shape update that must return the
identical input object. `transmute()` always returns an independent object and
accepts no ownership flag.

## Summarize And Count

```python
adata.ap.summarize(
    obs={"cells": ap.n()},
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by="cell_type",
)

adata.ap.count("cell_type", sort=True)
```

Canonical scalar matrix summaries are automatically chunked across features or
rows. Chunking is internal and preserves assignment order, result dtype, NA
semantics, and first-seen grouping order; there is no public tuning argument.

Add `max_matrix_values=` whenever a verb reads matrix-backed expressions and a
hard cumulative projection limit is needed. Invalid selectors, sources, axes,
sizes, budgets, and joins raise the typed errors listed in {doc}`../api` before
the operation mutates its input.
