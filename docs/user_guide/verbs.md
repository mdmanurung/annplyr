# Core verbs

The core verbs map familiar dataframe operations onto AnnData axes. They are
small building blocks: each call does one recognizable transformation and its
return type determines whether the pipeline continues as AnnData, grouped
AnnData, or pandas.

## Choose the verb that matches the task

| Task | Verb |
|---|---|
| Keep cells or features matching predicates | `filter()` |
| Keep metadata columns or matrix features | `select()` |
| Add derived metadata | `mutate()` |
| Keep only newly derived metadata | `transmute()` |
| Reorder or take an axis subset | `arrange()`, `slice_*()` |
| Reduce groups to a table | `summarize()`, `count()`, `tally()` |
| Pass the result to project code | `pipe()` |

Axis-changing verbs use positional AnnData slicing and return independent
objects by default.

## Filter cells with several sources

Predicates in one tuple are combined with `AND`. Different source arguments
can be used in the same call:

```python
qc_b_cells = adata.ap.filter(
    obs=(
        ap.col("pct_counts_mt") < 10,
        ap.col("cell_type") == "B cell",
    ),
    x=ap.col("MS4A1") > 0,
)
```

The result contains only cells meeting the metadata and expression criteria,
while every layer, embedding, and pairwise matrix is subset consistently.
Feature predicates use `var=`; pass `layer="counts"` when `x=` expressions
should read a named layer instead of `X`.

## Select the analysis schema

`select()` keeps metadata columns and matrix features in the requested order:

```python
markers = adata.ap.select(
    obs=["sample_id", "condition", "cell_type", ap.starts_with("qc_")],
    var=["gene_symbol", "highly_variable"],
    x=["MS4A1", "CD79A", "CD3D", "NKG7"],
)
```

Use `all_of()` for strict configurable names, `any_of()` when absent names are
acceptable, and helpers such as `starts_with()`, `contains()`, or `where()` for
schema-driven selection.

## Derive metadata with mutate

`mutate()` writes `obs` or `var` columns. Matrix, layer, `raw`, embedding, and
loading arguments are read-only sources for aligned derived values.

```python
annotated = adata.ap.mutate(
    obs={
        "qc_pass": (ap.col("total_counts") >= 1_000)
        & (ap.col("pct_counts_mt") < 10),
    },
    x={
        "B_markers_detected": (ap.col("MS4A1") > 0)
        & (ap.col("CD79A") > 0),
    },
)
```

Assignments are sequential, so a later expression may read a column created
earlier in the same call. `across()` applies the same transformation to a
selected set:

```python
fractions = annotated.ap.mutate(
    obs=ap.across(
        ap.starts_with("pct_counts_"),
        lambda name: ap.col(name) / 100,
        names="{col}_fraction",
    )
)
```

Use `inplace=True` only when the same-shape metadata update should modify and
return the exact input object. `transmute()` instead returns an independent
AnnData containing the grouping keys and newly produced metadata columns.

## Rename and relocate metadata

Rename mappings use new names as keys and existing names as values:

```python
renamed = adata.ap.rename(
    obs={"sample_id": "sample"},
    var={"gene_symbol": "symbol"},
)

front_loaded = renamed.ap.relocate(
    obs=["sample_id", "condition"],
    before="cell_type",
)
```

`rename_with()` applies a function to selected names, which is useful for
normalizing imported sample-sheet columns.

## Arrange and slice

```python
ordered = adata.ap.arrange(obs=ap.desc("total_counts"))
highest_count_cells = adata.ap.slice_max(ap.col("total_counts"), n=100)
highly_variable_genes = adata.ap.slice_max(
    ap.col("variances"),
    n=2_000,
    axis="var",
)
```

Grouped versions apply ordering and slicing within each group, making patterns
such as “top ten cells per sample” explicit.

## Summarize and count

`summarize()` reduces AnnData to a pandas DataFrame. Group keys can come from
`obs` or `var`, and reductions can read metadata and selected matrix features
together:

```python
cell_type_summary = adata.ap.summarize(
    obs={
        "cells": ap.n(),
        "donors": ap.n_distinct("donor_id"),
        "mean_counts": ap.mean("total_counts"),
    },
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by=["condition", "cell_type"],
)
```

`count()` and `tally()` cover common frequency tables; `add_count()` and
`add_tally()` write aligned counts back to metadata. Canonical matrix summaries
are chunked internally without changing dtype, missing-value, or group-order
semantics.

## Compose with pipe

`pipe()` passes the current object to a regular callable. It is useful at the
boundary between general wrangling and project-specific analysis:

```python
def save_analysis(data, path):
    data.write_h5ad(path)
    return path


output = adata.ap.filter(obs=ap.col("qc_pass")).ap.pipe(
    save_analysis,
    "results/qc_cells.h5ad",
)
```

Add a finite `max_matrix_values=` budget to any verb that reads matrix-backed
expressions when the projection size must be bounded. Invalid selectors,
sources, axes, sizes, budgets, and joins fail with the typed errors listed in
{doc}`../api`.
