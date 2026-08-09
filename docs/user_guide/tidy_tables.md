# Tidy tables

Use extraction helpers at a deliberate boundary where the next step is pandas,
a plotting library, statistical modelling, reporting, or inspection outside
AnnData. Until that boundary, keeping transformations on AnnData preserves
alignment and avoids copying unused matrix values.

## Wide tables for models and inspection

`to_df()` returns one row per observation.

```{testcode}
wide = adata.ap.to_df(
    obs=["cell_type", "batch"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)

assert wide.shape == (4, 4)
```

Each selected feature becomes one column. Add embedding coordinates with a
keyed mapping when a model or interactive table needs them:

```python
model_frame = adata.ap.to_df(
    obs=["sample_id", "condition", "cell_type"],
    x=["MS4A1", "CD79A", "NKG7"],
    obsm={"X_pca": ["0", "1"]},
    max_matrix_values=5 * adata.n_obs,
)
```

## Long tables for faceted plots

`to_tidy()` returns stable observation, feature, and value columns.

```{testcode}
long = adata.ap.to_tidy(
    obs=["cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)

assert long.shape == (8, 4)
```

The stable `obs_name`, `feature`, and `value` columns map directly onto common
plotting grammars:

```python
import seaborn as sns

sns.catplot(
    data=long,
    x="cell_type",
    y="value",
    col="feature",
    kind="violin",
    sharey=False,
)
```

Whole-matrix long exports require explicit opt-in:

```{testcode}
long = adata.ap.to_tidy(
    allow_all_features=True,
    max_matrix_values=12,
)

assert len(long) == 12
```

Whole-matrix opt-in is useful for genuinely small assays. For ordinary
single-cell matrices, select the genes represented in the plot instead.
Budgets are checked cumulatively across every requested matrix source before
the first read.

## Grouped tables for cohort figures

Sometimes the figure needs one row per biological group rather than one row per
cell. Summarize before leaving AnnData:

```python
sample_summary = adata.ap.summarize(
    obs={"cells": ap.n()},
    x={
        "mean_MS4A1": ap.mean("MS4A1"),
        "mean_NKG7": ap.mean("NKG7"),
    },
    by=["sample_id", "condition", "cell_type"],
)
```

This produces a compact pandas DataFrame without first constructing a
cell-by-gene long table.

## Extract other AnnData containers

Use `as_frame()` for controlled access to AnnData containers:

```{testcode}
raw = adata.ap.as_frame("raw", select=["MS4A1"], max_matrix_values=4)
assert raw.shape == (4, 1)
```

The same interface addresses `obs`, `var`, `x`, `raw`, `obsm`, `varm`,
`obsp`, `varp`, and tabular `uns`. Pairwise and `uns` sources are extraction
only.

## Reshape extracted pandas data

These helpers work on plain DataFrames after extraction. annplyr provides
`pivot_wider`, `nest`, `unnest`, `chop`, `unchop`, `pack`, `unpack`,
`separate`, `separate_rows`, `extract`, `unite`, `drop_na`, and `fill`.
