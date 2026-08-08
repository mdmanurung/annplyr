# Tidy Tables

Use extraction helpers when the next step is pandas, plotting, reporting, or
inspection outside AnnData.

## Wide Tables

`to_df()` returns one row per observation.

```{testcode}
wide = adata.ap.to_df(
    obs=["cell_type", "batch"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)

assert wide.shape == (4, 4)
```

## Long Tables

`to_tidy()` returns stable observation, feature, and value columns.

```{testcode}
long = adata.ap.to_tidy(
    obs=["cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)

assert long.shape == (8, 4)
```

Whole-matrix long exports require explicit opt-in:

```{testcode}
long = adata.ap.to_tidy(
    allow_all_features=True,
    max_matrix_values=12,
)

assert len(long) == 12
```

Budgets are checked cumulatively across every requested matrix source before
the first read. An opaque raw Narwhals expression receives a conservative
full-source charge.

## General Frame Extraction

Use `as_frame()` for controlled access to AnnData containers:

```{testcode}
raw = adata.ap.as_frame("raw", select=["MS4A1"], max_matrix_values=4)
assert raw.shape == (4, 1)
```

The same interface addresses `obs`, `var`, `x`, `raw`, `obsm`, `varm`,
`obsp`, `varp`, and tabular `uns`. Pairwise and `uns` sources are extraction
only.

## Pandas Rectangling Helpers

These work on plain DataFrames and complement the AnnData verbs above for
post-extraction wrangling.

`annplyr` also provides pandas helpers such as `pivot_wider`, `nest`, `unnest`,
`chop`, `unchop`, `pack`, `unpack`, `separate`, `separate_rows`, `extract`,
`unite`, `drop_na`, and `fill`.
