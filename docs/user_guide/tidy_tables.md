# Tidy Tables

Use extraction helpers when the next step is pandas, plotting, reporting, or
inspection outside AnnData.

## Wide Tables

`to_df()` returns one row per observation.

```python
wide = adata.ap.to_df(
    obs=["cell_type", "sample"],
    x=["MS4A1", "CD79A"],
    obsm={"X_pca": ["0", "1"]},
)
```

## Long Tables

`to_tidy()` returns stable observation, feature, and value columns.

```python
long = adata.ap.to_tidy(
    obs=["cell_type"],
    x=["MS4A1", "CD79A"],
)
```

Whole-matrix long exports require explicit opt-in:

```python
long = adata.ap.to_tidy(
    allow_all_features=True,
    max_matrix_values=1_000_000,
)
```

## General Frame Extraction

Use `as_frame()` for controlled access to AnnData containers:

```python
raw = adata.ap.as_frame("raw", select=["MS4A1"])
neighbors = adata.ap.as_frame("obsp", key="connectivities", select=adata.obs_names[:10])
qc = adata.ap.as_frame("uns", key="qc_metrics")
```

## Pandas Rectangling Helpers

`annplyr` also provides pandas helpers such as `pivot_wider`, `nest`, `unnest`,
`chop`, `unchop`, `pack`, `unpack`, `separate`, `separate_rows`, `extract`,
`unite`, `drop_na`, and `fill`.
