# Scverse Safety

`annplyr` is designed to make common AnnData wrangling safer rather than to
replace AnnData storage semantics.

## Sparse Matrices

Sparse matrix sources are converted to pandas sparse columns where possible.
Prefer selecting features explicitly before exporting to pandas.

```python
adata.ap.to_df(x=["MS4A1", "CD79A"])
```

## Backed AnnData

Read-only operations are supported where AnnData can provide the data. Mutating
verbs raise `AnnplyrError` for backed objects.

```python
adata = adata.to_memory()
adata = adata.ap.mutate(obs={"new": ap.col("old")})
```

## Raw

`raw` is an explicit read-only source:

```python
adata.ap.filter(raw=ap.col("MS4A1") > 0)
adata.ap.to_tidy(obs=["cell_type"], raw=["MS4A1"])
```

## Pairwise Matrices

`obsp` and `varp` are extraction-only sources through `as_frame()` or `pull()`.
They are not join keys or mutation targets.

```python
adata.ap.as_frame("obsp", key="connectivities", select=adata.obs_names[:5])
```

## Materialization Budgets

Use `max_matrix_values=` when exporting matrix data to pandas.

```python
adata.ap.to_tidy(
    allow_all_features=True,
    max_matrix_values=1_000_000,
)
```

This fails before materializing more matrix values than the budget allows.
