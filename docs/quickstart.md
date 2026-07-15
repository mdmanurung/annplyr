# Quickstart

`annplyr` exposes dataframe-style verbs on AnnData objects through `adata.ap`.

```python
import annplyr as ap

filtered = adata.ap.filter(
    obs=ap.col("batch") == "A",
    x=ap.col("MS4A1") > 0,
)
```

## Select Metadata And Features

```python
selected = adata.ap.select(
    obs=ap.col("sample", "cell_type"),
    var=ap.starts_with("gene_"),
    x=ap.col("MS4A1", "CD79A"),
)
```

## Add Metadata From Matrix Sources

`mutate()` writes only metadata columns; `x=`, `obsm=`, and other matrix
arguments provide read-only expression inputs.

```python
annotated = adata.ap.mutate(
    obs={"high_counts": ap.col("n_counts") > 10_000},
    x={"MS4A1": ap.col("MS4A1")},
)
```

## Summarize By Metadata

```python
summary = adata.ap.summarize(
    obs={"cells": ap.n(), "mean_counts": ap.mean("n_counts")},
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by="cell_type",
)
```

## Export Plot-Ready Tables

```python
wide = adata.ap.to_df(obs=["cell_type", "sample"], x=["MS4A1", "CD79A"])
long = adata.ap.to_tidy(obs=["cell_type"], x=["MS4A1", "CD79A"])
```

Next, read the {doc}`tutorials` or the full {doc}`notebooks/all_of_annplyr`
notebook.
