# Core Verbs

The core verbs mirror tidy dataframe workflows while respecting AnnData axes.

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

## Summarize And Count

```python
adata.ap.summarize(
    obs={"cells": ap.n()},
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by="cell_type",
)

adata.ap.count("cell_type", sort=True)
```
