# Quickstart

Importing annplyr registers the accessor:

```python
import annplyr as ap

filtered = adata.ap.filter(
    obs=ap.col("batch") == "A",
    x=ap.col("MS4A1") > 0,
)
```

Common patterns:

```python
selected = adata.ap.select(
    obs=ap.col("sample", "cell_type"),
    var=ap.starts_with("gene_"),
    x=ap.col("MS4A1", "CD79A"),
)

ordered = adata.ap.arrange(
    obs=ap.desc("n_counts"),
    var=ap.col("chromosome"),
)

annotated = adata.ap.mutate(
    obs={
        "high_counts": ap.col("n_counts") > 10_000,
        "ms4a1_positive": ap.col("MS4A1") > 0,
    },
    x={"MS4A1": ap.col("MS4A1")},
)

raw_wide = adata.ap.to_df(obs=["cell_type"], raw=["MS4A1"])
connectivity = adata.ap.as_frame("obsp", key="connectivities", select=adata.obs_names[:3])
```

Use grouped summaries for observation metadata and selected matrix sources:

```python
summary = adata.ap.summarize(
    obs={"cells": ap.n(), "mean_counts": ap.mean("n_counts")},
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by="cell_type",
)
```

Use grouped operations with `group_by`; most verbs operate group-locally on the result:

```python
# row number within each group
adata.ap.group_by(obs="batch").mutate(obs={"group_row": ap.row_number()})

# count observations per group, optionally weighted
adata.ap.count("batch")           # returns DataFrame
adata.ap.count("batch", wt="n_counts", name="total_counts")
adata.ap.add_count("batch")       # appends the count column to AnnData
```

Use plot-ready extraction when the next step is pandas, plotnine, seaborn, or a notebook table:

```python
wide = adata.ap.to_df(obs=["cell_type", "sample"], x=["MS4A1", "CD79A"])
long = adata.ap.to_tidy(obs=["cell_type"], x=["MS4A1", "CD79A"])
bounded = adata.ap.to_tidy(allow_all_features=True, max_matrix_values=1_000_000)
```
