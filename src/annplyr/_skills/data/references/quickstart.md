# Quickstart

Importing annplyr registers the accessor:

```python
import annplyr as ap

filtered = adata.ap.filter(
    ap.col("batch") == "A",
    ap.col("x:MS4A1") > 0,
)
```

## Naming the container

A bare name is `obs` metadata. A matrix always needs a qualifier, because no
default is safe: `x:` for `X`, `raw:` for `.raw`, a layer name such as
`counts:`, an `obsm`/`varm` key such as `X_pca:`, and `var:`/`varm:` for the
feature axis. One expression cannot mix containers; split it or precompute the
shared part with `mutate` first.

## Chaining

Every AnnData-returning verb returns a `TidyAnnData`, which is an `AnnData`
subclass that also carries the verbs, so only the first call needs `.ap`:

```python
summary = (
    adata.ap.filter(ap.col("n_genes") >= 500)
    .select(obs=ap.everything(), x=panel)
    .mutate({"score": (ap.col("raw:MS4A1") + ap.col("raw:CD79A")) / 2})
    .summarize({"cells": ap.n(), "mean_score": ap.mean("score")}, by="cell_type")
)
```

`isinstance(result, AnnData)` still holds, so Scanpy and AnnData keep working.
Use `result.to_anndata()` when a plain `AnnData` is required.

Each verb is also a plain function taking the AnnData object first, so a
pipeline can start without the accessor at all. Both spellings call the same
implementation.

Common patterns:

```python
selected = adata.ap.select(
    obs=ap.col("sample", "cell_type"),
    var=ap.starts_with("gene_"),
    x=ap.col("MS4A1", "CD79A"),
)

ordered = adata.ap.arrange(ap.desc("n_counts"), ap.col("var:chromosome"))

annotated = adata.ap.mutate(
    {
        "high_counts": ap.col("n_counts") > 10_000,
        "ms4a1_positive": ap.col("x:MS4A1") > 0,
    }
)

raw_wide = adata.ap.to_wide(obs=["cell_type"], raw=["MS4A1"])
connectivity = adata.ap.as_frame("obsp", key="connectivities", select=adata.obs_names[:3])
```

Use grouped summaries for observation metadata and selected matrix sources:

```python
summary = adata.ap.summarize(
    {"cells": ap.n(), "mean_counts": ap.mean("n_counts"), "mean_MS4A1": ap.mean("x:MS4A1")},
    by="cell_type",
)
```

Use grouped operations with `group_by`; most verbs operate group-locally on the result:

```python
# row number within each group; ungroup only at an AnnData boundary
ranked = adata.ap.group_by(obs="batch").mutate(obs={"group_row": ap.row_number()})
ranked_adata = ranked.ungroup()

# count observations per group, optionally weighted
adata.ap.count("batch")  # returns DataFrame
adata.ap.count("batch", wt="n_counts", name="total_counts")
adata.ap.add_count("batch")  # appends the count column to AnnData
```

Use plot-ready extraction when the next step is pandas, plotnine, seaborn, or a notebook table:

```python
wide = adata.ap.to_wide(obs=["cell_type", "sample"], x=["MS4A1", "CD79A"])
long = adata.ap.to_tidy(obs=["cell_type"], x=["MS4A1", "CD79A"])
bounded = adata.ap.to_tidy(allow_all_features=True, max_matrix_values=1_000_000)
```

Axis-changing verbs return independent AnnData by default. `copy=False` may be a view or materialized result. Same-shape writers use `inplace=True` only when exact input identity is required.
