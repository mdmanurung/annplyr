# Quickstart

Wrap the object once, import the containers you read, then call verbs:

```python
import annplyr as ap
from annplyr.data import X, obs, raw

adata = ap.tidy(adata)
filtered = adata.filter(
    obs.batch == "A",
    X.MS4A1 > 0,
)
```

## Naming the container

A column name always carries its container, because no default is safe:
`obs.n_genes` for metadata, `X.MS4A1` for `X`, `raw.MS4A1` for `.raw`,
`var.n_cells` for the feature axis, `counts.MS4A1` for a layer and `X_pca[0]` for
an `obsm` or `varm` key. Import layers and keys from `annplyr.data`, or use
`ap.layer("name")` for a layer whose name is reserved. Brackets cover names that
are not identifiers: `raw["HLA-DRA"]`. One expression cannot mix containers;
split it or precompute the shared part with `mutate` first. Strings such as
`"raw:MS4A1"` still work for names computed in a loop.

## Chaining

Every AnnData-returning verb returns a `TidyAnnData`, which is an `AnnData`
subclass that also carries the verbs, so the chain needs no `.ap`:

```python
summary = (
    adata.filter(obs.n_genes >= 500)
    .select(ap.starts_with("n_"), obs.cell_type, X.MS4A1, X.CD79A)
    .mutate(score=(raw.MS4A1 + raw.CD79A) / 2)
    .summarize(cells=ap.n(), mean_score=ap.mean(obs.score), by=obs.cell_type)
)
```

`isinstance(result, AnnData)` still holds, so Scanpy and AnnData keep working.
Use `result.to_anndata()` when a plain `AnnData` is required.

Each verb is also a plain function taking the AnnData object first
(`ap.filter(adata, ...)`), and `adata.ap.filter(...)` works on a plain AnnData.
All three spellings call the same implementation.

Common patterns:

```python
selected = adata.select(
    obs.sample,
    obs.cell_type,
    ap.starts_with("n_"),
    var.gene_ids,
    X.MS4A1,
    X.CD79A,
)

without_label = adata.select(-obs.cell_type)  # a minus sign drops an obs column

ordered = adata.arrange(ap.desc(obs.n_counts), obs.sample)

annotated = adata.mutate(
    high_counts=obs.n_counts > 10_000,
    ms4a1_positive=X.MS4A1 > 0,
)

renamed = adata.rename(cell_type=obs.louvain)
moved = renamed.relocate(obs.cell_type, before="n_genes")

raw_wide = adata.to_wide(obs.cell_type, raw.MS4A1)
umap = adata.to_wide(obs.cell_type, X_umap)  # a whole obsm key is allowed
first_pc = adata.pull(X_pca[0])
```

Use grouped summaries for observation metadata and selected matrix features:

```python
summary = adata.summarize(
    cells=ap.n(),
    mean_counts=ap.mean(obs.n_counts),
    mean_MS4A1=ap.mean(X.MS4A1),
    by=obs.cell_type,
)
```

Use grouped operations with `group_by`; most verbs operate group-locally on the result:

```python
# row number within each group; ungroup only at an AnnData boundary
ranked = adata.group_by(obs.batch).mutate(group_row=ap.row_number())
ranked_adata = ranked.ungroup()

# count observations per group, optionally weighted
adata.count(obs.batch)  # returns DataFrame
adata.count(obs.batch, wt=obs.n_counts, name="total_counts")
adata.add_count(obs.batch)  # appends the count column to AnnData
```

Use feature joins and extraction when the next step is pandas, plotnine, seaborn, or a notebook table:

```python
long = adata.join_features(raw.MS4A1, raw.CD79A, shape="long")  # .cell, .feature, .abundance_raw
with_genes = adata.join_features(raw.MS4A1, shape="wide")  # an object with MS4A1 in obs
wide = adata.to_wide(obs.cell_type, obs.sample, X.MS4A1, X.CD79A)
bounded = adata.to_wide(raw, allow_all_features=True, max_matrix_values=1_000_000)
```

Pool or split cells by group:

```python
pseudobulk = adata.aggregate_cells(obs.donor, obs.cell_type, source=counts)  # source is required
by_type = adata.nest(obs.cell_type)  # DataFrame with a `data` column of objects
stacked = ap.unnest(by_type)
adata.glimpse()  # obs head plus layer, obsm and raw names
```

Every verb returns a new object. Same-shape verbs share the matrices with the input; a subsetting verb on a backed input returns the view.
