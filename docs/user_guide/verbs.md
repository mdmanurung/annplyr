# Core verbs

The core verbs map familiar dataframe operations onto AnnData axes. Each call
does one recognizable transformation, and its return type tells you whether the
pipeline continues as AnnData, grouped AnnData, or pandas.

Every example runs against Scanpy's PBMC3K: 2,638 cells, 1,838 highly variable
genes in `X`, and 13,714 log-normalised genes in `.raw`. Names come from
`annplyr.data`, and `ap.tidy` makes the verbs methods of the object:

```{testcode}
from annplyr.data import X, obs, raw, var

adata = ap.tidy(adata)
```

## Choose the verb that matches the task

| Task | Verb |
|---|---|
| Keep cells or features matching predicates | `filter()` |
| Keep metadata columns or matrix features | `select()` |
| Add derived metadata | `mutate()` |
| Keep only newly derived metadata | `transmute()` |
| Reorder or take an axis subset | `arrange()`, `slice_*()` |
| Reduce groups to a table | `summarize()`, `count()` |
| Bring genes next to cell metadata | `join_features()` |
| Pseudobulk or split by group | `aggregate_cells()`, `nest()` |
| Pass the result to project code | `pipe()` |

Axis-changing verbs use positional AnnData slicing. Every verb returns a new
object.

## Filter across several sources at once

Predicates passed together are combined with `AND`, and each one names its own
source, so different sources can be mixed in the same call. Here metadata and expression are combined:

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
b_cells = adata.filter(
    obs.louvain == "B cells",
    obs.percent_mito < 0.03,
    raw.MS4A1 > 0,
)

print(b_cells.shape)
```

```{testoutput}
(250, 1838)
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
meta_mask = (adata.obs["louvain"] == "B cells") & (adata.obs["percent_mito"] < 0.03)
expr_mask = sc.get.obs_df(adata, keys=["MS4A1"], use_raw=True)["MS4A1"] > 0
b_cells_ref = adata[meta_mask & expr_mask.to_numpy()].copy()

print(b_cells_ref.shape)
```

```{testoutput}
(250, 1838)
```

:::

::::

Both baseline masks are pandas Series indexed by `obs_names`, but indexing an
AnnData with their conjunction wants a numpy array. Mixing metadata and
expression criteria is where hand-written masks go wrong.

Feature predicates use the `var` container, as in `var.n_cells > 2000`. A layer
imported from `annplyr.data`, such as `counts.MS4A1`, reads that layer instead
of `X`.

## Select the analysis schema

`select()` keeps metadata columns and matrix features in the requested order:

```{testcode}
markers = adata.select(
    obs.louvain,
    ap.starts_with("n_"),
    var.n_cells,
    X.MS4A1,
    X.CD79A,
    X.NKG7,
)

print(markers.shape)
print(markers.obs.columns.tolist())
print(markers.var_names.tolist())
```

```{testoutput}
(2638, 3)
['louvain', 'n_genes', 'n_counts']
['MS4A1', 'CD79A', 'NKG7']
```

One call subset both axes: three metadata columns and three of 1,838 genes.
Put a minus sign in front of an `obs` column to drop it: `select(-obs.louvain)`.
Use `all_of()` for strict configurable names, `any_of()` when absent names are
acceptable, and helpers such as `starts_with()`, `contains()`, or `where()` for
schema-driven selection of `obs` columns. Index a container with a selector to
pick features the same way: `select(X[ap.any_of(panel)])` keeps the panel genes
present in `X`, and `X[ap.starts_with("CD")]` keeps every gene starting with
`CD`. `any_of` drops absent genes without a message, and on processed data `X`
often holds only the highly variable genes, so check the panel first with
`ap.feature_present(adata, panel)`, or read `raw[...]` for every detected gene.

## Derive metadata with mutate

`mutate()` writes `obs` columns, or `var` columns with `into="var"`. Matrix,
layer, `raw`, embedding, and loading qualifiers are read-only sources for the
derived values:

```{testcode}
annotated = adata.mutate(
    deep=obs.n_counts >= 2000,
    b_markers=(raw.MS4A1 > 0) & (raw.CD79A > 0),
)

print(annotated.count(obs.louvain, obs.b_markers, sort=True).head())
```

```{testoutput}
           louvain  b_markers     n
0      CD4 T cells      False  1140
1  CD14+ Monocytes      False   480
2      CD8 T cells      False   315
3          B cells       True   277
4         NK cells      False   154
```

277 of 342 B cells co-express both markers, and essentially no other cell type
does. That is a two-line sanity check on the annotation.

Assignments are sequential, so a later expression may read a column created
earlier in the same call. `across()` applies one transformation to a selected
set:

```{testcode}
scaled = adata.mutate(
    scaled=ap.across(ap.starts_with("n_"), lambda name: ap.col(name) / 1000, names="{col}_k"),
)

print(scaled.obs.head(3).round(3))
```

```{testoutput}
                  n_genes  percent_mito  n_counts      louvain  n_genes_k  n_counts_k
index                                                                                
AAACATACAACCAC-1      781         0.030    2419.0  CD4 T cells      0.781       2.419
AAACATTGAGCTAC-1     1352         0.038    4903.0      B cells      1.352       4.903
AAACATTGATCAGC-1     1131         0.009    3147.0  CD4 T cells      1.131       3.147
```

`transmute()` instead returns an AnnData containing the grouping keys and the
newly produced columns only. New columns are keyword arguments. A dict is
accepted as a single positional argument for names that are not identifiers:
`mutate({"b-score": ...})`.

## Rename and relocate metadata

Rename takes the **new** name as the keyword and the existing column as the
value, like `rename()` in dplyr. A dict handles new names that are not
identifiers. `relocate()` takes names positionally, with `before=` or `after=`:

```{testcode}
renamed = adata.rename(cell_type=obs.louvain)

print(renamed.obs.columns.tolist())
print(renamed.relocate(obs.cell_type, before="n_genes").obs.columns.tolist())
```

```{testoutput}
['n_genes', 'percent_mito', 'n_counts', 'cell_type']
['cell_type', 'n_genes', 'percent_mito', 'n_counts']
```

`rename_with()` applies a function to selected names, which is useful for
normalizing imported sample-sheet columns.

## Arrange and slice

`arrange()` applies keys left to right and keeps ties in their original order.
A categorical column sorts by its category order, as in dplyr, not
alphabetically.

`slice_*()` works on either axis, so the same verb answers "deepest cells" and
"most widely detected genes":

```{testcode}
print(adata.slice_max(obs.n_counts, n=3).obs[["louvain", "n_counts"]])
print(adata.slice_max(var.n_cells, n=5, axis="var").var)
```

```{testoutput}
                          louvain  n_counts
index                                      
ACGAACTGGCTATG-1   Megakaryocytes    8875.0
GGGCCAACCTTGGA-1  Dendritic cells    8415.0
CAGGTTGAGGATCT-1          B cells    8011.0
        n_cells
index          
CYBA       2277
S100A6     2109
UBB        2094
CD37       2063
ARPC1B     2004
```

Grouped versions apply ordering and slicing within each group, which makes
"top ten cells per sample" a single expression.

## Summarize and count

`summarize()` reduces AnnData to a pandas DataFrame. Group keys can come from
`obs` or `var`, and reductions can read metadata and selected matrix features
in the same call:

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
summary = adata.summarize(
    cells=ap.n(),
    median_genes=ap.median(obs.n_genes),
    mean_mito=ap.mean(obs.percent_mito),
    mean_NKG7=ap.mean(raw.NKG7),
    by=obs.louvain,
)

print(summary.round(3))
```

```{testoutput}
             louvain  cells  median_genes  mean_mito  mean_NKG7
0        CD4 T cells   1144         809.0      0.019      0.119
1    CD14+ Monocytes    480         859.0      0.023      0.162
2            B cells    342         677.0      0.022      0.061
3        CD8 T cells    316         824.5      0.023      1.877
4           NK cells    154         890.0      0.021      3.102
5  FCGR3A+ Monocytes    150        1272.0      0.025      0.303
6    Dendritic cells     37        1544.0      0.020      0.303
7     Megakaryocytes     15         364.0      0.022      0.285
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
import pandas as pd

grouped = adata.obs.groupby("louvain", sort=True, observed=True)
nkg7 = sc.get.obs_df(adata, keys=["NKG7"], use_raw=True)["NKG7"]

summary_ref = pd.DataFrame(
    {
        "cells": grouped.size(),
        "median_genes": grouped["n_genes"].median(),
        "mean_mito": grouped["percent_mito"].mean(),
        "mean_NKG7": nkg7.groupby(adata.obs["louvain"], sort=True, observed=True).mean(),
    }
).reset_index()

print(summary_ref.round(3))
```

```{testoutput}
             louvain  cells  median_genes  mean_mito  mean_NKG7
0        CD4 T cells   1144         809.0      0.019      0.119
1    CD14+ Monocytes    480         859.0      0.023      0.162
2            B cells    342         677.0      0.022      0.061
3        CD8 T cells    316         824.5      0.023      1.877
4           NK cells    154         890.0      0.021      3.102
5  FCGR3A+ Monocytes    150        1272.0      0.025      0.303
6    Dendritic cells     37        1544.0      0.020      0.303
7     Megakaryocytes     15         364.0      0.022      0.285
```

:::

::::

`NKG7` concentrates in NK and CD8 T cells, which is the cytotoxic signature.

`count()` is the frequency table; `add_count()` writes the aligned counts back to
metadata. To count without a grouping key, use `summarize(n=ap.n())`:


```{testcode}
print(adata.add_count(obs.louvain).obs.head(3))
```

```{testoutput}
                  n_genes  percent_mito  n_counts      louvain     n
index                                                               
AAACATACAACCAC-1      781      0.030178    2419.0  CD4 T cells  1144
AAACATTGAGCTAC-1     1352      0.037936    4903.0      B cells   342
AAACATTGATCAGC-1     1131      0.008897    3147.0  CD4 T cells  1144
```

Canonical matrix summaries are chunked internally without changing dtype,
missing-value, or group-order semantics.

## Compose with pipe

`pipe()` passes the current object to a regular callable, which keeps the seam
between general wrangling and project-specific analysis explicit:

```python
def save_analysis(data, path):
    data.write_h5ad(path)
    return path


output = adata.filter(obs.percent_mito < 0.03).pipe(
    save_analysis,
    "results/qc_cells.h5ad",
)
```

Add a finite `max_matrix_values=` budget to any verb that reads matrix-backed
expressions when the projection size must be bounded. Invalid selectors,
sources, axes, sizes, budgets, and joins fail with the typed errors listed in
{doc}`../api`.
