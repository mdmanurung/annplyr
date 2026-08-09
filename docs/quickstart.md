# Quickstart

This page runs a complete wrangling pass over a real dataset: Scanpy's PBMC3K,
holding 2,638 peripheral blood mononuclear cells, 1,838 highly variable genes,
Louvain cell-type labels, and a `.raw` attribute with log-normalised counts for
all 13,714 detected genes.

Every example below executes during the documentation build, so every number on
this page is the number that code produces. Each step appears twice: once with
annplyr, once the way the same result is usually assembled today. Read the
second tab when you want to see what annplyr is saving you.

```{testcode}
import scanpy as sc

import annplyr as ap

adata = sc.datasets.pbmc3k_processed()
print(adata.obs["louvain"].value_counts())
```

```{testoutput}
louvain
CD4 T cells          1144
CD14+ Monocytes       480
B cells               342
CD8 T cells           316
NK cells              154
FCGR3A+ Monocytes     150
Dendritic cells        37
Megakaryocytes         15
Name: count, dtype: int64
```

```{note}
`pbmc3k_processed()` ships **scaled** values in `X`, the z-scores Scanpy's PCA
was run on. The interpretable log-normalised values live in `.raw`, so every
expression example on this page reads `raw=` rather than `x=`. annplyr makes you
state which one you mean.
```

## 1. Filter cells on QC metrics

Predicates in a tuple are combined with `AND`, so the call reads like the QC
rule itself and each threshold stays legible.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
qc = adata.ap.filter(
    obs=(
        ap.col("n_genes") >= 500,
        ap.col("percent_mito") < 0.03,
    )
)

print(qc.shape)
```

```{testoutput}
(2117, 1838)
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
mask = (adata.obs["n_genes"] >= 500) & (adata.obs["percent_mito"] < 0.03)
qc_ref = adata[mask].copy()

print(qc_ref.shape)
```

```{testoutput}
(2117, 1838)
```

:::

::::

The difference here is small but it compounds: you never build an intermediate
mask variable, never repeat `adata.obs[...]` once per threshold, and never have
to remember `.copy()`. Axis-changing verbs return an independent AnnData by
default, so `adata` is still intact for a second threshold.

## 2. Count what survived

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
print(qc.ap.count(by="louvain", sort=True))
```

```{testoutput}
             louvain    n
0        CD4 T cells  985
1    CD14+ Monocytes  349
2            B cells  267
3        CD8 T cells  229
4           NK cells  135
5  FCGR3A+ Monocytes  114
6    Dendritic cells   35
7     Megakaryocytes    3
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
counts_ref = qc_ref.obs["louvain"].value_counts().rename_axis("louvain").reset_index(name="n")

print(counts_ref)
```

```{testoutput}
             louvain    n
0        CD4 T cells  985
1    CD14+ Monocytes  349
2            B cells  267
3        CD8 T cells  229
4           NK cells  135
5  FCGR3A+ Monocytes  114
6    Dendritic cells   35
7     Megakaryocytes    3
```

:::

::::

Both give a tidy frame. `count()` gets there without the
`rename_axis`/`reset_index` dance, and takes the same `by=` selector every
other annplyr verb takes.

## 3. Score a marker pair into `obs`

This is where the gap widens. The score is a function of expression, but it
belongs in `obs`. The baseline has to pull a frame out of the matrix, do the
arithmetic, and write the column back with the cell order intact.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
scored = qc.ap.mutate(raw={"b_score": (ap.col("MS4A1") + ap.col("CD79A")) / 2})

print(scored.obs["b_score"].head().round(3))
```

```{testoutput}
index
AAACATTGATCAGC-1    0.000
AAACCGTGCTTCCG-1    0.000
AAACCGTGTATGCG-1    0.000
AAACGCACTGGTAC-1    0.000
AAACGCTGTAGCCA-1    0.347
Name: b_score, dtype: float32
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
expr = sc.get.obs_df(qc_ref, keys=["MS4A1", "CD79A"], use_raw=True)
qc_ref.obs["b_score"] = (expr["MS4A1"] + expr["CD79A"]) / 2

print(qc_ref.obs["b_score"].head().round(3))
```

```{testoutput}
index
AAACATTGATCAGC-1    0.000
AAACCGTGCTTCCG-1    0.000
AAACCGTGTATGCG-1    0.000
AAACGCACTGGTAC-1    0.000
AAACGCTGTAGCCA-1    0.347
Name: b_score, dtype: float32
```

:::

::::

The baseline mutates `qc_ref` in place; annplyr returns a new object and leaves
`qc` untouched. Only the two requested genes are read out of the 13,714 in
`.raw`. The matrix is a read-only source and is never densified.

## 4. Summarize metadata and expression together

One call, two sources, one grouping. The baseline needs a `groupby`, a separate
expression pull, a second `groupby` on the aligned series, and a manual
assembly. Every one of those steps is a chance to misalign a group.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
summary = scored.ap.summarize(
    obs={
        "cells": ap.n(),
        "median_genes": ap.median("n_genes"),
        "mean_b_score": ap.mean("b_score"),
    },
    raw={"mean_LYZ": ap.mean("LYZ")},
    by="louvain",
)

print(summary.round(3))
```

```{testoutput}
             louvain  cells  median_genes  mean_b_score  mean_LYZ
0        CD4 T cells    985         824.0         0.035     0.449
1    CD14+ Monocytes    349         923.0         0.035     3.747
2           NK cells    135         892.0         0.028     0.374
3        CD8 T cells    229         838.0         0.023     0.363
4  FCGR3A+ Monocytes    114        1272.0         0.049     2.227
5            B cells    267         713.0         1.254     0.403
6    Dendritic cells     35        1567.0         0.111     3.878
7     Megakaryocytes      3         984.0         0.116     1.030
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
import pandas as pd

lyz = sc.get.obs_df(qc_ref, keys=["LYZ"], use_raw=True)["LYZ"]
grouped = qc_ref.obs.groupby("louvain", sort=False, observed=True)

summary_ref = pd.DataFrame(
    {
        "cells": grouped.size(),
        "median_genes": grouped["n_genes"].median(),
        "mean_b_score": grouped["b_score"].mean(),
        "mean_LYZ": lyz.groupby(qc_ref.obs["louvain"], sort=False, observed=True).mean(),
    }
).reset_index()

print(summary_ref.round(3))
```

```{testoutput}
             louvain  cells  median_genes  mean_b_score  mean_LYZ
0        CD4 T cells    985         824.0         0.035     0.449
1    CD14+ Monocytes    349         923.0         0.035     3.747
2           NK cells    135         892.0         0.028     0.374
3        CD8 T cells    229         838.0         0.023     0.363
4  FCGR3A+ Monocytes    114        1272.0         0.049     2.227
5            B cells    267         713.0         1.254     0.403
6    Dendritic cells     35        1567.0         0.111     3.878
7     Megakaryocytes      3         984.0         0.116     1.030
```

:::

::::

B cells carry the marker score; `CD14+` monocytes and dendritic cells carry
`LYZ`.

When several steps share the same groups, group once and keep the grouping:

```{testcode}
by_type = scored.ap.group_by(obs="louvain")

print(by_type.summarize(raw={"mean_MS4A1": ap.mean("MS4A1")}).round(3))
```

```{testoutput}
             louvain  mean_MS4A1
0        CD4 T cells       0.036
1    CD14+ Monocytes       0.042
2           NK cells       0.041
3        CD8 T cells       0.031
4  FCGR3A+ Monocytes       0.061
5            B cells       1.016
6    Dendritic cells       0.059
7     Megakaryocytes       0.000
```

`MS4A1` is CD20, so an order of magnitude between B cells and everything else
is the expected answer.

## 5. Build the plot-ready long table

Plotting libraries want one row per cell per gene, with the metadata carried
along. `max_matrix_values=` bounds the whole request before any source is read,
so a mistyped selector fails instead of densifying 13,714 genes.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
tidy = scored.ap.to_tidy(
    obs=["louvain"],
    raw=["MS4A1", "CD79A", "LYZ", "NKG7"],
    max_matrix_values=4 * scored.n_obs,
)

print(tidy.head())
print(tidy.shape)
```

```{testoutput}
           obs_name feature     value          louvain
0  AAACATTGATCAGC-1   MS4A1       0.0      CD4 T cells
1  AAACATTGATCAGC-1   CD79A       0.0      CD4 T cells
2  AAACATTGATCAGC-1     LYZ  1.098612      CD4 T cells
3  AAACATTGATCAGC-1    NKG7       0.0      CD4 T cells
4  AAACCGTGCTTCCG-1   MS4A1       0.0  CD14+ Monocytes
(8468, 4)
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
long_ref = sc.get.obs_df(qc_ref, keys=["MS4A1", "CD79A", "LYZ", "NKG7"], use_raw=True)
long_ref["louvain"] = qc_ref.obs["louvain"]
long_ref = long_ref.reset_index(names="obs_name").melt(
    id_vars=["obs_name", "louvain"],
    var_name="feature",
    value_name="value",
)

print(long_ref.shape)
```

```{testoutput}
(8468, 4)
```

:::

::::

`tidy` is an ordinary pandas DataFrame, ready for seaborn, plotnine, Altair, or
a statistical model. `scored` is still a fully aligned AnnData object for
whatever comes next in the analysis.

```{note}
The `value` column comes back as a pandas sparse column, which keeps a large
export cheap but does not support every pandas method, `.round()` among them.
Call `tidy["value"] = tidy["value"].sparse.to_dense()` if you need a plain
float column.
```

Next: work through the {doc}`tutorials`, learn the storage and ownership rules
in {doc}`user_guide/concepts`, or look up an exact signature in {doc}`api`.
