# Tidy tables

Use extraction helpers at a deliberate boundary, where the next step is pandas,
a plotting library, statistical modelling, reporting, or inspection outside
AnnData. Until that boundary, keeping transformations on AnnData
preserves alignment and avoids copying matrix values nobody asked for.

The examples extract from Scanpy's PBMC3K. Its `X` holds scaled values, so the
log-normalised counts come from `raw=`.

## Wide tables for models and inspection

`to_df()` returns one row per observation, one column per requested value.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
wide = adata.ap.to_df(
    obs=["louvain", "n_genes"],
    raw=["MS4A1", "CD79A"],
    max_matrix_values=2 * adata.n_obs,
)

print(wide.head())
print(wide.shape)
```

```{testoutput}
                          louvain  n_genes  raw_MS4A1  raw_CD79A
index                                                           
AAACATACAACCAC-1      CD4 T cells      781        0.0        0.0
AAACATTGAGCTAC-1          B cells     1352    1.94591   1.386294
AAACATTGATCAGC-1      CD4 T cells     1131        0.0        0.0
AAACCGTGCTTCCG-1  CD14+ Monocytes      960        0.0        0.0
AAACCGTGTATGCG-1         NK cells      522        0.0        0.0
(2638, 4)
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
wide_ref = sc.get.obs_df(adata, keys=["MS4A1", "CD79A"], use_raw=True)
wide_ref.insert(0, "n_genes", adata.obs["n_genes"])
wide_ref.insert(0, "louvain", adata.obs["louvain"])

print(wide_ref.head())
```

```{testoutput}
                          louvain  n_genes    MS4A1     CD79A
index                                                        
AAACATACAACCAC-1      CD4 T cells      781  0.00000  0.000000
AAACATTGAGCTAC-1          B cells     1352  1.94591  1.386294
AAACATTGATCAGC-1      CD4 T cells     1131  0.00000  0.000000
AAACCGTGCTTCCG-1  CD14+ Monocytes      960  0.00000  0.000000
AAACCGTGTATGCG-1         NK cells      522  0.00000  0.000000
```

:::

::::

Two differences. annplyr prefixes matrix columns with their
source (`raw_MS4A1`), so a gene can never be confused with an `obs` column of
the same name. And its columns stay pandas *sparse*, which is why the zeros
print as `0.0`. A wide export of many genes stays cheap.

Embeddings come from the same call through a keyed mapping:

```{testcode}
model_frame = adata.ap.to_df(
    obs=["louvain"],
    obsm={"X_pca": ["0", "1"]},
    max_matrix_values=2 * adata.n_obs,
)

print(model_frame.head().round(3))
```

```{testoutput}
                          louvain  X_pca_0  X_pca_1
index                                              
AAACATACAACCAC-1      CD4 T cells    5.556    0.258
AAACATTGAGCTAC-1          B cells    7.210    7.482
AAACATTGATCAGC-1      CD4 T cells    2.694   -1.584
AAACCGTGCTTCCG-1  CD14+ Monocytes  -10.143   -1.369
AAACCGTGTATGCG-1         NK cells   -1.113   -8.153
```

## Long tables for faceted plots

`to_tidy()` returns stable `obs_name`, `feature`, and `value` columns, with the
requested metadata carried along.

```{testcode}
long = adata.ap.to_tidy(
    obs=["louvain"],
    raw=["MS4A1", "CD79A"],
    max_matrix_values=2 * adata.n_obs,
)

print(long.head())
print(long.shape)
```

```{testoutput}
           obs_name feature     value      louvain
0  AAACATACAACCAC-1   MS4A1       0.0  CD4 T cells
1  AAACATACAACCAC-1   CD79A       0.0  CD4 T cells
2  AAACATTGAGCTAC-1   MS4A1   1.94591      B cells
3  AAACATTGAGCTAC-1   CD79A  1.386294      B cells
4  AAACATTGATCAGC-1   MS4A1       0.0  CD4 T cells
(5276, 4)
```

Those column names map straight onto a plotting grammar:

```python
import seaborn as sns

sns.catplot(data=long, x="louvain", y="value", col="feature", kind="violin", sharey=False)
```

## Budgets are checked before the first read

`max_matrix_values=` bounds the whole request, cumulatively across every
requested matrix source, and the check happens **before** any adapter read:

```{testcode}
try:
    adata.ap.to_tidy(obs=["louvain"], raw=["MS4A1", "CD79A"], max_matrix_values=100)
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (to_tidy) would materialize 5276 matrix values, which exceeds max_matrix_values=100
```

The error arrives with the arithmetic already done, so the number is a real
count of what the request would have cost, not an estimate discovered halfway
through.

Whole-matrix long exports need an explicit opt-in *and* survive the same
budget. PBMC3K's `.raw` is 2,638 × 13,714, so the full melt is nearly five
million values:

```{testcode}
try:
    adata.ap.to_tidy(allow_all_features=True, max_matrix_values=10_000)
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (to_tidy) would materialize 4848644 matrix values, which exceeds max_matrix_values=10000
```

That number is what the opt-in exists to surface. For ordinary single-cell
matrices, select the genes the figure shows.

## Grouped tables for cohort figures

When the figure wants one row per group rather than one row per cell,
summarize before leaving AnnData. No cell-by-gene intermediate is built:

```{testcode}
print(
    adata.ap.summarize(
        obs={"cells": ap.n()},
        raw={"mean_MS4A1": ap.mean("MS4A1")},
        by="louvain",
    ).round(3)
)
```

```{testoutput}
             louvain  cells  mean_MS4A1
0        CD4 T cells   1144       0.034
1            B cells    342       0.993
2    CD14+ Monocytes    480       0.036
3           NK cells    154       0.036
4        CD8 T cells    316       0.040
5  FCGR3A+ Monocytes    150       0.060
6    Dendritic cells     37       0.056
7     Megakaryocytes     15       0.046
```

## Extract other AnnData containers

`as_frame()` gives controlled access to any container, including the ones that
have no verb-level shortcut, such as PCA loadings in `varm`:

```{testcode}
loadings = adata.ap.as_frame("varm", key="PCs", select=["0", "1"])

print(loadings.head().round(3))
print(loadings.shape)
```

```{testoutput}
             0      1
index                
TNFRSF4 -0.026  0.003
CPSF3L  -0.008  0.009
ATAD3C  -0.003  0.003
C1orf86  0.011 -0.000
RER1     0.014  0.027
(1838, 2)
```

The same interface addresses `obs`, `var`, `x`, `raw`, `obsm`, `varm`, `obsp`,
`varp`, and tabular `uns`. Pairwise and `uns` sources are extraction-only.

## Reshape extracted pandas data

These helpers work on plain DataFrames after extraction: `pivot_wider`, `nest`,
`unnest`, `chop`, `unchop`, `pack`, `unpack`, `separate`, `separate_rows`,
`extract`, `unite`, `drop_na`, and `fill`.
