# Tidy tables

Use extraction helpers at a deliberate boundary, where the next step is pandas,
a plotting library, statistical modelling, reporting, or inspection outside
AnnData. Until that boundary, keeping transformations on AnnData
preserves alignment and avoids copying matrix values nobody asked for.

The examples extract from Scanpy's PBMC3K. Its `X` holds scaled values, so the
log-normalised counts come from `raw`.

```{testcode}
from annplyr.data import PCs, X_pca, obs, raw

adata = ap.tidy(adata)
```

## Wide tables for models and inspection

`to_wide()` returns one row per observation, one column per requested name.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
wide = adata.to_wide(
    obs.louvain,
    obs.n_genes,
    raw.MS4A1,
    raw.CD79A,
    max_matrix_values=2 * adata.n_obs,
)

print(wide.head())
print(wide.shape)
```

```{testoutput}
                          louvain  n_genes    MS4A1     CD79A
index                                                        
AAACATACAACCAC-1      CD4 T cells      781      0.0       0.0
AAACATTGAGCTAC-1          B cells     1352  1.94591  1.386294
AAACATTGATCAGC-1      CD4 T cells     1131      0.0       0.0
AAACCGTGCTTCCG-1  CD14+ Monocytes      960      0.0       0.0
AAACCGTGTATGCG-1         NK cells      522      0.0       0.0
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

Genes keep their own names, as in tidySingleCellExperiment, and a gene that
shares its name with an `obs` column raises `DuplicateNameError` instead of
being confused with it. The columns stay pandas *sparse*, which is why the zeros
print as `0.0`. A wide export of many genes stays cheap.

Embeddings come from the same call. `X_pca[0]` is the first component, named `PC_1` as in Seurat, tidySCE and ggann:

```{testcode}
model_frame = adata.to_wide(
    obs.louvain,
    X_pca[0],
    X_pca[1],
    max_matrix_values=2 * adata.n_obs,
)

print(model_frame.head().round(3))
```

```{testoutput}
                          louvain    PC_1   PC_2
index                                           
AAACATACAACCAC-1      CD4 T cells   5.556  0.258
AAACATTGAGCTAC-1          B cells   7.210  7.482
AAACATTGATCAGC-1      CD4 T cells   2.694 -1.584
AAACCGTGCTTCCG-1  CD14+ Monocytes -10.143 -1.369
AAACCGTGTATGCG-1         NK cells  -1.113 -8.153
```

## Long tables for faceted plots

`join_features()` returns every `obs` column plus `.cell`, `.feature` and one
`.abundance_<container>` column per container, so `raw` gives `.abundance_raw`.
The dotted names cannot collide with a metadata column.

```{testcode}
long = adata.join_features(
    raw.MS4A1,
    raw.CD79A,
    shape="long",
    max_matrix_values=2 * adata.n_obs,
)

print(long.drop(columns=["percent_mito", "n_counts"]).head())
print(long.shape)
```

```{testoutput}
              .cell .feature  .abundance_raw  n_genes      louvain
0  AAACATACAACCAC-1    MS4A1        0.000000      781  CD4 T cells
1  AAACATACAACCAC-1    CD79A        0.000000      781  CD4 T cells
2  AAACATTGAGCTAC-1    MS4A1        1.945910     1352      B cells
3  AAACATTGAGCTAC-1    CD79A        1.386294     1352      B cells
4  AAACATTGATCAGC-1    MS4A1        0.000000     1131  CD4 T cells
(5276, 7)
```

Those column names map straight onto a plotting grammar:

```python
import seaborn as sns

sns.catplot(data=long, x="louvain", y=".abundance_raw", col=".feature", kind="violin", sharey=False)
```

Pass `shape="wide"` to get the object back with the genes added to `obs`
instead, the form tidyseurat users know as `join_features(shape = "wide")`:

```{testcode}
with_genes = adata.join_features(raw.MS4A1, raw.CD79A, shape="wide")

print(with_genes.obs[["louvain", "MS4A1", "CD79A"]].head(3))
```

```{testoutput}
                      louvain    MS4A1     CD79A
index                                           
AAACATACAACCAC-1  CD4 T cells  0.00000  0.000000
AAACATTGAGCTAC-1      B cells  1.94591  1.386294
AAACATTGATCAGC-1  CD4 T cells  0.00000  0.000000
```


## Budgets are checked before the first read

`max_matrix_values=` bounds the whole request, cumulatively across every
requested matrix source, and the check happens **before** any adapter read:

```{testcode}
try:
    adata.join_features(raw.MS4A1, raw.CD79A, max_matrix_values=100)
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (join_features) would materialize 5276 matrix values, which exceeds max_matrix_values=100
```

The error arrives with the arithmetic already done, so the number is a real
count of what the request would have cost, not an estimate discovered halfway
through.

Naming a whole matrix needs an explicit opt-in *and* survives the same
budget. PBMC3K's `.raw` is 2,638 × 13,714, so the full melt is 36 million
values:

```{testcode}
try:
    adata.pivot_longer(obs.louvain, raw, allow_all_features=True, max_matrix_values=10_000)
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (pivot_longer) would materialize 36177532 matrix values, which exceeds max_matrix_values=10000
```

That number is what the opt-in exists to surface. For ordinary single-cell
matrices, select the genes the figure shows.

## Grouped tables for cohort figures

When the figure wants one row per group rather than one row per cell,
summarize before leaving AnnData. No cell-by-gene intermediate is built:

```{testcode}
print(
    adata.summarize(
        cells=ap.n(),
        mean_MS4A1=ap.mean(raw.MS4A1),
        by=obs.louvain,
    ).round(3)
)
```

```{testoutput}
             louvain  cells  mean_MS4A1
0        CD4 T cells   1144       0.034
1    CD14+ Monocytes    480       0.036
2            B cells    342       0.993
3        CD8 T cells    316       0.040
4           NK cells    154       0.036
5  FCGR3A+ Monocytes    150       0.060
6    Dendritic cells     37       0.056
7     Megakaryocytes     15       0.046
```

## Extract other AnnData containers

`pull()` reads one column from any container, including the ones that have no
verb-level shortcut, such as PCA loadings in `varm`. Join the columns you want
with pandas:

```{testcode}
import pandas as pd

loadings = pd.concat([adata.pull(PCs[0]), adata.pull(PCs[1])], axis=1)

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

`pull()` addresses `obs`, `var`, `X`, `raw`, `obsm`, `varm`, `obsp`, `varp`, and
tabular `uns`. Pairwise and `uns` sources are extraction-only.

## Split by group and pool cells

`nest()` returns one row per group, with that group's object in a `data`
column. `ap.unnest()` stacks the objects back into one AnnData:

```{testcode}
by_type = adata.nest(obs.louvain)

print(by_type["louvain"].tolist()[:3])
print(by_type["data"].iloc[0].shape)
print(ap.unnest(by_type).shape)
```

```{testoutput}
['CD4 T cells', 'B cells', 'CD14+ Monocytes']
(1144, 1838)
(2638, 1838)
```

`aggregate_cells()` sums a matrix over each group and returns a pseudobulk
AnnData: one row per group, the keys and `n_cells` in `obs`. Name the matrix
with `source=`, because `X` is often scaled:

```{testcode}
pseudobulk = adata.aggregate_cells(obs.louvain, source=raw)

print(pseudobulk.shape)
print(pseudobulk.obs["n_cells"].tolist())
```

```{testoutput}
(8, 13714)
[1144, 342, 480, 154, 316, 150, 37, 15]

```

## Reshape extracted pandas data

These helpers work on plain DataFrames after extraction. They live in
`ap.tidyr`: `pivot_wider`, `nest`, `unnest`, `chop`, `unchop`, `pack`, `unpack`,
`separate`, `separate_rows`, `extract`, `unite`, `drop_na`, `fill`, `hoist`,
`unnest_longer`, and `unnest_wider`. For example,
`ap.tidyr.pivot_wider(long, names_from=".feature", values_from=".abundance_raw")`.
