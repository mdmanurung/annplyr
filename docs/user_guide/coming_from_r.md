# Coming from R

annplyr 0.5 follows the verbs of tidyseurat, tidySingleCellExperiment (tidySCE)
and tidytof. The one change from R is that every name carries its container:
`obs.percent_mito` where R writes `percent_mito`. Python has no bare names, so
the container prefix is the shortest honest spelling, and it removes the clash
between a gene and a metadata column of the same name.

## Syntax table

Column one is the annplyr 0.4 spelling that 0.5 replaces. Column two is the
spelling now. Column three is the R code for the same task.

| Task | annplyr 0.4 | annplyr now | tidyseurat / tidySCE / tidytof |
|---|---|---|---|
| Set up | `adata.ap.filter(...)` | `adata = ap.tidy(adata)`, then `adata.filter(...)` | `pbmc %>% filter(...)` |
| Print as table | `adata.obs.head()` | `adata.glimpse()` | printing `pbmc` shows a tibble |
| Filter on metadata | `.filter(ap.col("percent_mito") < 0.03)` | `.filter(obs.percent_mito < 0.03)` | `filter(percent_mito < 0.03)` |
| Filter on a gene | `.filter(ap.col("counts:MS4A1") > 0)` | `.filter(counts.MS4A1 > 0)` | `join_features("MS4A1", shape = "wide") %>% filter(MS4A1 > 0)` |
| Select | `.select(obs=["louvain", "n_genes"])` | `.select(obs.louvain, obs.n_genes)` | `select(louvain, n_genes)` |
| Drop a column | none | `.select(-obs.louvain)` | `select(-louvain)` |
| Select helpers | `.select(obs=ap.starts_with("n_"))` | `.select(ap.starts_with("n_"))` | `select(starts_with("n_"))` |
| Mutate | `.mutate({"b_score": (ap.col("raw:MS4A1") + ap.col("raw:CD79A")) / 2})` | `.mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)` | `join_features(...) %>% mutate(b_score = (MS4A1 + CD79A) / 2)` |
| if_else | `.mutate({"cls": ap.if_else(ap.col("louvain").is_in(M), "myeloid", "lymphoid")})` | `.mutate(cls=ap.if_else(obs.louvain.is_in(M), "myeloid", "lymphoid"))` | `mutate(cls = if_else(louvain %in% m, "myeloid", "lymphoid"))` |
| Rename | `.rename(obs={"cell_type": "louvain"})` | `.rename(cell_type=obs.louvain)` | `rename(cell_type = louvain)` |
| Arrange | `.arrange(ap.desc("n_genes"))` | `.arrange(ap.desc(obs.n_genes))` | `arrange(desc(n_genes))` |
| Count | `.count(by=["sample", "louvain"])` | `.count(obs.sample, obs.louvain)` | `count(sample, louvain)` |
| Add count | `.add_count(by="louvain")` | `.add_count(obs.louvain)` | `add_count(louvain)` |
| Grouped summary | `.group_by("louvain").summarize({"n": ap.n(), "ms4a1": ap.mean("raw:MS4A1")})` | `.group_by(obs.louvain).summarize(n=ap.n(), ms4a1=ap.mean(raw.MS4A1))` | `group_by(louvain) %>% summarise(n = n())` |
| Summary with `.by` | `.summarize({"n": ap.n()}, by="louvain")` | `.summarize(n=ap.n(), by=obs.louvain)` | `summarise(n = n(), .by = louvain)` |
| Top n per group | `.group_by("louvain").slice_max(ap.col("n_genes"), n=10)` | `.group_by(obs.louvain).slice_max(obs.n_genes, n=10)` | `group_by(louvain) %>% slice_max(n_genes, n = 10)` |
| Join metadata | `.left_join(meta, by="sample")` | unchanged | `left_join(meta, by = "sample")` |
| Genes, long table | `.to_tidy(obs=["louvain"], raw=["MS4A1", "LYZ"])` | `.join_features(raw.MS4A1, raw.LYZ, shape="long")`, giving `.cell`, `.feature`, `.abundance_raw` | `join_features(c("MS4A1", "LYZ"), shape = "long")` |
| Genes, wide | `.to_wide(obs=["louvain"], raw=["MS4A1"])`, a DataFrame | `.join_features(raw.MS4A1, shape="wide")`, an object with `MS4A1` in obs | `join_features("MS4A1", shape = "wide")` |
| Embeddings as columns | `.to_wide(obs=["louvain"], obsm={"X_umap": ...})` | `.to_wide(obs.louvain, X_umap)` | `as_tibble()` or `select(contains("UMAP"))` |
| Long metadata | `.pivot_longer(obs=["louvain", "leiden"], ...)` | `.pivot_longer(obs.louvain, obs.leiden, names_to=..., values_to=...)` | `pivot_longer(cols = c(louvain, leiden), ...)` |
| Nest per group | `.nest_by(by="cls")`, metadata only | `.nest(obs.cls)`, one object per group; `ap.unnest(df)` reverses it | `nest(data = -cls)` |
| Pseudobulk | none | `.aggregate_cells(obs.donor, obs.louvain, source=counts)`, giving an object | `aggregate_cells(.sample = c(donor, louvain), assays = "counts")` |
| Downsample (tidytof) | `.group_by("donor").slice_sample(n=100, random_state=0)` | `.group_by(obs.donor).slice_sample(n=100, random_state=0)` | `tof_downsample(group_cols = donor, num_cells = 100, method = "constant")` |

A whole container passed where names are expected means all of its columns.
`obs` and `obsm` keys such as `X_umap` are allowed as they are. Matrix
containers (`X`, `raw` and layers) need `allow_all_features=True`, because
`to_wide(raw)` would otherwise read all 13,714 genes of PBMC3K by accident.

## Where parity is impossible

- **Bare names.** Python has none, so the floor is a container prefix such as
  `obs.`. Strings such as `"percent_mito"` still work for names computed in a
  loop.
- **`.by` and `.sample`.** They become the `by=` keyword and positional keys.
- **`c(a, b)` and `%>%`.** They become positional arguments and method
  chaining.
- **Out of scope.** Clustering, dimension reduction and preprocessing stay with
  Scanpy.

## Try it

This block runs the R-style steps on PBMC3K. The names come from
`annplyr.data`:

```{testcode}
from annplyr.data import obs, raw

adata = ap.tidy(adata)
myeloid = ["CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic cells"]

labelled = (
    adata.filter(obs.percent_mito < 0.03)
    .mutate(cls=ap.if_else(obs.louvain.is_in(myeloid), "myeloid", "lymphoid"))
    .select(-obs.n_counts)
)

print(labelled.count(obs.cls, sort=True))
print(
    labelled.summarize(
        cells=ap.n(),
        mean_MS4A1=ap.mean(raw.MS4A1),
        by=obs.cls,
    ).round(3)
)
```

```{testoutput}
        cls     n
0  lymphoid  1712
1   myeloid   545
        cls  cells  mean_MS4A1
0  lymphoid   1712       0.199
1   myeloid    545       0.043
```

`join_features` with `shape="wide"` adds genes to `obs`, as tidySCE does:

```{testcode}
with_gene = labelled.join_features(raw.MS4A1, shape="wide")

print(with_gene.obs.columns.tolist())
print(with_gene.nest(obs.cls)["cls"].tolist())
```

```{testoutput}
['n_genes', 'percent_mito', 'louvain', 'cls', 'MS4A1']
['lymphoid', 'myeloid']
```

