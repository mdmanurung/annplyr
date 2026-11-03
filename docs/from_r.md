---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Coming from R

annplyr follows the verbs of tidySingleCellExperiment (tidySCE), tidyseurat and tidytof. One thing changes: every name carries its container. Where R writes `percent_mito`, annplyr writes `obs.percent_mito`. Python has no bare names, and the prefix also keeps a gene and a metadata column of the same name apart.

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

## Syntax table

`adata` is `ap.tidy(adata)`. Names such as `obs`, `raw` and `counts` come from `annplyr.data`.

| Task | tidySCE, tidyseurat, tidytof | annplyr |
|---|---|---|
| Set up | `pbmc \|> filter(...)` | `adata = ap.tidy(adata)`, then `adata.filter(...)` |
| Print as a table | printing `pbmc` | `adata.glimpse()` |
| Filter on metadata | `filter(percent_mito < 0.03)` | `.filter(obs.percent_mito < 0.03)` |
| Filter on a gene | `join_features("MS4A1", shape = "wide") \|> filter(MS4A1 > 0)` | `.filter(raw.MS4A1 > 0)` |
| Select | `select(louvain, n_genes)` | `.select(obs.louvain, obs.n_genes)` |
| Drop a column | `select(-louvain)` | `.select(-obs.louvain)` |
| Select helpers | `select(starts_with("n_"))` | `.select(ap.starts_with("n_"))` |
| Gene with a dash | `` `HLA-DRA` `` | `raw["HLA-DRA"]` |
| Mutate | `mutate(b_score = (MS4A1 + CD79A) / 2)` | `.mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)` |
| if_else | `if_else(louvain %in% m, "myeloid", "lymphoid")` | `ap.if_else(obs.louvain.is_in(m), "myeloid", "lymphoid")` |
| Rename | `rename(cell_type = louvain)` | `.rename(cell_type=obs.louvain)` |
| Arrange | `arrange(desc(n_genes))` | `.arrange(ap.desc(obs.n_genes))` |
| Count | `count(donor, louvain)` | `.count(obs.donor, obs.louvain)` |
| Add count | `add_count(louvain)` | `.add_count(obs.louvain)` |
| Grouped summary | `group_by(louvain) \|> summarise(n = n())` | `.group_by(obs.louvain).summarize(n=ap.n())` |
| Summary with `.by` | `summarise(n = n(), .by = louvain)` | `.summarize(n=ap.n(), by=obs.louvain)` |
| Top n per group | `slice_max(n_genes, n = 10, by = louvain)` | `.slice_max(obs.n_genes, n=10, by=obs.louvain)` |
| Join a sample sheet | `left_join(sheet, by = "donor")` | `.left_join(sheet, by=obs.donor)` |
| Genes, long table | `join_features(c("MS4A1", "LYZ"), shape = "long")` | `.join_features(raw.MS4A1, raw.LYZ, shape="long")` |
| Genes, wide | `join_features("MS4A1", shape = "wide")` | `.join_features(raw.MS4A1, shape="wide")` |
| Embedding as columns | `as_tibble()` | `.to_wide(obs.louvain, X_umap)` |
| Long metadata | `pivot_longer(cols = c(louvain, leiden), ...)` | `.pivot_longer(obs.louvain, obs.leiden, names_to=..., values_to=...)` |
| Nest per group | `nest(data = -cls)` | `.nest(obs.cls)`; `ap.unnest(table)` reverses it |
| Pseudobulk | `aggregate_cells(.sample = c(donor, louvain), assays = "counts")` | `.aggregate_cells(obs.donor, obs.louvain, source=counts)` |
| Downsample (tidytof) | `tof_downsample(group_cols = donor, num_cells = 100, method = "constant")` | `.group_by(obs.donor).slice_sample(n=100, random_state=0)` |

Python cannot match R in four places:

- **Bare names.** The floor is a container prefix. Strings such as `"percent_mito"` still work for names built in a loop.
- **`.by` and `.sample`.** They become the `by=` keyword and positional keys.
- **`c(a, b)` and the pipe.** They become positional arguments and method chaining.
- **Preprocessing.** Clustering, dimension reduction and normalization stay in Scanpy.

## Three worked examples

The R code is shown as text. The annplyr code runs on PBMC3k.

### Label cells and filter

```r
pbmc |>
  filter(percent_mito < 0.03) |>
  mutate(cls = if_else(louvain %in% myeloid, "myeloid", "lymphoid"))
```

```{code-cell} ipython3
myeloid = ["CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic cells"]
labelled = adata.filter(obs.percent_mito < 0.03).mutate(
    cls=ap.if_else(obs.louvain.is_in(myeloid), "myeloid", "lymphoid")
)
labelled.count(obs.cls)
```

```{code-cell} ipython3
:tags: [remove-cell]

kept = adata.obs[adata.obs["percent_mito"] < 0.03]
expected = kept["louvain"].isin(myeloid).map({True: "myeloid", False: "lymphoid"})
pd.testing.assert_series_equal(
    labelled.count(obs.cls).set_index("cls")["n"],
    expected.value_counts().sort_index().rename("n").rename_axis("cls"),
    check_dtype=False,
)
```

```{code-cell} ipython3
:tags: [ggann]

ag.plot_embedding(labelled, "umap", color="cls", title="myeloid or lymphoid")
```

*Myeloid cells fill the left lobe. Lymphoid cells fill the other clusters.*

### Put a gene next to the metadata

```r
pbmc |>
  join_features(c("MS4A1", "CD79A"), shape = "wide") |>
  mutate(b_score = (MS4A1 + CD79A) / 2) |>
  ggplot(aes(louvain, b_score)) + geom_boxplot()
```

`join_features(shape="wide")` adds each gene to `obs` under its own name, as tidySCE does.

```{code-cell} ipython3
scored = adata.join_features(raw.MS4A1, raw.CD79A, shape="wide").mutate(
    b_score=(obs.MS4A1 + obs.CD79A) / 2
)
scored.obs[["MS4A1", "CD79A", "b_score"]].head(3).round(2)
```

```{code-cell} ipython3
:tags: [remove-cell]

by_hand = sc.get.obs_df(adata, ["MS4A1", "CD79A"], use_raw=True).mean(axis=1)
np.testing.assert_allclose(scored.obs["b_score"], by_hand, rtol=1e-5)
```

```{code-cell} ipython3
:tags: [ggann]

from plotnine import element_text, geom_boxplot, theme

(
    ag.gganndata(scored, ag.aes("louvain", "b_score"))
    + geom_boxplot(outlier_size=0.5)
    + theme(axis_text_x=element_text(rotation=45, ha="right"))
)
```

*The B-cell box sits well above the others.*

### Join a sample sheet and count

```r
pbmc |>
  left_join(sheet, by = "donor", relationship = "many-to-one") |>
  count(donor, louvain)
```

PBMC3k has one donor. This example adds three **synthetic** donors, assigned in turn to each cell, and a sample sheet of made-up ages.

```{code-cell} ipython3
cells = ap.tidy(adata.copy())
cells.obs["donor"] = np.array(["d1", "d2", "d3"])[np.arange(cells.n_obs) % 3]
sheet = pd.DataFrame({"donor": ["d1", "d2", "d3"], "age": [31, 47, 62]})

with_age = cells.left_join(sheet, by=obs.donor, relationship="many-to-one")
with_age.count(obs.donor, obs.age)
```

```{code-cell} ipython3
:tags: [remove-cell]

merged = cells.obs.merge(sheet, on="donor", how="left", validate="m:1")
assert (merged["age"].to_numpy() == with_age.obs["age"].to_numpy()).all()
```

```{code-cell} ipython3
:tags: [ggann]

ag.plot_proportions(with_age, "louvain", split_by="donor")
```

*Each donor carries nearly the same cluster mix, because the donors are assigned in turn and not drawn from real people.*

## R parity

Ten tasks from the tidySCE vignette give the same answer in tidySCE and annplyr. The test used the PBMC3k `.raw` matrix with synthetic donors, R 4.5.2 and Bioconductor 3.22. All ten match at `rtol` 1e-6:

| Task | What matched |
|---|---|
| filter | 2628 cells kept by both, same order |
| count | 24 groups, same row order |
| mutate | largest absolute difference 1.08e-07 |
| summarise | n, median and mean over 8 clusters |
| `slice_max(by=)` | 24 cells, same order |
| `join_features` long | 7914 cell by feature rows |
| rename and `if_else` | lymphoid 1971, myeloid 667 |
| `left_join` | 2638 cells |
| `aggregate_cells` | 24 donor by cluster groups, largest relative difference 4.4e-15 |
| `arrange` | first 50 cells, same order |

These numbers were recorded when the scripts last ran. The build does not rerun them, because R is not part of the docs environment. To rerun, from the repository root:

```bash
uvx hatch -e integration run python scripts/r_parity/01_annplyr.py
cd scripts/r_parity/out && Rscript ../02_tidysce.R && cd ../../..
uvx hatch -e integration run python scripts/r_parity/03_compare.py
```

The last step writes `scripts/r_parity/out/parity_summary.csv`.
