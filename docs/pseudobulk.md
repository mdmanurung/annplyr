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

# Pseudobulk

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

`aggregate_cells` sums counts over donor and cell type and returns an AnnData with one row per pair. Donors are the replicates, so that table is what a differential-expression model should see. Cells are not replicates.

This page builds the table from PBMC3k, checks it, plots it, and shows the object you hand to a DE tool.

```{note}
PBMC3k comes from one donor. The six donors and the two conditions below are **synthetic**: cells get a random donor label and nothing else. Any result on this page shows how the table is built, not biology.
```

## Raw counts with processed labels

The processed PBMC3k object has labels but scaled values. The raw object has counts but no labels. Join them by cell name and check that every processed cell is found.

```{code-cell} ipython3
processed = sc.datasets.pbmc3k_processed()
cells = sc.datasets.pbmc3k()

# The check: every labelled cell exists in the raw object.
assert processed.obs_names.isin(cells.obs_names).all()

print(cells.n_obs, "raw cells;", processed.n_obs, "labelled cells")
cells = cells[processed.obs_names].copy()
cells.obs["louvain"] = processed.obs["louvain"]
cells.layers["counts"] = cells.X.copy()
counts = cells.layers["counts"]
print("integer counts:", bool((counts.data % 1 == 0).all()))
```

Cells that failed QC in the processed object are dropped. The counts are integers, so they are safe to sum.

## Add synthetic donors and condition

Donors go in as a random label. The sample sheet carries the condition, one row per donor. `relationship="many-to-one"` makes the join fail if a donor appears twice.

```{code-cell} ipython3
rng = np.random.default_rng(1)
donors = [f"D{i}" for i in range(1, 7)]
cells.obs["donor"] = pd.Categorical(rng.choice(donors, cells.n_obs), categories=donors)

sheet = pd.DataFrame({"donor": donors, "condition": ["control"] * 3 + ["treated"] * 3})
adata = ap.tidy(cells).left_join(sheet, by=obs.donor, relationship="many-to-one")
adata.count(obs.donor, obs.condition)
```

## Aggregate cells to pseudobulks

`source=layer("counts")` names the matrix to sum. `X` here is also counts, but naming the layer keeps the call right when `X` is scaled.

```{code-cell} ipython3
from annplyr import layer

pb = adata.aggregate_cells(obs.donor, obs.louvain, source=layer("counts"))
print(pb.shape)
pb.obs.head()
```

```{code-cell} ipython3
:tags: [remove-cell]

# scanpy's own aggregate gives the same sums
ref = sc.get.aggregate(adata, by=["donor", "louvain"], func="sum", layer="counts")
ref = ref[pb.obs_names].layers["sum"]
expected = pb.obs["donor"].astype(str).map(sheet.set_index("donor")["condition"])
assert (pb.obs["condition"].astype(str) == expected).all()
dense = lambda m: m.toarray() if hasattr(m, "toarray") else np.asarray(m)
np.testing.assert_allclose(dense(pb.X), dense(ref))
```

`aggregate_cells` keeps the obs columns that are constant within every group, so `condition` comes through with the donor. Small groups give noisy pseudobulks, so keep the ones built from at least 10 cells.

```{code-cell} ipython3
n_all = pb.n_obs
pb = pb.filter(obs.n_cells >= 10)
pb.obs["lib_size"] = np.asarray(pb.X.sum(axis=1)).ravel()
print(f"kept {pb.n_obs} of {n_all} pseudobulks;", pb.obs["n_cells"].min(), "to", pb.obs["n_cells"].max(), "cells each")
```

## Cells per sample

Every pseudobulk should rest on enough cells, and no donor should dominate a cell type.

```{code-cell} ipython3
from plotnine import aes, element_text, geom_col, labs, position_dodge, theme

(
    ag.gganndata(pb, aes("louvain", "n_cells", fill="donor"))
    + geom_col(position=position_dodge())
    + labs(x=None, y="cells in pseudobulk")
    + theme(figure_size=(5, 3), axis_text_x=element_text(rotation=45, ha="right"))
)
```

Look for cell types with short bars: those pseudobulks are the noisiest.

## Library size per pseudobulk

Summing makes library size track cell number. The model has to account for that, which is what normalization factors do.

```{code-cell} ipython3
from plotnine import geom_point, scale_y_log10

(
    ag.gganndata(pb, aes("n_cells", "lib_size", color="louvain"))
    + geom_point(size=2)
    + scale_y_log10()
    + labs(x="cells in pseudobulk", y="total counts")
    + theme(figure_size=(4, 3))
)
```

Look at the trend: total counts rise with cells per pseudobulk, and the cell types sit at different cell numbers.

## PCA of log-CPM pseudobulks

Scale to counts per million and take logs on a copy, so the counts stay intact. Then run PCA on the pseudobulks.

```{code-cell} ipython3
logcpm = pb.copy()
sc.pp.normalize_total(logcpm, target_sum=1e6)
sc.pp.log1p(logcpm)
sc.pp.highly_variable_genes(logcpm, n_top_genes=2000)
sc.pp.pca(logcpm, n_comps=5, mask_var="highly_variable")
print(logcpm.uns["pca"]["variance_ratio"].round(2))
```

```{code-cell} ipython3
from plotnine import theme

ag.compose(
    [
        ag.plot_embedding(logcpm, "pca", color="louvain", equal_aspect=False)
        + theme(aspect_ratio=0.8),
        ag.plot_embedding(logcpm, "pca", color="donor", equal_aspect=False)
        + theme(aspect_ratio=0.8),
    ],
    ncol=2,
) + theme(figure_size=(9, 3.5))
```

Left, cell type separates the samples. Right, the synthetic donors mix, as they should: donor labels were random.

## Hand-off for a DE tool

A DE tool wants a samples-by-genes count matrix and a table of sample metadata. Test one cell type at a time, with donor as the replicate and condition as the contrast.

```{code-cell} ipython3
b = pb.filter(obs.louvain == "B cells")
meta = b.obs[["donor", "condition", "n_cells"]]
mat = pd.DataFrame(b.X.toarray(), index=b.obs_names, columns=b.var_names)
mat = mat.loc[:, (mat > 0).sum(axis=0) >= 3]
print(mat.shape)
meta
```

Pass `mat` and `meta` to a pseudobulk DE package, for example with the design `~ condition`. This page does not run DE: with random donor labels there is nothing to find.
