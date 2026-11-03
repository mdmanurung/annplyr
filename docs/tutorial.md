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

# Tutorial

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

This page takes PBMC3k from QC to a finished figure in six steps. Each step is one or two annplyr calls on an AnnData object, and each result is still an AnnData object or a plain table.

```{note}
PBMC3k is one donor. The three donors on this page are **synthetic**: each cell gets a random label. The cell types, QC values and expression are real.
```

## Look at the object

`glimpse` shows the shape and the first rows of `obs`, and lists what else the object holds.

```{code-cell} ipython3
adata.glimpse()
```

## Filter cells on QC

Name the container in front of each column: `obs.n_genes`. Predicates passed together must all hold. The result is a new object, and `adata` keeps every cell.

```{code-cell} ipython3
qc = adata.filter(obs.n_genes >= 500, obs.percent_mito < 0.03)
print(f"kept {qc.n_obs} of {adata.n_obs} cells")
```

```{code-cell} ipython3
:tags: [remove-cell]

keep = (adata.obs["n_genes"] >= 500) & (adata.obs["percent_mito"] < 0.03)
assert qc.obs_names.equals(adata.obs_names[keep])
```

## Check the marker panel

Check a panel before you use it. The processed `X` holds only the 1,838 highly variable genes, so `LYZ` is missing from it. It is in `raw`.

```{code-cell} ipython3
panel = ["MS4A1", "CD79A", "NKG7", "LYZ", "CD4"]
print(ap.feature_present(qc, panel).to_frame())
print(ap.feature_present(qc.raw.to_adata(), panel).to_frame()["found"].all())
```

The second line prints `True`: all five are in `raw`, so the rest of the page reads panel genes from `raw`.

## Score the B-cell cluster

A score is a new column in `obs`. `raw.MS4A1` reads one gene from the `raw` container, and the matrix itself is not changed.

```{code-cell} ipython3
scored = qc.mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)
scored.obs.groupby("louvain", observed=True)["b_score"].mean().round(2)
```

```{code-cell} ipython3
ag.plot_violin(scored, ["b_score"], group_by="louvain")
```

Look at the B cells: the score is high there, and the other clusters sit near zero.

## Composition by donor

Add the synthetic donor labels, then count cells per donor and cluster. `count` returns a table.

```{code-cell} ipython3
rng = np.random.default_rng(0)
scored.obs["donor"] = pd.Categorical(rng.choice(["D1", "D2", "D3"], scored.n_obs))
comp = scored.count(obs.donor, obs.louvain)
comp.head()
```

```{code-cell} ipython3
ag.plot_proportions(scored, group_by="louvain", split_by="donor")
```

Look at the segment heights across donors. Labels were random, so each donor should show a similar mix.

## Summaries per cluster

`summarize` with `by=` gives one row per cluster. Mix metadata and expression in one call.

```{code-cell} ipython3
summary = scored.summarize(
    cells=ap.n(),
    median_genes=ap.median(obs.n_genes),
    mean_b_score=ap.mean(obs.b_score),
    mean_LYZ=ap.mean(raw.LYZ),
    by=obs.louvain,
)
summary.round(2)
```

```{code-cell} ipython3
:tags: [remove-cell]

df = sc.get.obs_df(scored, keys=["n_genes", "b_score", "LYZ"], use_raw=True)
g = df.groupby(scored.obs["louvain"], observed=True)
ref = pd.DataFrame(
    {
        "cells": g.size(),
        "median_genes": g["n_genes"].median(),
        "mean_b_score": g["b_score"].mean(),
        "mean_LYZ": g["LYZ"].mean(),
    }
)
pd.testing.assert_frame_equal(
    summary.set_index("louvain").astype(float), ref.astype(float), check_names=False, check_dtype=False
)
```

## One figure

`ag.compose` arranges plots into a tagged grid. Every panel reads the same `scored` object.

```{code-cell} ipython3
from plotnine import theme

ag.compose(
    [
        ag.plot_embedding(scored, "umap", color="louvain")
        + theme(legend_position="none"),
        ag.plot_embedding(scored, "umap", color="b_score"),
        ag.plot_proportions(scored, group_by="louvain", split_by="donor"),
        ag.plot_violin(scored, ["LYZ"], group_by="louvain", use_raw=True, rotation=90)
        + theme(legend_position="none"),
    ],
    ncol=2,
) + theme(figure_size=(8, 6))
```

Top left, clusters on the UMAP (colours as in the key beside panel C). Top right, the B-cell score marks one cluster. Bottom left, composition per donor. Bottom right, `LYZ` is high in the monocyte clusters.

Next: [From scanpy](from_scanpy.md) maps each scanpy idiom to its annplyr call, and [Pseudobulk](pseudobulk.md) turns this object into donor-level counts.
