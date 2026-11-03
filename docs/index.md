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

# annplyr

annplyr gives AnnData the verbs of dplyr: `filter`, `mutate`, `summarize`, `join`. Each verb takes a name that says where the data lives (`obs.percent_mito`, `raw.MS4A1`) and returns an AnnData that stays aligned, or a pandas table. Scanpy, scverse tools and your own code keep working on the result.

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

```bash
pip install annplyr
```

```{code-cell} ipython3
import annplyr as ap
from annplyr.data import obs, raw

adata = ap.tidy(sc.datasets.pbmc3k_processed())
scored = adata.filter(obs.percent_mito < 0.03).mutate(
    b_score=(raw.MS4A1 + raw.CD79A) / 2
)
scored.summarize(cells=ap.n(), b_score=ap.mean(obs.b_score), by=obs.louvain).round(2)
```

`filter` drops the cells above the mitochondrial cut-off and subsets `X`, `raw` and the embeddings with them. `mutate` writes the score to `obs`. `summarize` returns a pandas table. The plot reads the score straight from `obs`.

```{code-cell} ipython3
:tags: [ggann]

ag.plot_embedding(scored, "umap", color="b_score", title="B-cell score (mean of MS4A1 and CD79A)")
```

*One cluster lights up, and the table above shows B cells have the highest score.*

## Start here

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item-card} I know Scanpy
:link: from_scanpy
:link-type: doc

Each Scanpy task next to its annplyr form, on the same data.
:::

:::{grid-item-card} I know R
:link: from_r
:link-type: doc

tidySingleCellExperiment, tidyseurat and tidytof syntax next to annplyr.
:::

:::{grid-item-card} I know neither
:link: tutorial
:link-type: doc

One PBMC3k analysis from QC to a final figure.
:::
::::

Pseudobulk by donor is on {doc}`pseudobulk`. The ideas behind the API are on {doc}`concepts`, and every function is listed on {doc}`api`.

Read the [changelog](https://github.com/mdmanurung/annplyr/blob/main/CHANGELOG.md) for what changed between versions. To cite annplyr, use [CITATION.cff](https://github.com/mdmanurung/annplyr/blob/main/CITATION.cff) or GitHub's "Cite this repository" button.

```{toctree}
:hidden:
:maxdepth: 1

from_scanpy
from_r
concepts
tutorial
pseudobulk
api
```
