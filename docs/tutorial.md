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

Stub: one end-to-end PBMC3k analysis.

```{code-cell} ipython3
ag.plot_embedding(adata, "umap", color="louvain")
```
