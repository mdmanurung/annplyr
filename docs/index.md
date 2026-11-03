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

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

Stub: one-paragraph pitch, install, a 10-line example, and a "Start here" box pointing to
{doc}`from_scanpy`. Link the [changelog](https://github.com/mdmanurung/annplyr/blob/main/CHANGELOG.md).

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

```{code-cell} ipython3
ag.plot_embedding(adata, "umap", color="louvain")
```
