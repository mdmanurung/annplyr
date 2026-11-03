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

# Concepts

Five ideas explain how annplyr behaves: every name carries its container, verbs never change their input, matrix reads have a budget, failures have types, and groups come out in key order.

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

## Containers and namespaces

An AnnData holds several matrices that can share a gene name. annplyr makes you say which one you mean. `obs.n_genes` is a metadata column, `X.CST3` is a gene in `X`, and `raw.CST3` is the same gene in `.raw`.

The containers are importable names from `annplyr.data`. `obs`, `var`, `X` and `raw` are fixed. Any other name is a layer or an `obsm` or `varm` key, so `from annplyr.data import counts, X_umap` works. Brackets cover names that are not Python identifiers: `raw["HLA-DRA"]`.

The choice matters. In PBMC3k, `X` is scaled and centred, and `raw` holds log-normalized values. The same gene gives two different columns:

```{code-cell} ipython3
both = adata.mutate(CST3_raw=raw.CST3, CST3_X=X.CST3)
both.obs[["CST3_raw", "CST3_X"]].describe().loc[["mean", "min", "max"]].round(2)
```

```{code-cell} ipython3
:tags: [remove-cell]

expected = sc.get.obs_df(adata, "CST3", use_raw=True)["CST3"]
np.testing.assert_allclose(both.obs["CST3_raw"], expected)
np.testing.assert_allclose(both.obs["CST3_X"], sc.get.obs_df(adata, "CST3")["CST3"], rtol=1e-6)
```

```{code-cell} ipython3
:tags: [ggann]

panels = [
    ag.plot_embedding(both, "umap", color=name, title=title)
    for name, title in [("CST3_raw", "raw.CST3"), ("CST3_X", "X.CST3")]
]
ag.compose(panels, ncol=2)
```

*Same gene, same cells, two scales. `raw.CST3` is zero or positive. `X.CST3` is centred on zero, so low cells are negative.*

Pick the container that matches the claim you make. Use `raw` for expression you report, and `X` for what a model saw.

`X` can also hold fewer genes than `raw`:

```{code-cell} ipython3
print("genes in X:", adata.n_vars, "| genes in raw:", adata.raw.n_vars)
print("LYZ in X:", "LYZ" in adata.var_names, "| LYZ in raw:", "LYZ" in adata.raw.var_names)
```

LYZ is a classic monocyte marker and it lives only in `raw`. Asking `X` for it raises `UnknownColumnError`:

```{code-cell} ipython3
:tags: [raises-exception]

adata.mutate(LYZ=X.LYZ)
```

## Ownership

Every verb returns a new object, and the input is never changed. There is no `inplace=` or `copy=`. To keep a result, assign it.

A verb that keeps the shape, such as `mutate`, copies `obs`, `var` and `uns` and shares the matrices with its input. That is why it is cheap, and why editing a matrix in place shows in both objects.

```{code-cell} ipython3
flagged = adata.mutate(flag=ap.lit(1))
print("same object:", flagged is adata)
print("flag in input obs:", "flag" in adata.obs)
print("shares X:", np.shares_memory(flagged.X, adata.X))
```

A subsetting verb on a backed file returns a view. Call `.to_memory()` before a verb that changes metadata.

## Budgets

Table outputs copy matrix values into pandas. `max_matrix_values=` puts a ceiling on that copy. annplyr plans every source first, then refuses before it reads anything.

```{code-cell} ipython3
:tags: [raises-exception]

adata.join_features(raw.LYZ, raw.CD74, shape="long", max_matrix_values=100)
```

Reading a whole matrix needs the same explicit choice. `adata.to_wide(raw)` is refused until you name features or pass `allow_all_features=True`.

## Typed errors

Failures raise subclasses of `AnnplyrError`, so you can catch what you mean. The two above are `UnknownColumnError` for a name that is not there and a plain `AnnplyrError` for a budget. A join that would duplicate cells raises `JoinRelationshipError`:

```{code-cell} ipython3
:tags: [raises-exception]

sheet = pd.DataFrame({"louvain": ["B cells", "B cells"], "lineage": ["B", "B2"]})
adata.left_join(sheet, by=obs.louvain, relationship="many-to-one")
```

`relationship=` is checked before any cell is added. The full list of error classes is on the {doc}`api` page.

## Group order

Tables from `summarize` and `count` list groups in key order, as dplyr does. A categorical key follows its level order, and missing keys come last. The order does not depend on which cell comes first in the object.

`group_keys()` is different. It lists groups in the order they first appear:

```{code-cell} ipython3
by_type = adata.group_by(obs.louvain)
print(by_type.group_keys().louvain.tolist()[:4])
print(adata.count(obs.louvain).louvain.tolist()[:4])
print(adata.obs.louvain.cat.categories.tolist()[:4])
```
