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

# Coming from scanpy

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

annplyr gives AnnData the dplyr verbs. Each task below runs twice on the same object, once as scanpy and pandas code and once as annplyr. A hidden check confirms that both give the same result. Then one ggann figure draws it.

Every task starts from this object. It is Scanpy's processed PBMC3k.

```python
import numpy as np
import pandas as pd
import scanpy as sc

import annplyr as ap
from annplyr.data import X, obs, raw

adata = ap.tidy(sc.datasets.pbmc3k_processed())
```

`ap.tidy` returns the same AnnData class with the verbs attached as methods. The `obs`, `raw` and `X` names say which container a column or gene comes from. Nothing is converted or copied.

## 1. Look at the object

**scanpy**

```{code-cell} ipython3
print(adata)
adata.obs.head()
```

**annplyr**

```{code-cell} ipython3
adata.glimpse()
```

```{code-cell} ipython3
ag.plot_embedding(adata, "umap", color="louvain", label=True)
```

*The eight Louvain clusters that `obs.louvain` lists, on the UMAP stored in `obsm`.*

`glimpse()` shows the shape, the first rows of `obs`, and which layers, embeddings and `raw` exist, in one call.

## 2. Filter cells on QC

**scanpy**

```{code-cell} ipython3
keep = (adata.obs["n_genes"] >= 500) & (adata.obs["percent_mito"] < 0.03)
qc_sc = adata[keep.to_numpy()].copy()
print(adata.n_obs - qc_sc.n_obs, "cells removed")
```

**annplyr**

```{code-cell} ipython3
qc_ap = adata.filter(obs.n_genes >= 500, obs.percent_mito < 0.03)
print(adata.n_obs - qc_ap.n_obs, "cells removed")
```

```{code-cell} ipython3
:tags: [remove-cell]

assert qc_sc.obs_names.equals(qc_ap.obs_names)
np.testing.assert_allclose(qc_sc.X, qc_ap.X)
np.testing.assert_allclose(qc_sc.obsm["X_umap"], qc_ap.obsm["X_umap"])
```

```{code-cell} ipython3
flagged = adata.mutate(
    qc=ap.if_else((obs.n_genes >= 500) & (obs.percent_mito < 0.03), "kept", "removed")
)
ag.plot_embedding(flagged, "umap", color="qc", palette={"kept": "#bbbbbb", "removed": "#d62728"})
```

```{code-cell} ipython3
:tags: [remove-cell]

assert (flagged.obs["qc"] == "kept").sum() == qc_ap.n_obs
```

*Red cells fail the filter. They fall in every large cluster, so the filter does not remove one cell type.*

`filter()` takes the predicates as arguments, joins them with AND, and keeps `X`, `obsm` and the graphs in step. The scanpy mask goes through `to_numpy()` because AnnData indexing wants an array, not a Series.

## 3. Pick a gene panel

The panel has two genes that `X` does not hold. `X` keeps only the highly variable genes.

**scanpy**

```{code-cell} ipython3
panel = ["MS4A1", "CD79A", "CD3E", "NKG7", "GNLY", "LYZ"]
have = [g for g in panel if g in adata.var_names]
panel_sc = adata[:, have].copy()
print(panel_sc.var_names.tolist())
```

**annplyr**

```{code-cell} ipython3
panel_ap = adata.select(X[ap.any_of(panel)])
print(panel_ap.var_names.tolist())
print(ap.feature_present(adata, panel).missing_features)
```

```{code-cell} ipython3
:tags: [remove-cell]

assert panel_sc.var_names.equals(panel_ap.var_names)
np.testing.assert_allclose(panel_sc.X, panel_ap.X)
assert ap.feature_present(adata, panel).missing_features == sorted(set(panel) - set(have))
```

```{code-cell} ipython3
ag.plot_dotplot(adata, panel, "louvain")
```

*Each dot is a cluster and a gene. Size is the share of cells expressing it.*

`any_of` drops absent genes without a message, like the list comprehension. `ap.feature_present` names what was dropped. The dot plot reads `raw`, so `CD3E` and `LYZ` are back there.

## 4. Put expression in a table

**scanpy**

```{code-cell} ipython3
wide_sc = sc.get.obs_df(adata, keys=["louvain", "n_genes", "MS4A1", "NKG7"], use_raw=True)
wide_sc.head()
```

**annplyr**

```{code-cell} ipython3
wide_ap = adata.to_wide(obs.louvain, obs.n_genes, raw.MS4A1, raw.NKG7)
wide_ap.head()
```

```{code-cell} ipython3
:tags: [remove-cell]

pd.testing.assert_frame_equal(wide_ap, wide_sc, check_dtype=False)
```

```{code-cell} ipython3
ag.plot_violin(adata, ["MS4A1", "NKG7"], "louvain")
```

*`MS4A1` marks B cells. `NKG7` marks NK and CD8 T cells.*

`sc.get.obs_df` needs `use_raw=True` to read `raw`. In annplyr the container is part of the name, so `raw.MS4A1` cannot be confused with `X.MS4A1`.

## 5. Add a score to `obs`

The score is the mean of two B cell genes.

**scanpy**

```{code-cell} ipython3
scored_sc = adata.copy()
b_expr = sc.get.obs_df(adata, keys=["MS4A1", "CD79A"], use_raw=True)
scored_sc.obs["b_score"] = b_expr.mean(axis=1)
scored_sc.obs.groupby("louvain", observed=True)["b_score"].mean().round(2)
```

**annplyr**

```{code-cell} ipython3
scored_ap = adata.mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)
scored_ap.summarize(b_score=ap.mean(obs.b_score), by=obs.louvain).round(2)
```

```{code-cell} ipython3
:tags: [remove-cell]

np.testing.assert_allclose(scored_sc.obs["b_score"], scored_ap.obs["b_score"], rtol=1e-6)
```

```{code-cell} ipython3
ag.plot_embedding(scored_ap, "umap", color="b_score")
```

*The score is high in one cluster only.*

`mutate()` returns a new object and leaves `adata` alone. The scanpy version needs an explicit `copy()` for the same guarantee. This is a plain mean. `sc.tl.score_genes` also subtracts a control gene set, so the two scores differ.

## 6. Summarize each cluster

**scanpy**

```{code-cell} ipython3
by_cluster = adata.obs.groupby("louvain", observed=True)
nkg7 = sc.get.obs_df(adata, keys=["NKG7"], use_raw=True)["NKG7"]
summary_sc = pd.DataFrame(
    {
        "cells": by_cluster.size(),
        "median_genes": by_cluster["n_genes"].median(),
        "mean_NKG7": nkg7.groupby(adata.obs["louvain"], observed=True).mean(),
    }
).reset_index()
summary_sc.round(2)
```

**annplyr**

```{code-cell} ipython3
summary_ap = adata.summarize(
    cells=ap.n(),
    median_genes=ap.median(obs.n_genes),
    mean_NKG7=ap.mean(raw.NKG7),
    by=obs.louvain,
)
summary_ap.round(2)
```

```{code-cell} ipython3
:tags: [remove-cell]

pd.testing.assert_frame_equal(summary_ap, summary_sc, check_dtype=False, check_exact=False)
```

```{code-cell} ipython3
import plotnine as p9

(
    p9.ggplot(summary_ap, p9.aes("louvain", "mean_NKG7"))
    + p9.geom_col()
    + p9.coord_flip()
    + p9.labs(x="", y="mean NKG7 (log-normalised)")
    + ag.theme_ggann()
    + p9.theme(figure_size=(4, 3))
)
```

*`summarize()` returns a pandas table, so a plain plotnine call draws it.*

Reductions read `obs` columns and genes in the same call, and the result keeps the cluster order stored on `obs.louvain`.

## 7. Count cells per group

PBMC3k is one sample. The donors below are synthetic: random labels, drawn once with a fixed seed.

```{code-cell} ipython3
rng = np.random.default_rng(0)
adata.obs["donor"] = pd.Categorical(rng.choice(["D1", "D2", "D3"], size=adata.n_obs))
```

**scanpy**

```{code-cell} ipython3
table_sc = pd.crosstab(adata.obs["donor"], adata.obs["louvain"])
table_sc
```

**annplyr**

```{code-cell} ipython3
table_ap = adata.count(obs.donor, obs.louvain)
table_ap.head()
```

```{code-cell} ipython3
:tags: [remove-cell]

stacked = table_sc.stack().reset_index(name="n")
pd.testing.assert_frame_equal(
    table_ap.reset_index(drop=True),
    stacked[stacked["n"] > 0].reset_index(drop=True),
    check_dtype=False,
)
```

```{code-cell} ipython3
ag.plot_proportions(adata, "louvain", "donor")
```

*Random donors give flat bars. Real donors would not.*

`count()` returns a long table with one row per observed pair. `pd.crosstab` returns a wide grid. The long form is what `ag.plot_proportions`, and every plotting grammar, wants.

## 8. Join a sample sheet

The sheet is invented. It gives each donor a condition and an age.

```{code-cell} ipython3
sheet = pd.DataFrame(
    {"donor": ["D1", "D2", "D3"], "condition": ["control", "control", "treated"], "age": [34, 51, 47]}
)
```

**scanpy**

```{code-cell} ipython3
merged = adata.obs.merge(sheet, on="donor", how="left")
print(merged.index[:2].tolist(), "<- the barcodes are gone")
joined_sc = adata.obs.join(sheet.set_index("donor"), on="donor", validate="m:1")
print(joined_sc.index.equals(adata.obs_names))
```

**annplyr**

```{code-cell} ipython3
joined_ap = adata.left_join(sheet, by=obs.donor, relationship="many-to-one")
print(joined_ap.obs_names[:2].tolist())
joined_ap.count(obs.donor, obs.condition, obs.age)
```

```{code-cell} ipython3
:tags: [remove-cell]

assert joined_ap.n_obs == adata.n_obs
assert joined_ap.obs_names.equals(adata.obs_names)
cols = ["donor", "condition", "age"]
pd.testing.assert_frame_equal(joined_ap.obs[cols].astype(str), joined_sc[cols].astype(str))
```

A repeated key in the sheet is the common mistake.

```{code-cell} ipython3
bad_sheet = pd.concat([sheet, sheet.iloc[[0]]])
print(len(adata.obs.merge(bad_sheet, on="donor")) - adata.n_obs, "extra rows from pandas")
try:
    adata.left_join(bad_sheet, by=obs.donor)
except ap.JoinRelationshipError as error:
    print(f"{type(error).__name__}: {error}")
```

```{code-cell} ipython3
:tags: [remove-cell]

assert len(adata.obs.merge(bad_sheet, on="donor")) > adata.n_obs
try:
    adata.left_join(bad_sheet, by=obs.donor)
except ap.JoinRelationshipError:
    pass
else:
    raise AssertionError("a repeated sheet key was accepted")
```

```{code-cell} ipython3
ag.plot_embedding(joined_ap, "umap", color="condition")
```

*The sheet's `condition` column now colours the cells.*

`left_join()` keeps the cell order and `obs_names`, and its default `relationship="many-to-one"` stops a repeated sheet row. The pandas merge adds rows silently, and an AnnData cannot hold them.

## 9. Relabel clusters

```{code-cell} ipython3
lineage = {
    "CD4 T cells": "T", "CD8 T cells": "T", "NK cells": "NK", "B cells": "B",
    "CD14+ Monocytes": "Myeloid", "FCGR3A+ Monocytes": "Myeloid",
    "Dendritic cells": "Myeloid", "Megakaryocytes": "Other",
}
```

**scanpy**

```{code-cell} ipython3
relabel_sc = adata.copy()
relabel_sc.obs["lineage"] = relabel_sc.obs["louvain"].astype(str).map(lineage).astype("category")
relabel_sc.obs["lineage"].value_counts()
```

**annplyr**

```{code-cell} ipython3
relabel_ap = adata.mutate(lineage=ap.recode(obs.louvain, lineage))
relabel_ap.count(obs.lineage)
```

```{code-cell} ipython3
:tags: [remove-cell]

assert (relabel_sc.obs["lineage"].astype(str) == relabel_ap.obs["lineage"].astype(str)).all()
```

```{code-cell} ipython3
ag.plot_embedding(relabel_ap, "umap", color="lineage", label=True)
```

*Eight clusters collapse to the lineages in the table above.*

`recode()` takes the mapping as data, so the same dict can come from a file. The result is a new `obs` column, and `louvain` stays as it was.

## 10. Keep the top cells per cluster

The goal is the three cells with the highest `n_counts` in each cluster.

**scanpy**

```{code-cell} ipython3
top_names = (
    adata.obs.sort_values("n_counts", ascending=False)
    .groupby("louvain", observed=True)
    .head(3)
    .index
)
top_sc = adata[top_names].copy()
top_sc.shape
```

**annplyr**

```{code-cell} ipython3
top_ap = adata.slice_max(obs.n_counts, n=3, by=obs.louvain)
top_ap.shape
```

```{code-cell} ipython3
:tags: [remove-cell]

assert sorted(top_sc.obs_names) == sorted(top_ap.obs_names)
```

```{code-cell} ipython3
(
    ag.gganndata(top_ap, p9.aes("louvain", "n_counts"))
    + p9.geom_point(size=2)
    + p9.coord_flip()
    + p9.labs(x="", y="UMI counts")
    + p9.theme(figure_size=(4, 3))
)
```

*Three points per cluster. `ag.gganndata` reads `louvain` and `n_counts` from the AnnData.*

`slice_max(by=)` is one call where pandas needs a sort, a groupby, a head and a second subset. The result is a smaller AnnData, so the embedding and `X` come along.

## 11. Make a long table for plotting

**scanpy**

```{code-cell} ipython3
genes = ["MS4A1", "NKG7"]
expr = sc.get.obs_df(adata, keys=genes, use_raw=True).rename_axis(".cell").reset_index()
long_sc = (
    expr.melt(id_vars=".cell", var_name=".feature", value_name=".abundance_raw")
    .merge(adata.obs, left_on=".cell", right_index=True, validate="m:1")
)
long_sc.head()
```

**annplyr**

```{code-cell} ipython3
long_ap = adata.join_features(raw.MS4A1, raw.NKG7, shape="long")
long_ap.head()
```

```{code-cell} ipython3
:tags: [remove-cell]

keys = [".cell", ".feature"]
a = long_ap.sort_values(keys).reset_index(drop=True)
b = long_sc[a.columns].sort_values(keys).reset_index(drop=True)
pd.testing.assert_frame_equal(a, b, check_dtype=False)
```

```{code-cell} ipython3
(
    p9.ggplot(long_ap, p9.aes("louvain", ".abundance_raw"))
    + p9.geom_boxplot(outlier_size=0.3)
    + p9.facet_wrap(".feature", ncol=1)
    + p9.labs(x="", y="log-normalised expression")
    + ag.theme_ggann()
    + p9.theme(axis_text_x=p9.element_text(rotation=45, ha="right"), figure_size=(4, 4))
)
```

*One panel per gene, drawn straight from the long table.*

The dotted column names (`.cell`, `.feature`, `.abundance_raw`) cannot clash with an `obs` column. Pass `shape="wide"` to add the genes to `obs` instead.

## 12. Use embedding coordinates as columns

```{code-cell} ipython3
from annplyr.data import X_umap
```

**scanpy**

```{code-cell} ipython3
coords = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"])
coords_sc = adata.obs[["louvain"]].join(coords)
coords_sc.head()
```

**annplyr**

```{code-cell} ipython3
coords_ap = adata.to_wide(obs.louvain, X_umap[0], X_umap[1])
coords_ap.head()
```

```{code-cell} ipython3
:tags: [remove-cell]

pd.testing.assert_frame_equal(coords_ap, coords_sc, check_dtype=False)
```

```{code-cell} ipython3
(
    p9.ggplot(coords_ap, p9.aes("UMAP_1", "UMAP_2", color="louvain"))
    + p9.geom_point(size=0.6)
    + ag.theme_ggann()
    + p9.theme(figure_size=(4.5, 3))
)
```

*The same UMAP as `ag.plot_embedding`, built by hand from the table.*

`X_umap[0]` is the first column of `obsm["X_umap"]`, and it is named `UMAP_1`.

## 13. Pseudobulk by donor and cluster

Pseudobulk needs raw counts. `X` is scaled, so the counts come from the unprocessed PBMC3k and are matched by barcode. The check fails if any cell has no counts.

```{code-cell} ipython3
from annplyr.data import counts

raw_counts = sc.datasets.pbmc3k()
if not adata.obs_names.isin(raw_counts.obs_names).all():
    raise ValueError("some cells have no raw counts")
adata.layers["counts"] = raw_counts[adata.obs_names, adata.var_names].X
```

**scanpy**

```{code-cell} ipython3
pb_sc = sc.get.aggregate(adata, by=["donor", "louvain"], func="sum", layer="counts")
pb_sc.shape
```

**annplyr**

```{code-cell} ipython3
pb_ap = adata.aggregate_cells(obs.donor, obs.louvain, source=counts)
pb_ap.obs.head()
```

```{code-cell} ipython3
:tags: [remove-cell]

key = lambda o: (o["donor"].astype(str) + "|" + o["louvain"].astype(str)).to_numpy()
order_sc = np.argsort(key(pb_sc.obs))
order_ap = np.argsort(key(pb_ap.obs))
assert (key(pb_sc.obs)[order_sc] == key(pb_ap.obs)[order_ap]).all()
sums_sc = np.asarray(pb_sc.layers["sum"].todense() if hasattr(pb_sc.layers["sum"], "todense") else pb_sc.layers["sum"])
sums_ap = np.asarray(pb_ap.X.todense() if hasattr(pb_ap.X, "todense") else pb_ap.X)
np.testing.assert_allclose(sums_sc[order_sc], sums_ap[order_ap])
```

```{code-cell} ipython3
lib = pb_ap.obs.assign(library_size=np.asarray(pb_ap.X.sum(axis=1)).ravel())
(
    p9.ggplot(lib, p9.aes("n_cells", "library_size", color="donor"))
    + p9.geom_point()
    + p9.scale_x_log10()
    + p9.scale_y_log10()
    + p9.labs(x="cells per pseudobulk", y="total counts")
    + ag.theme_ggann()
    + p9.theme(figure_size=(4, 3))
)
```

*Library size grows with the number of pooled cells. Each point is one donor and cluster.*

`aggregate_cells()` sums the layer you name in `source=`, one row per observed `donor` and `louvain` pair, and stores `n_cells` in `obs`. The donors here are random, so this shows the mechanics and not a biological result. The {doc}`Pseudobulk page <pseudobulk>` takes it through to a design table for a DE tool.

## What stays in scanpy

`sc.pp`, `sc.tl` and `sc.pl` are untouched. annplyr verbs return AnnData, so their results feed straight into scanpy.

```{code-cell} ipython3
sc.tl.rank_genes_groups(qc_ap, "louvain", method="t-test")
top = sc.get.rank_genes_groups_df(qc_ap, group=None).groupby("group", observed=True).head(1)
top[["group", "names", "scores"]].round(1)
```

```{code-cell} ipython3
:tags: [remove-cell]

import anndata

assert isinstance(qc_ap, anndata.AnnData)
assert top["group"].nunique() == qc_ap.obs["louvain"].nunique()
```

```{code-cell} ipython3
ag.plot_dotplot(qc_ap, top["names"].tolist(), "louvain")
```

*The top-ranked gene of each cluster on the QC-filtered object from task 2. A gene that tops two clusters appears once.*

Use annplyr for the table work around a scanpy workflow: filter, select, score, join, count, summarize. Use scanpy for normalisation, neighbours, embeddings and clustering.
