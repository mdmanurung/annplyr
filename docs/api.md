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

# API

This page lists every public name with its signature. For task-led examples, start with {doc}`from_scanpy` or {doc}`tutorial`.

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

| If you want to... | Start with... |
|---|---|
| wrap an object so verbs are methods | `tidy` |
| name a container in an expression | `annplyr.data`, `layer`, `Source` |
| subset or reorder AnnData | `filter`, `select`, `arrange`, `slice*` |
| derive aligned metadata | `mutate`, `transmute`, `add_count` |
| reuse group-local operations | `group_by` and `GroupedAnnData` |
| attach sample or feature tables | the six join methods |
| produce a pandas table | `summarize`, `count`, `to_wide`, `pivot_longer`, `pull`, `join_features` |
| pool or split cells by group | `aggregate_cells`, `nest`, `unnest` |
| reshape an extracted DataFrame | `annplyr.tidyr` |
| construct expressions or selectors | `col`, `all_of`, `any_of`, `across` |

## Verb families

Six families cover the verbs. Each example runs on PBMC3k and draws what the verb returned.

### Subset

`filter`, `select`, `arrange`, `distinct` and the `slice*` verbs choose cells or features. The result is an AnnData with `X`, `raw` and the embeddings cut to match.

```{code-cell} ipython3
sub = adata.filter(obs.louvain.is_in(["B cells", "NK cells"]), obs.n_genes > 500)
print(adata.shape, "->", sub.shape, "| umap:", sub.obsm["X_umap"].shape)
```

```{code-cell} ipython3
:tags: [ggann]

ag.plot_embedding(sub, "umap", color="louvain")
```

*Two clusters remain, and the embedding kept their coordinates.*

### Derive

`mutate`, `transmute`, `rename`, `relocate` and `add_count` write metadata. The matrices are read and never changed.

```{code-cell} ipython3
scored = adata.mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)
scored.obs[["louvain", "b_score"]].head(3).round(2)
```

```{code-cell} ipython3
:tags: [ggann]

ag.plot_embedding(scored, "umap", color="b_score")
```

*The new `b_score` column colours the map like any gene.*

### Summarize

`summarize`, `count` and `group_by` reduce cells to groups. The result is a pandas table, in key order.

```{code-cell} ipython3
by_cluster = adata.summarize(cells=ap.n(), median_genes=ap.median(obs.n_genes), by=obs.louvain)
by_cluster
```

```{code-cell} ipython3
:tags: [ggann]

from plotnine import aes, coord_flip, geom_col, ggplot

ggplot(by_cluster, aes("louvain", "cells")) + geom_col() + coord_flip() + ag.theme_ggann()
```

*Cells per cluster, from the `cells` column above.*

### Join

The six joins attach a table to `obs` or `var`. They raise `JoinRelationshipError` when a join would add or duplicate cells.

```{code-cell} ipython3
lineage = pd.DataFrame(
    {
        "louvain": adata.obs["louvain"].cat.categories,
        "lineage": ["T", "myeloid", "B", "T", "NK", "myeloid", "myeloid", "platelet"],
    }
)
joined = adata.left_join(lineage, by=obs.louvain, relationship="many-to-one")
joined.count(obs.lineage)
```

```{code-cell} ipython3
:tags: [ggann]

ag.plot_embedding(joined, "umap", color="lineage")
```

*The lineage came from the table, and every cell kept its place on the map.*

### Reshape

`join_features`, `to_wide`, `pivot_longer`, `pull` and `nest` turn AnnData into tables. `join_features` returns one row per cell and feature.

```{code-cell} ipython3
long = adata.join_features(raw.MS4A1, raw.CD79A, raw.NKG7, shape="long")
long = long.rename(columns={".feature": "gene", ".abundance_raw": "expr"})
mean_expr = long.groupby(["gene", "louvain"], observed=True)["expr"].mean().reset_index()
mean_expr.head(3).round(2)
```

```{code-cell} ipython3
:tags: [ggann]

from plotnine import element_text, geom_tile, theme

(
    ggplot(mean_expr, aes("louvain", "gene", fill="expr"))
    + geom_tile()
    + theme(axis_text_x=element_text(rotation=45, ha="right"))
    + ag.theme_ggann()
)
```

*Mean expression per cluster, straight from the long table.*

### Pseudobulk

`aggregate_cells` pools cells into one row per group. The processed PBMC3k `.raw` holds log values, so this example takes the mean. Sum true counts for differential expression, as on {doc}`pseudobulk`. The donors here are **synthetic**.

```{code-cell} ipython3
cells = ap.tidy(adata.copy())
cells.obs["donor"] = np.array(["d1", "d2", "d3"])[np.arange(cells.n_obs) % 3]
pb = cells.aggregate_cells(obs.donor, obs.louvain, source=raw, fun="mean")
print(pb.shape)
pb.obs.head(3)
```

```{code-cell} ipython3
:tags: [ggann]

from plotnine import geom_col

ggplot(pb.obs, aes("louvain", "n_cells", fill="donor")) + geom_col(position="dodge") + coord_flip() + ag.theme_ggann()
```

*Cells behind each pseudobulk sample. Megakaryocytes and dendritic cells have only a handful per donor, so those samples are noisy.*

## Accessor

The `adata.ap` accessor is registered when `annplyr` is imported, and
`ap.tidy(adata)` returns a `TidyAnnData` whose verbs are methods. AnnData-returning
verbs preserve axis alignment for `X`, layers, `raw`, `obsm`, `varm`, `obsp`,
and `varp` through positional slicing. Expression sources include metadata,
selected `X`/layer features, `raw`, `obsm`, and `varm`; `pull()` also reaches
`obsp`, `varp`, and tabular `uns` values.

Joins are metadata joins for `obs` or `var`. They may enrich or subset an axis,
but they raise `JoinRelationshipError` when a requested join would add or
duplicate cells/features. Long matrix exports materialize data into pandas and
therefore require explicit feature selection unless `allow_all_features=True`;
use `max_matrix_values=` on export helpers when a hard materialization budget is
needed. Every verb returns a new object. A subsetting verb on a backed input
returns the view; call `.to_memory()` before a metadata-changing verb.

## Container names

`annplyr.data` holds the container names that expressions and positional
selections use. `obs`, `var`, `X` and `raw` are fixed. Any other name is a layer
or an `obsm`/`varm` key, so `from annplyr.data import counts, X_umap` works.
`ap.layer("name")` reads a layer whose name is reserved or not importable.
Strings such as `"raw:MS4A1"` and `ap.col` remain available for computed names.

```{eval-rst}
.. autoclass:: annplyr.Source


.. autofunction:: annplyr.layer

.. autofunction:: annplyr.tidy
```

## Ownership, materialization, and failures by family

| Public family | Return and ownership | Materialization | Principal typed failures |
|---|---|---|---|
| `filter`, `select`, `arrange`, `distinct`, `slice*` | New independent AnnData; a backed input returns the view | Projects requested positions/sources | `SelectionError`, `UnknownColumnError`, `UnknownSourceError`, `IncompatibleAxisError`, `AnnplyrError` |
| six accessor joins | New AnnData; grouped calls return a grouped wrapper | Never manufactures right-only aligned records | `JoinRelationshipError`, `DuplicateNameError`, `IncompatibleAxisError` |
| `rename`, `rename_with`, `relocate`, `mutate`, `add_count` | New object that copies `obs`/`var`/`uns` and shares the matrices | Metadata-only writes; matrix sources are read-only | `SelectionError`, `DuplicateNameError`, `SizeMismatchError`, `AnnplyrError` |
| `transmute` | New object holding the keys and the new columns | Replaces metadata columns while reading projected sources | selection/source/size/budget errors |
| `aggregate_cells`, `unnest` | New AnnData (pseudobulk, or the stacked groups) | Reads the named `source` once per group | selection/source/axis/budget errors |
| grouped AnnData verbs | New grouped wrapper; `ungroup()` returns AnnData | One shared positional group plan; keys retained or updated | selection, axis, join, size, and budget errors |
| summaries, `count`, and `group_*` inspection | pandas DataFrame/list-like metadata; no AnnData mutation | Canonical matrix reducers use deterministic bounded chunks; metadata and opaque expressions read their planned sources | selection/source/axis/budget errors |
| `pull`, `to_wide`, `join_features`, `pivot_longer`, `nest` | pandas/Series, `shape="wide"` AnnData, or a table of group objects; input unchanged | Explicit pandas materialization with cumulative budgets | source, selection, name-repair, and budget errors |
| expressions, selectors, and aggregators | `AnnplyrExpr`, selector, or expansion metadata; input unchanged | Lazy until evaluated; raw Narwhals is opaque and fully charged | `SelectionError`, `UnknownColumnError`, `UnknownSourceError` |
| `annplyr.tidyr` helpers | New pandas objects; input unchanged | Operates after explicit tabular extraction | `SelectionError`, `DuplicateNameError`, `SizeMismatchError` |
| single-cell utilities | Readers return tables/reports; writers return a new object | Metadata/name/palette operations only | `SelectionError`, `DuplicateNameError`, `SizeMismatchError`, `JoinRelationshipError` |

All cumulative `max_matrix_values=` checks plan every source before the first
adapter read. The mechanically checked callable-level contract is recorded in
[`api-contract-v0.3.md`](https://github.com/mdmanurung/annplyr/blob/main/docs/development/api-contract-v0.3.md).

The reduction chunk target is an internal implementation detail, currently
25,165,824 logical values. It is deliberately not a public tuning flag.

```{eval-rst}
.. autoclass:: annplyr.AnnplyrAccessor
   :members:

.. autoclass:: annplyr.GroupedAnnData
   :members:

.. autoclass:: annplyr.TidyAnnData
   :members: to_anndata, wrap
```

`TidyAnnData` is the `AnnData` subclass that AnnData-returning verbs give back.
It carries every accessor verb as a method, so `adata.ap.filter(...).select(...)`
chains without re-entering `.ap`.

## Consumer typing

`AnnplyrAccessor`, `GroupedAnnData`, and `AnnplyrExpr` are exported from the
package root for annotations. Construct expressions with helpers such as
`col()`, and obtain grouped objects through `group_by()`; the exported classes
are types, not alternative public constructors.

AnnData namespaces are registered dynamically, so a static checker cannot
discover `AnnData.ap` from the upstream class definition. Cast once at a typed
boundary, after importing annplyr, and all AnnData-returning accessor methods
remain chainable:

```python
from typing import cast

from anndata import AnnData
import annplyr as ap
from annplyr.data import X, obs
from annplyr.typing import AnnDataWithAnnplyr


def select_markers(adata: AnnData) -> AnnDataWithAnnplyr:
    typed = cast(AnnDataWithAnnplyr, adata)
    return typed.ap.filter(obs.qc_pass).select(X.MS4A1, X.CD79A)
```

`AnnDataWithAnnplyr` is an alias of the concrete `TidyAnnData` class. The cast
does not copy or wrap the object. `isinstance` checks against `AnnData` still
hold. `group_by()` returns `GroupedAnnData` whenever keys are supplied, and
returns the typed AnnData when called with no keys.

The supported aliases live in `annplyr.typing`:

| Alias | Contract |
|---|---|
| `Selector` | strings, annplyr selectors/expressions, raw Narwhals expressions, or sequences of those values |
| `Expression` | a column name, `AnnplyrExpr`, or raw Narwhals expression |
| `Source` | the container literals `pull()` reads (`"obs"`, `"var"`, `"x"`, `"raw"`, `"obsm"`, and so on); not the `annplyr.Source` class |
| `SourceSelectors` | keyed selectors for aligned mapping sources |
| `JoinInput`, `JoinBy`, and join-option aliases | the exact dataframe/mapping and validated literal choices accepted by joins |
| `GroupedReturn` | typed AnnData or `GroupedAnnData` |

The generated [`public-typing-contract.md`](https://github.com/mdmanurung/annplyr/blob/main/docs/development/public-typing-contract.md) inventories every root
symbol, public method signature, overload, return annotation, and alias.

## Flat verb functions

`ap.filter(adata, ...)` is the same call as `adata.ap.filter(...)`.

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.add_count
   annplyr.aggregate_cells
   annplyr.anti_join
   annplyr.arrange
   annplyr.count
   annplyr.distinct
   annplyr.filter
   annplyr.full_join
   annplyr.glimpse
   annplyr.group_by
   annplyr.inner_join
   annplyr.join_features
   annplyr.left_join
   annplyr.mutate
   annplyr.nest
   annplyr.pipe
   annplyr.pivot_longer
   annplyr.pull
   annplyr.relocate
   annplyr.rename
   annplyr.rename_with
   annplyr.right_join
   annplyr.select
   annplyr.semi_join
   annplyr.slice
   annplyr.slice_head
   annplyr.slice_max
   annplyr.slice_min
   annplyr.slice_sample
   annplyr.slice_tail
   annplyr.summarise
   annplyr.summarize
   annplyr.to_wide
   annplyr.transmute
   annplyr.unnest
```

## Expressions and selectors

Expression helpers return `AnnplyrExpr`. Operators, namespaces, and method
chains propagate conservative dependency/cardinality metadata; call
`to_narwhals()` for explicit raw-Narwhals interoperation. `where()` is a
schema/dtype selector and must not inspect values.

```{eval-rst}
.. autoclass:: annplyr.AnnplyrExpr
   :members: to_narwhals
```

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.all_of
   annplyr.any_of
   annplyr.across
   annplyr.average_rank
   annplyr.between
   annplyr.case_when
   annplyr.case_match
   annplyr.col
   annplyr.coalesce
   annplyr.contains
   annplyr.cum_max
   annplyr.cum_min
   annplyr.cum_prod
   annplyr.cum_sum
   annplyr.cumall
   annplyr.cumany
   annplyr.cume_dist
   annplyr.cummean
   annplyr.dense_rank
   annplyr.desc
   annplyr.ends_with
   annplyr.everything
   annplyr.first
   annplyr.if_all
   annplyr.if_any
   annplyr.if_else
   annplyr.is_na
   annplyr.lag
   annplyr.last
   annplyr.last_col
   annplyr.lead
   annplyr.lit
   annplyr.matches
   annplyr.max_rank
   annplyr.min_rank
   annplyr.n_distinct
   annplyr.na_if
   annplyr.near
   annplyr.nth
   annplyr.ntile
   annplyr.num_range
   annplyr.obs_names
   annplyr.percent_rank
   annplyr.pick
   annplyr.recode
   annplyr.replace_na
   annplyr.row_number
   annplyr.starts_with
   annplyr.var_names
   annplyr.where
```

## Aggregation helpers

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.max
   annplyr.mean
   annplyr.median
   annplyr.min
   annplyr.n
   annplyr.sd
   annplyr.sum
```

## DataFrame helpers in `annplyr.tidyr`

These helpers reshape the pandas DataFrames that annplyr returns. They take a
DataFrame, not an AnnData, so they live in their own module: call
`ap.tidyr.pivot_wider(df, ...)`. The module's `nest` and `unnest` work on
DataFrame list columns. `ap.nest` and `ap.unnest` work on AnnData objects.

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.tidyr.chop
   annplyr.tidyr.drop_na
   annplyr.tidyr.extract
   annplyr.tidyr.fill
   annplyr.tidyr.hoist
   annplyr.tidyr.nest
   annplyr.tidyr.pack
   annplyr.tidyr.pivot_wider
   annplyr.tidyr.separate
   annplyr.tidyr.separate_rows
   annplyr.tidyr.unchop
   annplyr.tidyr.unite
   annplyr.tidyr.unnest
   annplyr.tidyr.unnest_longer
   annplyr.tidyr.unnest_wider
   annplyr.tidyr.unpack
```

## Single-cell utility helpers

These helpers cover narrow AnnData metadata ergonomics inspired by common
single-cell workflows. They do not implement biological QC scoring, plotting
wrappers, species gene registries, or object conversion workflows.

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.FeaturePresence
   annplyr.add_name_prefix
   annplyr.feature_present
   annplyr.get_palette
   annplyr.name_duplicates
   annplyr.rename_obs_names
   annplyr.rename_var_names
   annplyr.replace_name_suffix
   annplyr.sample_meta
   annplyr.store_palette

```

## Errors

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.AnnplyrError
   annplyr.SelectionError
   annplyr.UnknownColumnError
   annplyr.UnknownSourceError
   annplyr.DuplicateNameError
   annplyr.NameRepairError
   annplyr.IncompatibleAxisError
   annplyr.SizeMismatchError
   annplyr.JoinRelationshipError
```
