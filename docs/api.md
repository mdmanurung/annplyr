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

```{code-cell} ipython3
:tags: [remove-cell]

from _setup import *
```

```{code-cell} ipython3
ag.plot_embedding(adata, "umap", color="louvain")
```

Use this page for exact signatures and return types. For task-led examples,
start with {doc}`from_scanpy` or {doc}`tutorial`.

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

## Consumer Typing

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

## Expressions And Selectors

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

## Aggregation Helpers

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

## Dataframe Helpers In `annplyr.tidyr`

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

## Single-Cell Utility Helpers

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
