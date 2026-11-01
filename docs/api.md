# API

Use this page for exact signatures and return types. For task-led examples,
start with {doc}`quickstart` or {doc}`tutorials`.

| If you want to... | Start with... |
|---|---|
| subset or reorder AnnData | `filter`, `select`, `arrange`, `slice*` |
| derive aligned metadata | `mutate`, `transmute`, `add_count` |
| reuse group-local operations | `group_by` and `GroupedAnnData` |
| attach sample or feature tables | the six join methods |
| produce a pandas table | `summarize`, `to_wide`, `to_tidy`, `as_frame` |
| construct expressions or selectors | `col`, `all_of`, `any_of`, `across` |

## Accessor

The `adata.ap` accessor is registered when `annplyr` is imported. AnnData-returning
verbs preserve axis alignment for `X`, layers, `raw`, `obsm`, `varm`, `obsp`,
and `varp` through positional slicing. Expression sources include metadata,
selected `X`/layer features, `raw`, `obsm`, and `varm`; controlled extraction
also exposes `obsp`, `varp`, and tabular `uns` values through `as_frame()`.

Joins are metadata joins for `obs` or `var`. They may enrich or subset an axis,
but they raise `JoinRelationshipError` when a requested join would add or
duplicate cells/features. Long matrix exports materialize data into pandas and
therefore require explicit feature selection unless `allow_all_features=True`;
use `max_matrix_values=` on export helpers when a hard materialization budget is
needed. Backed `copy=True` subsetting returns an independent in-memory object;
call `.to_memory()` before same-shape mutation.

## Ownership, materialization, and failures by family

| Public family | Return and ownership | Materialization | Principal typed failures |
|---|---|---|---|
| `filter`, `select`, `arrange`, `distinct`, `slice*` | Independent AnnData by default; `copy=False` may be a view or materialized | Projects requested positions/sources; backed `copy=True` loads the selection | `SelectionError`, `UnknownColumnError`, `UnknownSourceError`, `IncompatibleAxisError`, `AnnplyrError` |
| six accessor joins | Independent AnnData by default; grouped calls return a grouped wrapper | Never manufactures right-only aligned records | `JoinRelationshipError`, `DuplicateNameError`, `IncompatibleAxisError` |
| `rename`, `rename_with`, `relocate`, `mutate`, `add_count`, `add_tally` | Independent by default; `inplace=True` returns exact input identity | Metadata-only writes; matrix sources are read-only | `SelectionError`, `DuplicateNameError`, `SizeMismatchError`, `AnnplyrError` |
| `transmute` | Always independent; no ownership argument | Replaces metadata columns while reading projected sources | selection/source/size/budget errors |
| grouped AnnData verbs | New grouped wrapper, or identical wrapper for explicit in-place calls; `ungroup()` returns AnnData | One shared positional group plan; keys retained or updated | selection, axis, join, size, and budget errors |
| summaries, count/tally, and `group_*` inspection | pandas DataFrame/list-like metadata; no AnnData mutation | Canonical matrix reducers use deterministic bounded chunks; metadata and opaque expressions read their planned sources | selection/source/axis/budget errors |
| `pull`, `to_wide`, `to_tidy`, `pivot_longer`, `as_frame`, `nest_by` | pandas/Series or nested extraction result; input unchanged | Explicit pandas materialization with cumulative budgets | source, selection, name-repair, and budget errors |
| expressions, selectors, and aggregators | `AnnplyrExpr`, selector, or expansion metadata; input unchanged | Lazy until evaluated; raw Narwhals is opaque and fully charged | `SelectionError`, `UnknownColumnError`, `UnknownSourceError` |
| dataframe rectangling helpers | New pandas objects; input unchanged | Operates after explicit tabular extraction | `SelectionError`, `DuplicateNameError`, `SizeMismatchError` |
| single-cell utilities | Readers return tables/reports; writers are independent unless `inplace=True` | Metadata/name/palette operations only | `SelectionError`, `DuplicateNameError`, `SizeMismatchError`, `JoinRelationshipError` |

All cumulative `max_matrix_values=` checks plan every source before the first
adapter read. The mechanically checked callable-level contract is recorded in
{doc}`development/api-contract-v0.3`.

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
from annplyr.typing import AnnDataWithAnnplyr


def select_markers(adata: AnnData) -> AnnDataWithAnnplyr:
    typed = cast(AnnDataWithAnnplyr, adata)
    return typed.ap.filter(ap.col("qc_pass")).select(x=["MS4A1", "CD79A"])
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
| `Source` | the exact `as_frame()` source literals |
| `SourceSelectors` | keyed selectors for aligned mapping sources |
| `JoinInput`, `JoinBy`, and join-option aliases | the exact dataframe/mapping and validated literal choices accepted by joins |
| `GroupedReturn` | typed AnnData or `GroupedAnnData` |

The generated {doc}`development/public-typing-contract` inventories every root
symbol, public method signature, overload, return annotation, and alias.

## Flat verb functions

`ap.filter(adata, ...)` is the same call as `adata.ap.filter(...)`.

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.add_count
   annplyr.add_tally
   annplyr.anti_join
   annplyr.arrange
   annplyr.as_frame
   annplyr.count
   annplyr.distinct
   annplyr.filter
   annplyr.filter_exp
   annplyr.filter_obs
   annplyr.full_join
   annplyr.group_by
   annplyr.inner_join
   annplyr.left_join
   annplyr.mutate
   annplyr.nest_by
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
   annplyr.tally
   annplyr.to_tidy
   annplyr.to_wide
   annplyr.transmute
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

## Dataframe Extraction Helpers

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.chop
   annplyr.drop_na
   annplyr.extract
   annplyr.fill
   annplyr.hoist
   annplyr.nest
   annplyr.pack
   annplyr.pivot_wider
   annplyr.separate
   annplyr.separate_rows
   annplyr.unchop
   annplyr.unite
   annplyr.unpack
   annplyr.unnest
   annplyr.unnest_longer
   annplyr.unnest_wider
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
   annplyr.add_sample_meta
   annplyr.feature_present
   annplyr.get_palette
   annplyr.name_duplicates
   annplyr.rename_obs_names
   annplyr.rename_var_names
   annplyr.replace_name_suffix
   annplyr.sample_meta
   annplyr.sample_summary
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
