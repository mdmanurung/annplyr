# Migrating from v0.2 to v0.3

Version 0.3 makes ownership, grouping, positional identity, expression
inspection, and matrix-read budgets explicit. Review this guide before removing
version pins from code written for v0.2.

## Copy and in-place arguments

Verbs that can change an AnnData axis or its order now default to an independent
result with `copy=True`. This includes `filter`, `select`, `arrange`, `distinct`,
`slice` and every `slice_*` variant, plus all six joins. Existing code that
relied on a view must request `copy=False` explicitly:

```python
independent = adata.ap.filter(obs=ap.col("batch") == "A")
possibly_deferred = adata.ap.filter(obs=ap.col("batch") == "A", copy=False)
```

`copy=False` permits either an AnnData view or a materialized object. Do not
test `is_view` as part of application logic. Test values and ordering, and use
`copy=True` when the result must not share mutable state with its input.

Same-shape operations now use `inplace=False` instead of `copy=`. The affected
accessor families are metadata mutation, metadata renaming/reordering, and
add-count/tally operations. The same rule applies to `add_sample_meta`,
`rename_obs_names`, `rename_var_names`, `add_name_prefix`,
`replace_name_suffix`, and `store_palette`:

```python
# v0.2
out = adata.ap.mutate(obs={"qc": ap.col("n_counts") > 1000}, copy=True)

# v0.3
out = adata.ap.mutate(obs={"qc": ap.col("n_counts") > 1000})
adata.ap.mutate(obs={"qc": ap.col("n_counts") > 1000}, inplace=True)
```

An in-place operation validates before writing and returns the identical input
object. `relocate(..., inplace=True)` changes metadata column order on that
object. Sorted `add_count` and `add_tally` may reorder the relevant axis while
preserving every aligned AnnData container. `transmute()` has no `copy` or
`inplace` argument and always returns an independent object.

On backed AnnData, `copy=True` axis operations produce an independent in-memory
selection. Same-shape mutation remains unsupported until the object is loaded
with `to_memory()`. A `copy=False` axis operation remains non-mutating but may
return a view or a materialized selection.

## Duplicate axis labels are positional

Axis labels no longer determine row reconstruction. Filtering, arranging,
distinct selection, slicing, and joins track integer positions, so duplicate
`obs_names` and `var_names` preserve the selected records and the alignment of
`X`, layers, `raw`, `obsm`, `varm`, `obsp`, and `varp`.

Code that used duplicate labels as if they were unique keys must still provide
an explicit unique metadata join key. The positional change protects AnnData
identity; it does not redefine join cardinality.

## Grouping persists

`group_by()` returns `GroupedAnnData`. Every AnnData-returning grouped verb now
returns another grouped wrapper, so pipelines remain grouped until
`ungroup()`:

```python
result = (
    adata.ap.group_by(obs="batch")
    .filter(obs=ap.col("n_counts") > 1000)
    .mutate(obs={"within_batch": ap.row_number()})
    .select(obs=["batch", "within_batch"])
    .ungroup()
)
```

Use `ungroup()` when code needs the wrapped AnnData. Iteration still yields
`(key, group)` pairs. Grouped non-in-place calls return a new wrapper; grouped
in-place calls return the identical wrapper.

Grouping keys are resolved once from existing metadata columns. Computed,
virtual, and empty group specifications are rejected. Groups follow first-seen
row order, include NA keys, retain categorical dtype/order, observe only groups
that occur, and drop unobserved categorical levels. `group_data()[".rows"]`
now contains zero-based integer positions rather than axis labels.

Grouped `select`, `relocate`, and `transmute` retain grouping columns
automatically and place them deterministically. `rename` and `rename_with`
update the group specification. When grouped `mutate` changes a key, that call
uses the old grouping plan and later calls regroup from the final values.

All six joins are available on a grouped wrapper. They execute globally rather
than once per group, update suffixed key names when needed, and preserve
grouping only when the result still contains valid keys. Right/full joins still
cannot add right-only cells or features; the default `unmatched="error"` makes
that divergence from unrestricted dataframe joins explicit.

## Expressions wrap Narwhals expressions

Expression helpers now return an `AnnplyrExpr`. The wrapper carries dependency,
output-width, and cardinality information used for safe projection and
batching. Operators, methods, namespaces, aliases, and virtual axis-name
expressions preserve the wrapper:

```python
expr = ap.col("n_counts") + 1
raw_expr = expr.to_narwhals()
```

Pass annplyr expressions to accessor methods as before. Use `to_narwhals()`
only for explicit interoperation with an API that requires a raw
`narwhals.Expr`. A raw Narwhals expression passed into annplyr is treated as
opaque: it remains supported, but budget planning charges the whole relevant
source because dependencies cannot be proven from public Narwhals metadata.

`where()` is now a schema/dtype selector. Its predicate receives a zero-length
typed pandas Series, so predicates may inspect dtype or schema but must not
depend on data values:

```python
import pandas as pd

numeric = ap.where(lambda series: pd.api.types.is_numeric_dtype(series.dtype))
```

Move value-dependent predicates into `filter`, `if_any`, or `if_all`.

## Matrix projection and read budgets

Version 0.3 projects selected rows and columns before conversion for dense,
CSR/CSC, sparse-array, DataFrame, and supported backed sources. Wide sparse
exports preserve pandas sparse columns, and AnnData outputs preserve the
original sparse format. Sorted reads required by a backed store are restored to
the requested order before return.

`max_matrix_values=` is now available on `filter`, `arrange`, `distinct`,
`mutate`, `transmute`, `summarize`/`summarise`, `pull`, `to_df`, `to_tidy`,
`pivot_longer`, `as_frame`, and grouped counterparts. The planner validates all
requests before the first adapter read. It charges projected matrix cells
cumulatively across sources; opaque expressions conservatively charge their
full source. Negative budgets and requests above the limit fail before reading:

```python
table = adata.ap.to_df(
    obs=["batch"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=2 * adata.n_obs,
)
```

`to_tidy` and matrix-backed `pivot_longer` still require explicit features
unless `allow_all_features=True` is supplied.

## Joins and typed failures

Joins never add or duplicate AnnData axis records. Left joins enrich metadata;
inner, semi, and anti joins may subset by position. Right and full joins accept
only results representable on the existing AnnData axis and raise
`JoinRelationshipError` for right-only records. Relationship, `multiple`,
`unmatched`, `na_matches`, and suffix rules are validated before construction.

Public contract violations use package errors:

- `SelectionError`, `UnknownColumnError`, or `UnknownSourceError` for invalid
  selection/expression requests;
- `DuplicateNameError` and `NameRepairError` for ambiguous output schemas;
- `IncompatibleAxisError` and `SizeMismatchError` for axis/length mismatches;
- `JoinRelationshipError` for unsafe join cardinality or right-only records;
- `AnnplyrError` for invalid ownership, backed mutation, and matrix-budget
  operations not covered by a more specific subclass.

The complete mechanically checked contract is recorded in
{doc}`development/api-contract-v0.3`.

## Typing before 1.0

The pre-1.0 typing surface now exports `AnnplyrAccessor`, `GroupedAnnData`, and
`AnnplyrExpr`. Code that imported those names from private modules should move
to the package root:

```python
# old, unsupported private import
from annplyr._grouped import GroupedAnnData

# supported
from annplyr import AnnplyrAccessor, AnnplyrExpr, GroupedAnnData
```

Static checkers do not see dynamically registered attributes on the upstream
`AnnData` class. Use one explicit typing-only cast at the boundary where an
AnnData enters an annplyr pipeline:

```python
from typing import cast

from annplyr.typing import AnnDataWithAnnplyr

typed = cast(AnnDataWithAnnplyr, adata)
result = typed.ap.group_by(obs="batch").filter(obs=ap.col("qc_pass")).ungroup()
```

The cast is a no-op: `AnnDataWithAnnplyr` is `AnnData` at runtime. The public
selector, expression, source, join-input, and grouped-return aliases are in
`annplyr.typing`. No runtime call forms were added or widened.
