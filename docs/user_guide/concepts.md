# Core Concepts

`annplyr` treats AnnData as an aligned object rather than a loose collection of
pandas tables. Verbs that return AnnData must preserve alignment across
coordinated containers.

## Axes

`obs` is the observation axis, usually cells. `var` is the feature axis, usually
genes or measurements.

Most row-like operations default to the `obs` axis:

```python
adata.ap.slice_head(n=5)
```

Pass `axis="var"` for feature-axis operations:

```python
adata.ap.slice_head(n=10, axis="var")
```

## Ownership

Axis-changing and ordering verbs default to `copy=True`, returning an
independent AnnData object. `copy=False` is a performance request, not a promise
of `is_view=True`: implementations may return a view or a materialized result.
Same-shape metadata operations use `inplace=False`; an explicit in-place call
validates first and returns the identical input object. `transmute()` is always
independent.

All axis operations track integer positions. Duplicate `obs_names` or
`var_names` therefore remain distinct and every aligned AnnData container keeps
the same positional identity.

## Sources and planning

Expressions can read from:

- `obs` metadata
- `var` metadata
- selected `X` or layer columns through `x=`
- `raw` through `raw=`
- `obsm` and `varm` matrices through keyed mappings

Controlled extraction also supports `obsp`, `varp`, and tabular `uns` values
through `as_frame()`.

The v0.3 planner resolves every requested source and its projected rows/columns
before the first adapter read. `max_matrix_values=` is cumulative across
matrix sources. A raw Narwhals expression has opaque dependencies and is
charged against its full source; an annplyr expression exposes conservative
dependency metadata for narrower projection.

## Expressions

Helpers such as `col`, `lit`, `if_else`, selectors, and virtual axis names
return `AnnplyrExpr` wrappers. Operators and method chains preserve the wrapper.
Call `to_narwhals()` only when another API explicitly requires a raw
`narwhals.Expr`. The `where()` selector receives a zero-length typed Series and
must inspect schema/dtype rather than data values.

## Alignment

AnnData-returning verbs use AnnData-native slicing so `X`, layers, `obsm`,
`varm`, `obsp`, and `varp` remain aligned after subsetting or reordering.

## Typed failures

`annplyr` raises typed package errors for invalid selectors, missing sources,
unsafe joins, duplicate names, and incompatible axis operations. See
{doc}`../api` for the full error reference.

- `SelectionError`, `UnknownColumnError`, and `UnknownSourceError` identify
  invalid selection or expression requests.
- `DuplicateNameError` and `NameRepairError` identify ambiguous output schemas.
- `IncompatibleAxisError` and `SizeMismatchError` identify invalid axis or
  length relationships.
- `JoinRelationshipError` identifies a join that would duplicate or add AnnData
  axis records.
- `AnnplyrError` covers ownership, backed mutation, and budget violations that
  do not have a narrower subclass.

## Design lineage

`annplyr` draws direct inspiration from
[annsel](https://github.com/srivarra/annsel), which introduced
predicate-based selection on AnnData objects. `annplyr` extends that idea to
the full `dplyr`/`tidyr` verb set (`mutate`, `summarize`, `group_by`, joins,
and tidy extraction) for R tidyverse users moving to Python single-cell
analysis.
