# Core Concepts

`annplyr` treats AnnData as an aligned object, not just a collection of pandas
tables. Verbs that return AnnData must preserve alignment across coordinated
containers.

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

## Sources

Expressions can read from:

- `obs` metadata
- `var` metadata
- selected `X` or layer columns through `x=`
- `raw` through `raw=`
- `obsm` and `varm` matrices through keyed mappings

Controlled extraction also supports `obsp`, `varp`, and tabular `uns` values
through `as_frame()`.

## Alignment

AnnData-returning verbs use AnnData-native slicing so `X`, layers, `obsm`,
`varm`, `obsp`, and `varp` remain aligned after subsetting or reordering.

## Errors

`annplyr` raises typed package errors for invalid selectors, missing sources,
unsafe joins, duplicate names, and incompatible axis operations. See
{doc}`../api` for the full error reference.

## Design lineage

`annplyr` draws direct inspiration from
[annsel](https://github.com/srivarra/annsel), which introduced
predicate-based selection on AnnData objects. `annplyr` extends that idea to
the full `dplyr`/`tidyr` verb set — `mutate`, `summarize`, `group_by`, joins,
and tidy extraction — for R tidyverse users moving to Python single-cell
analysis.
