# API

## Accessor

The `adata.ap` accessor is registered when `annplyr` is imported. AnnData-returning
verbs preserve axis alignment for `X`, layers, `obsm`, `varm`, `obsp`, and
`varp` through AnnData-native slicing.

Joins are metadata joins for `obs` or `var`. They may enrich or subset an axis,
but they raise `JoinRelationshipError` when a requested join would add or
duplicate cells/features. Long matrix exports materialize data into pandas and
therefore require explicit feature selection unless `allow_all_features=True`.
Mutating verbs raise `AnnplyrError` on backed AnnData objects; call
`.to_memory()` first when mutation is intentional.

```{eval-rst}
.. autoclass:: annplyr._accessor.AnnplyrAccessor
   :members:

.. autoclass:: annplyr._grouped.GroupedAnnData
   :members:
```

## Expressions And Selectors

```{eval-rst}
.. autosummary::
   :toctree: generated

   annplyr.all_of
   annplyr.any_of
   annplyr.between
   annplyr.case_when
   annplyr.col
   annplyr.coalesce
   annplyr.contains
   annplyr.cum_max
   annplyr.cum_min
   annplyr.cum_prod
   annplyr.cum_sum
   annplyr.cume_dist
   annplyr.dense_rank
   annplyr.desc
   annplyr.ends_with
   annplyr.everything
   annplyr.first
   annplyr.if_else
   annplyr.is_na
   annplyr.lag
   annplyr.last
   annplyr.last_col
   annplyr.lead
   annplyr.lit
   annplyr.matches
   annplyr.min_rank
   annplyr.n_distinct
   annplyr.na_if
   annplyr.nth
   annplyr.num_range
   annplyr.obs_names
   annplyr.percent_rank
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

   annplyr.pivot_wider
   annplyr.unnest
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
