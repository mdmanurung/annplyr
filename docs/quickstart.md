# Quickstart

`annplyr` exposes dataframe-style verbs through `adata.ap`. The examples on
this page run against a four-cell AnnData fixture during the documentation
doctest build.

## Filter observations

Importing annplyr registers the accessor. Axis-changing verbs return an
independent AnnData object by default.

```{testcode}
import annplyr as ap

filtered = adata.ap.filter(
    obs=ap.col("batch") == "A",
    x=ap.col("MS4A1") > 0,
)

assert filtered.n_obs == 2
assert filtered.n_vars == 3
assert not filtered.is_view
```

Pass `copy=False` only when either a view or a materialized result is acceptable.

## Select metadata and features

```{testcode}
selected = adata.ap.select(
    obs=["batch", "cell_type"],
    x=["MS4A1", "CD79A"],
)

assert selected.shape == (4, 2)
assert list(selected.obs.columns) == ["batch", "cell_type"]
assert list(selected.var_names) == ["MS4A1", "CD79A"]
```

## Add metadata from matrix sources

`mutate()` writes only metadata columns. Matrix arguments provide read-only
expression sources. Same-shape operations return an independent object unless
`inplace=True` is explicit.

```{testcode}
annotated = adata.ap.mutate(
    obs={"high_counts": ap.col("n_counts") > 1_000},
    x={"MS4A1_value": ap.col("MS4A1")},
)

assert "high_counts" in annotated.obs
assert "MS4A1_value" in annotated.obs
assert "high_counts" not in adata.obs
```

## Keep grouping through a pipeline

AnnData-returning grouped verbs preserve grouping. Call `ungroup()` at the
boundary where an ordinary AnnData object is required.

```{testcode}
grouped = (
    adata.ap.group_by(obs="batch")
    .filter(obs=ap.col("n_counts") >= 1_000)
    .mutate(obs={"within_batch": ap.row_number()})
)

assert grouped.group_vars() == ["batch"]
assert grouped.ungroup().n_obs == 3

summary = grouped.summarize(obs={"cells": ap.n()})
assert summary.to_dict("list") == {"batch": ["A", "B"], "cells": [2, 1]}
```

Groups use first-seen order, include NA keys, observe only categorical values
that occur, and expose zero-based integer positions through
`group_data()[".rows"]`.

## Export bounded tables

Matrix-backed exports project selected features before conversion. The budget
is cumulative across matrix sources and is checked before the first read.

```{testcode}
wide = adata.ap.to_df(
    obs=["cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)
long = adata.ap.to_tidy(
    obs=["cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)

assert wide.shape == (4, 3)
assert long.shape == (8, 4)
assert list(long.columns) == ["obs_name", "feature", "value", "cell_type"]
```

Continue with {doc}`user_guide/concepts`, review the breaking changes in
{doc}`migration-v0.3`, or look up exact signatures in {doc}`api`.
