# Persistent Grouping

`group_by()` resolves existing metadata keys for one AnnData axis and returns a
`GroupedAnnData` wrapper.

```{testcode}
grouped = adata.ap.group_by(obs="batch")

assert grouped.group_vars() == ["batch"]
assert grouped.group_keys().to_dict("list") == {"batch": ["A", "B"]}
assert grouped.group_data().to_dict("list") == {
    "batch": ["A", "B"],
    ".rows": [[0, 2], [1, 3]],
}
```

The `.rows` entries are zero-based integer positions, not axis labels. Groups
follow first-seen order, include NA values, preserve categorical dtype/order,
and omit unobserved categories such as the fixture's `unused` level.

## Iteration and local verbs

Iteration yields each key and an aligned AnnData subset:

```{testcode}
observed = [(key, group.n_obs) for key, group in grouped]
assert observed == [({"batch": "A"}, 2), ({"batch": "B"}, 2)]
```

AnnData-returning grouped verbs return another grouped wrapper. Summary/count
methods return pandas tables.

```{testcode}
ranked = grouped.mutate(obs={"within_batch": ap.row_number()})
top = ranked.slice_max(ap.col("n_counts"), n=1)

assert top.group_vars() == ["batch"]
assert top.ungroup().n_obs == 2
assert ranked.summarize(obs={"cells": ap.n()})["cells"].tolist() == [2, 2]
```

Call `ungroup()` to recover an ordinary AnnData object.

## Key retention and updates

Grouped `select`, `relocate`, and `transmute` retain grouping keys
automatically. `rename` and `rename_with` update the group specification.

```{testcode}
selected = grouped.select(obs=["cell_type"])
assert list(selected.ungroup().obs.columns) == ["batch", "cell_type"]

renamed = grouped.rename(obs={"sample": "batch"})
assert renamed.group_vars() == ["sample"]
```

When grouped `mutate` changes a key, that call uses the old plan; the returned
wrapper resolves future calls from the final key values. All six grouped joins
run once globally and update suffixed key names when necessary.

Grouping is axis-specific. Pass `var=` for feature groups. Computed,
virtual-name, cross-axis, and empty grouping specifications raise typed
selection or axis errors.
