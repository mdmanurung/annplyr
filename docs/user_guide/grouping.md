# Persistent grouping

Use `group_by()` when several steps should use the same sample, condition, cell
type, or feature groups. It resolves existing metadata keys for one AnnData
axis and returns a `GroupedAnnData` wrapper.

For a single summary, `summarize(..., by=...)` is shorter. Persistent grouping
becomes useful for a sequence such as rank within sample, retain the top cells,
then compute sample-level counts.

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
follow first-seen order, include NA values, preserve categorical dtype and
order, and omit unobserved categories such as the fixture's `unused` level.

## Inspect groups before transforming

Iteration yields each key and an aligned AnnData subset:

```{testcode}
observed = [(key, group.n_obs) for key, group in grouped]
assert observed == [({"batch": "A"}, 2), ({"batch": "B"}, 2)]
```

Iteration is useful for handing each aligned subset to project code. Most
workflows can remain vectorized through grouped verbs: AnnData-returning methods
return another grouped wrapper, while summary and count methods return pandas
tables.

```{testcode}
ranked = grouped.mutate(obs={"within_batch": ap.row_number()})
top = ranked.slice_max(ap.col("n_counts"), n=1)

assert top.group_vars() == ["batch"]
assert top.ungroup().n_obs == 2
assert ranked.summarize(obs={"cells": ap.n()})["cells"].tolist() == [2, 2]
```

Call `ungroup()` at the point where an ordinary AnnData object is needed by
Scanpy, serialization, or another accessor pipeline.

## A sample-level ranking pattern

The pattern below keeps the highest-count cell from each donor while retaining
grouping until the selection is complete:

```{testcode}
top_per_donor = (
    adata.ap.group_by(obs="donor_id")
    .mutate(obs={"within_donor_rank": ap.min_rank("n_counts", descending=True)})
    .filter(obs=ap.col("within_donor_rank") == 1)
    .ungroup()
)

assert list(top_per_donor.obs_names) == ["cell_2", "cell_3"]
```

## Keep group keys valid

Grouped `select`, `relocate`, and `transmute` retain grouping keys
automatically. `rename` and `rename_with` update the group specification.

```{testcode}
selected = grouped.select(obs=["cell_type"])
assert list(selected.ungroup().obs.columns) == ["batch", "cell_type"]

renamed = grouped.rename(obs={"sample": "batch"})
assert renamed.group_vars() == ["sample"]
```

If grouped `mutate()` changes a key, expressions in that call use the groups
established at its start. The returned wrapper resolves subsequent verbs from
the final key values. All six grouped joins run once globally and update
suffixed key names when necessary.

Grouping is axis-specific. Pass `var=` for feature groups such as chromosome or
feature type. Computed, virtual-name, cross-axis, and empty grouping
specifications raise typed selection or axis errors.
