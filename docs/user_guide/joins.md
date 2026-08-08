# AnnData-Safe Joins

Joins enrich or subset `obs` or `var` metadata without manufacturing aligned
matrix rows. They track the AnnData axis by integer position, so duplicate axis
labels remain distinct.

## Enrich metadata

```{testcode}
sample_table = pd.DataFrame(
    {"batch": ["A", "B"], "condition": ["treated", "control"]}
)
joined = adata.ap.left_join(sample_table, by="batch", axis="obs")

assert joined.n_obs == adata.n_obs
assert joined.obs["condition"].tolist() == ["treated", "control", "treated", "control"]
assert "condition" not in adata.obs
```

Joins return independent AnnData by default. Pass `copy=False` only when either
a view or a materialized result is acceptable.

## Filter an axis

Inner, semi, and anti joins may subset an axis while preserving original
left-side order.

```{testcode}
allowed = pd.DataFrame({"batch": ["A"]})
subset = adata.ap.semi_join(allowed, by="batch", axis="obs")
removed = adata.ap.anti_join(allowed, by="batch", axis="obs")

assert list(subset.obs_names) == ["cell_0", "cell_2"]
assert list(removed.obs_names) == ["cell_1", "cell_3"]
```

## Cardinality and right-only records

Use `relationship`, `multiple`, `unmatched`, `na_matches`, and `suffixes` to
make join assumptions explicit. The default relationship is `many-to-one`.

```{testcode}
checked = adata.ap.left_join(
    sample_table,
    by="batch",
    axis="obs",
    relationship="many-to-one",
)
assert checked.n_obs == 4
```

A right or full join cannot add a right-only cell or feature because no aligned
`X`, layer, embedding, or pairwise record exists for it. Such requests raise
`JoinRelationshipError` rather than returning a partly fabricated AnnData.
Duplicate right keys, unmatched policies, invalid suffixes, and cross-axis
requests also fail with typed package errors.

Grouped joins execute globally and return a grouped wrapper when the resulting
metadata still contains valid grouping keys. Suffixing a key updates the stored
group specification deterministically.
