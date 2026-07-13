# Grouping

Use `group_by()` when operations should happen within observation or feature
groups.

```python
grouped = adata.ap.group_by(obs="cell_type")
```

Grouped objects support iteration:

```python
for key, group in grouped:
    print(key, group.n_obs)
```

They also support group-local verbs:

```python
grouped.summarize(obs={"cells": ap.n(), "mean_counts": ap.mean("n_counts")})
grouped.mutate(obs={"row": ap.row_number()})
grouped.slice_max(ap.col("n_counts"), n=3)
```

Pass `var=` to group the feature axis:

```python
adata.ap.group_by(var="feature_type").summarize(var={"features": ap.n()})
```

Grouping is intentionally axis-aware. A grouped operation cannot silently mix
observation and feature grouping.
