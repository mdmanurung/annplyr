# Scverse Safety

`annplyr` preserves AnnData storage relationships; it does not replace AnnData
storage semantics.

## Ownership and alignment

Axis-changing verbs return independent objects by default and subset every
aligned container by integer position. This keeps `X`, layers, `raw`, `obsm`,
`varm`, `obsp`, and `varp` consistent even when axis names are duplicated.
`copy=False` remains non-mutating but may produce a view or a materialized
selection.

Same-shape mutations use `inplace=True` only when exact input identity is
required. Validation completes before the first write.

## Sparse matrices

Selected sparse sources stay sparse during projection. Wide pandas exports use
pandas sparse columns where possible, and AnnData-returning projections retain
the input CSR/CSC format.

```{testcode}
selected = adata.ap.to_df(x=["MS4A1", "CD79A"], max_matrix_values=8)
assert selected.shape == (4, 2)
```

## Backed AnnData

Read-only projection, filtering, ordering, and export are supported for backed
dense, CSR, and CSC sources. `copy=True` axis operations return independent
in-memory AnnData. Same-shape mutation is rejected while the object is backed;
load it explicitly when mutation is intended:

```python
memory = backed.to_memory()
memory.ap.mutate(obs={"new": ap.col("old")}, inplace=True)
```

Backed benchmarks are warm-cache measurements unless the caller controls the
filesystem cache.

## Raw and pairwise sources

`raw` is an explicit, read-only expression or extraction source. Pairwise
`obsp`/`varp` matrices and tabular `uns` entries are extraction-only through
`pull()` or `as_frame()`; they are not mutation targets or join keys.

```{testcode}
raw_table = adata.ap.to_df(raw=["MS4A1"], max_matrix_values=4)
assert raw_table.shape == (4, 1)
```

## Cumulative materialization budgets

Use `max_matrix_values=` on filtering, ordering, distinct selection, mutation,
summaries, extraction, and grouped counterparts. The planner charges projected
cells cumulatively across sources and rejects the complete request before its
first adapter read.

```{testcode}
bounded = adata.ap.to_tidy(
    x=["MS4A1", "CD79A"],
    max_matrix_values=8,
)
assert len(bounded) == 8
```

Whole-matrix long exports still require `allow_all_features=True`. Set both the
opt-in and a finite budget when full materialization is deliberate.
