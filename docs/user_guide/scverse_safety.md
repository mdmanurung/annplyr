# Scverse safety

`annplyr` works with AnnData storage relationships rather than replacing them.
The practical rules are short:

1. keep transformations as AnnData until a table is actually needed;
2. select matrix features explicitly;
3. use a finite budget for potentially large projections;
4. call `.to_memory()` before mutating backed objects;
5. use `inplace=True` only when shared mutation is intentional.

## Ownership and alignment

Axis-changing verbs return independent objects by default and subset every
aligned container by integer position. PBMC3K carries `X`, `raw`, four
embeddings, PCA loadings, and two pairwise graphs. One filter keeps all of them
consistent:

```{testcode}
subset = adata.ap.filter(obs=ap.col("louvain") == "B cells")

print(subset.shape, subset.raw.shape)
print(subset.obsm["X_umap"].shape, subset.obsp["connectivities"].shape)
print(adata.n_obs)
```

```{testoutput}
(342, 1838) (342, 13714)
(342, 2) (342, 342)
2638
```

The pairwise graph was subset on **both** axes, `raw` kept all 13,714 genes for
the 342 retained cells, and the original object still has its 2,638 cells.
`copy=False` remains non-mutating but may produce a view or a materialized
selection, so do not branch on `is_view` after it.

Same-shape mutations use `inplace=True` only when exact input identity is
required. Validation completes before the first write.

## Sparse matrices stay sparse

Selected sparse sources stay sparse during projection. Wide pandas exports use
pandas sparse columns, and AnnData-returning projections retain the input
CSR/CSC format:

```{testcode}
print(adata.ap.to_df(raw=["MS4A1", "CD79A"], max_matrix_values=2 * adata.n_obs).dtypes)
```

```{testoutput}
raw_MS4A1    Sparse[float32, 0.0]
raw_CD79A    Sparse[float32, 0.0]
dtype: object
```

Expressions that combine two matrix columns densify only those projected
columns for the duration of the evaluation, because pandas cannot do sparse
arithmetic in every subtype. The result matches the dense path exactly, and the
exported frames above keep their subtype.

## Backed AnnData

Read-only projection, filtering, ordering, and export are supported for backed
dense, CSR, and CSC sources. `copy=True` axis operations return an independent
in-memory AnnData. Same-shape mutation is rejected while the object is backed;
load it explicitly when mutation is intended:

```python
memory = backed.to_memory()
memory.ap.mutate(obs={"new": ap.col("old")}, inplace=True)
```

This keeps the move from disk-backed analysis to in-memory mutation visible in
the workflow rather than implicit.

## Raw and pairwise sources

`raw` is an explicit, read-only expression and extraction source. Pairwise
`obsp`/`varp` matrices and tabular `uns` entries are extraction-only through
`pull()` or `as_frame()`; they are not mutation targets or join keys.

```{testcode}
print(adata.ap.to_df(raw=["MS4A1"], max_matrix_values=adata.n_obs).shape)
```

```{testoutput}
(2638, 1)
```

## Cumulative materialization budgets

Use `max_matrix_values=` on filtering, ordering, distinct selection, mutation,
summaries, extraction, and their grouped counterparts. The planner charges
projected cells cumulatively across sources and rejects the complete request
before its first adapter read:

```{testcode}
try:
    adata.ap.to_tidy(raw=["MS4A1", "CD79A"], max_matrix_values=100)
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (to_tidy) would materialize 5276 matrix values, which exceeds max_matrix_values=100
```

Canonical matrix summaries then read deterministic row or feature chunks rather
than the complete selected source. Sparse chunks stay sparse; there is no
whole-source dense fallback. Exact `median` and `n_distinct` summaries may
retain state proportional to one reduction vector, but feature batches are
finalized before the next batch is read. Opaque raw Narwhals or cross-row
custom summary expressions keep their eager compatibility path, so give those a
finite `max_matrix_values=` budget.

Whole-matrix long exports additionally require `allow_all_features=True`. Set
both the opt-in and a finite budget when full materialization is deliberate.
