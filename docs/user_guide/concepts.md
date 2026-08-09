# Core concepts

`annplyr` treats AnnData as one aligned data structure. `obs` and `var` are
convenient metadata tables, but their rows are coordinated with `X`, layers,
`raw`, embeddings, loadings, and pairwise matrices. Every verb is designed
around that relationship.

## Read a pipeline by its return type

Most annplyr workflows alternate between two kinds of result:

| Operation | Typical result | Use it for |
|---|---|---|
| filter, select, mutate, arrange, slice, join | AnnData | Continue analysis or call `.ap` again |
| grouped AnnData verbs | `GroupedAnnData` | Reuse group-local semantics across steps |
| summarize, count, extraction | pandas object | Plot, report, model, or export |

That distinction makes materialization visible. A cell filter remains AnnData;
a long expression table is deliberately pandas.

## Axes are explicit

`obs` is the observation axis, usually cells. `var` is the feature axis,
usually genes, proteins, peaks, or other measurements. Row-like operations
default to cells:

```python
first_cells = adata.ap.slice_head(n=100)
```

Use `axis="var"` when the same verb should operate on features:

```python
top_features = adata.ap.slice_max(ap.col("highly_variable_rank"), n=2_000, axis="var")
```

Arguments also name their source. This call filters cells by `obs`, retains
feature annotations from `var`, and projects genes from `X`:

```python
subset = adata.ap.filter(obs=ap.col("qc_pass")).ap.select(
    var=["gene_symbol", "highly_variable"],
    x=["MS4A1", "CD79A", "CD3D"],
)
```

## Sources stay attached to AnnData

Expressions can read from:

- `obs` and `var` metadata;
- selected features in `X` or a named layer through `x=`;
- the immutable raw snapshot through `raw=`;
- keyed `obsm` and `varm` matrices.

Controlled extraction through `as_frame()` also covers `obsp`, `varp`, and
tabular `uns` values. The source is always explicit, which prevents a gene name
from being confused with an `obs` column of the same name.

```python
annotated = adata.ap.mutate(
    obs={"high_quality": ap.col("pct_counts_mt") < 10},
    x={
        "B_markers_detected": (ap.col("MS4A1") > 0)
        & (ap.col("CD79A") > 0),
    },
    obsm={"X_pca": {"PC1": ap.col("0")}},
)
```

The matrix-like inputs above are read-only. `mutate()` writes the three derived
values to metadata and leaves each source matrix unchanged.

## Ownership is predictable

Axis-changing and ordering verbs return independent AnnData objects by
default. `copy=False` asks annplyr to avoid a guaranteed copy, but the result
may be either an AnnData view or a materialized object depending on what is
safe. Do not branch logic on `is_view` after `copy=False`.

Same-shape metadata operations use `inplace=True` when exact input identity is
required:

```python
adata.ap.mutate(obs={"qc_pass": ap.col("pct_counts_mt") < 10}, inplace=True)
```

Validation completes before the first write. Without `inplace=True`, the
original AnnData remains unchanged.

## Positional identity protects duplicate names

Real objects can contain repeated `obs_names` after concatenation or repeated
feature labels in alternate identifier columns. annplyr tracks integer
positions for filtering, ordering, joins, and grouped operations, so duplicate
labels do not collapse distinct cells or features. All aligned containers use
the same positional subset.

## Expressions and selectors compose

`col()`, `lit()`, rank helpers, null helpers, and selectors return
`AnnplyrExpr` objects. Operators preserve those objects, so a domain rule can
be built once and reused:

```python
qc_rule = (ap.col("total_counts") >= 1_000) & (ap.col("pct_counts_mt") < 10)
marker_genes = ap.any_of(["MS4A1", "CD79A", "CD3D"])

analysis = adata.ap.filter(obs=qc_rule).ap.select(x=marker_genes)
```

`where()` is a schema selector: its predicate receives a zero-length typed
Series and should inspect dtype, not values. Call `to_narwhals()` only when an
external API explicitly needs a raw Narwhals expression.

## Grouping can persist

Use `by=` for one summary and `group_by()` when several operations share the
same groups:

```python
top_cells = (
    adata.ap.group_by(obs="sample_id")
    .mutate(obs={"within_sample_rank": ap.min_rank("total_counts", descending=True)})
    .filter(obs=ap.col("within_sample_rank") <= 10)
    .ungroup()
)
```

The wrapper records group keys, not a detached copy of the data. See
{doc}`grouping` for group order, missing values, categorical keys, and key
updates.

## Projection and budgets make reads deliberate

annplyr resolves selected rows and columns before reading matrix sources. Pass
`max_matrix_values=` when a workflow needs a hard cumulative bound:

```python
plot_data = adata.ap.to_tidy(
    obs=["sample_id", "cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=2 * adata.n_obs,
)
```

The complete request is rejected before the first matrix read if it exceeds the
budget. Canonical scalar summaries use deterministic internal chunks for dense,
CSR, CSC, and backed sources; users specify the scientific projection, not an
implementation-specific chunk size.

## Typed errors expose invalid assumptions

Invalid selectors, missing sources, incompatible axes, unsafe joins, duplicate
outputs, and over-budget reads raise annplyr-specific exceptions. Catch a
specific error when the workflow can recover; otherwise let the message expose
the failed assumption.

The complete error list is in {doc}`../api`.

## Design lineage

`annplyr` is inspired by
[annsel](https://github.com/srivarra/annsel), which introduced
predicate-based selection for AnnData. annplyr extends that idea with
tidyverse-style mutation, summaries, grouping, joins, and rectangling for users
who want a familiar dataframe grammar inside Python single-cell workflows.
