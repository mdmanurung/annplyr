# Core concepts

`annplyr` treats AnnData as one aligned data structure. `obs` and `var` are
convenient metadata tables, but their rows are coordinated with `X`, layers,
`raw`, embeddings, loadings, and pairwise matrices. Every verb is built around
that relationship.

Examples run against Scanpy's PBMC3K. Each container has a name in
`annplyr.data`, and `ap.tidy(adata)` makes the verbs methods of the object:

```{testcode}
from annplyr.data import X, X_pca, obs, raw, var

adata = ap.tidy(adata)
```

Coming from R's tidyseurat or tidySCE? See {doc}`coming_from_r`.

## Read a pipeline by its return type

Most annplyr workflows alternate between two kinds of result:

| Operation | Typical result | Use it for |
|---|---|---|
| filter, select, mutate, arrange, slice, join | AnnData | Continue the chain, or pass it to Scanpy |
| grouped AnnData verbs | `GroupedAnnData` | Reuse group-local semantics across steps |
| summarize, count, extraction | pandas object | Plot, report, model, or export |

```{testcode}
print(type(adata.filter(obs.louvain == "B cells")).__name__)
print(type(adata.group_by(obs.louvain)).__name__)
print(type(adata.count(obs.louvain)).__name__)
```

```{testoutput}
TidyAnnData
GroupedAnnData
DataFrame
```

That distinction makes materialization visible: a cell filter stays AnnData; a
long expression table is deliberately pandas.

## Axes are explicit

`obs` is the observation axis, usually cells. `var` is the feature axis, usually
genes, proteins, or peaks. Row-like operations default to cells, and `axis="var"`
switches the same verb to features.

A column name carries its **container**: `obs.percent_mito` is `obs` metadata,
`X.MS4A1` is a gene in `X`, and `raw.MS4A1` the same gene in `.raw`. One call can
filter cells by metadata while selecting genes from the matrix:

```{testcode}
subset = adata.filter(obs.percent_mito < 0.03).select(
    var.n_cells,
    X.MS4A1,
    X.CD79A,
    X.NKG7,
)

print(subset.shape)
```

```{testoutput}
(2257, 3)
```

Naming the container is what prevents a gene called `n_counts` from colliding
with the `obs` column of the same name. Strings still work for names computed in
a loop: `"percent_mito"` is an `obs` column and `"raw:MS4A1"` a gene in `.raw`.

## Containers are importable names

| Name | Container | Example |
|---|---|---|
| `obs`, `var` | metadata tables | `obs.louvain` |
| `X` | the main matrix | `X.MS4A1` |
| `raw` | the `.raw` snapshot | `raw["HLA-DRA"]` |
| any other name | a layer, an `obsm` key or a `varm` key | `X_pca[0]` |

Import the layers and keys you use: `from annplyr.data import counts, X_umap`.
Brackets cover names that are not identifiers. A mistyped layer fails at the
verb with `UnknownSourceError`, which lists the layers and keys that exist. For
a layer named like an annplyr function or a container, use `ap.layer("mean")`.

## Sources stay attached to AnnData

Expressions can read from `obs` and `var` metadata, selected features in `X` or
a named layer, the immutable snapshot in `raw`, and keyed `obsm`/`varm`
matrices. `pull()` additionally reaches `obsp`, `varp`, and tabular `uns`.

```{testcode}
annotated = adata.mutate(
    deep=obs.n_counts >= 2000,
    ms4a1_detected=raw.MS4A1 > 0,
    pc1=X_pca[0],
)

print(annotated.obs[["louvain", "deep", "ms4a1_detected", "pc1"]].head(3).round(3))
```

```{testoutput}
                      louvain  deep  ms4a1_detected    pc1
index                                                     
AAACATACAACCAC-1  CD4 T cells  True           False  5.556
AAACATTGAGCTAC-1      B cells  True            True  7.210
AAACATTGATCAGC-1  CD4 T cells  True           False  2.694
```

Three derived columns from three different containers in one call. Every
matrix-like input is read-only: `mutate()` writes to metadata and leaves `X`,
`raw`, and `X_pca` untouched.

## Ownership is predictable

Every verb returns a new object, and the input is never changed:

```{testcode}
flagged = adata.mutate(flag=ap.lit(1))
print(flagged is adata, "flag" in adata.obs, "flag" in flagged.obs)
```

```{testoutput}
False False True
```

There are no `copy=` or `inplace=` keywords. To keep a result, assign it:
`adata = adata.mutate(...)`. A verb that keeps the shape, such as `mutate()`,
copies `obs`, `var` and `uns` but shares `X`, layers and the aligned matrices
with its input, so an in-place edit of a matrix shows in both objects. A
subsetting verb on a backed file returns a view, so call `.to_memory()` before
you change the result.


## Positional identity protects duplicate names

Concatenated objects routinely carry repeated `obs_names`. annplyr tracks
integer positions for filtering, ordering, joins, and grouped operations, so
duplicate labels never collapse distinct cells:

```{testcode}
duplicated = ap.tidy(sc.datasets.pbmc3k_processed())
names = duplicated.obs_names.to_numpy().copy()
names[1] = names[0]
duplicated.obs_names = names

print(duplicated.obs_names[:3].tolist())
print(duplicated.filter(obs.louvain == "CD4 T cells").n_obs)
```

```{testoutput}
['AAACATACAACCAC-1', 'AAACATACAACCAC-1', 'AAACATTGATCAGC-1']
1144
```

Label-based selection would have returned the wrong number of cells here, or
raised. All aligned containers use the same positional subset.

## Expressions and selectors compose

`col()`, `lit()`, rank helpers, null helpers, and selectors return
`AnnplyrExpr` objects, and operators preserve them, so a domain rule can be
built once and reused:

```{testcode}
qc_rule = (obs.n_genes >= 500) & (obs.percent_mito < 0.03)

analysis = adata.filter(qc_rule).select(ap.starts_with("n_"), X.MS4A1, X.CD79A)

print(analysis.shape, analysis.obs.columns.tolist(), analysis.var_names.tolist())
```

```{testoutput}
(2117, 2) ['n_genes', 'n_counts'] ['MS4A1', 'CD79A']
```

Selectors such as `ap.starts_with()`, `ap.any_of()` and `ap.where()` choose
`obs` columns. `any_of()` tolerates an absent name; `all_of()` would raise.
`where()` is a schema selector. Its predicate receives a zero-length typed Series
and should inspect dtype, not values. Call `to_narwhals()` only when an external
API needs a raw Narwhals expression.

## Grouping can persist

Use `by=` for one summary, and `group_by()` when several operations share the
same groups. The wrapper records group keys rather than a detached copy of the
data. See {doc}`grouping` for group order, missing values, categorical keys, and
key updates.

## Projection and budgets make reads deliberate

annplyr resolves selected rows and columns before reading matrix sources, and
`max_matrix_values=` sets a hard cumulative bound that is checked before the
first read:

```{testcode}
try:
    adata.join_features(raw.MS4A1, raw.CD79A, max_matrix_values=100)

except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (join_features) would materialize 5276 matrix values, which exceeds max_matrix_values=100
```

Canonical scalar summaries then use deterministic internal chunks for dense,
CSR, CSC, and backed sources. You specify the scientific projection, not an
implementation-specific chunk size.

## Typed errors expose invalid assumptions

Invalid selectors, missing sources, incompatible axes, unsafe joins, duplicate
outputs, and over-budget reads raise annplyr-specific exceptions. Catch a
specific error when the workflow can recover; otherwise let the message expose
the failed assumption. The complete list is in {doc}`../api`.

## Design lineage

`annplyr` is inspired by [annsel](https://github.com/srivarra/annsel), which
introduced predicate-based selection for AnnData. annplyr extends that idea with
tidyverse-style mutation, summaries, grouping, joins, and rectangling for users
who want a familiar dataframe grammar inside Python single-cell workflows.
