# Core concepts

`annplyr` treats AnnData as one aligned data structure. `obs` and `var` are
convenient metadata tables, but their rows are coordinated with `X`, layers,
`raw`, embeddings, loadings, and pairwise matrices. Every verb is built around
that relationship.

Examples run against Scanpy's PBMC3K.

## Read a pipeline by its return type

Most annplyr workflows alternate between two kinds of result:

| Operation | Typical result | Use it for |
|---|---|---|
| filter, select, mutate, arrange, slice, join | AnnData | Continue analysis or call `.ap` again |
| grouped AnnData verbs | `GroupedAnnData` | Reuse group-local semantics across steps |
| summarize, count, extraction | pandas object | Plot, report, model, or export |

```{testcode}
print(type(adata.ap.filter(obs=ap.col("louvain") == "B cells")).__name__)
print(type(adata.ap.group_by(obs="louvain")).__name__)
print(type(adata.ap.count(by="louvain")).__name__)
```

```{testoutput}
AnnData
GroupedAnnData
DataFrame
```

That distinction makes materialization visible: a cell filter stays AnnData; a
long expression table is deliberately pandas.

## Axes are explicit

`obs` is the observation axis, usually cells. `var` is the feature axis, usually
genes, proteins, or peaks. Row-like operations default to cells, and `axis="var"`
switches the same verb to features.

Arguments also name their **source**, so one call can filter cells by metadata
while selecting genes from the matrix:

```{testcode}
subset = adata.ap.filter(obs=ap.col("percent_mito") < 0.03).ap.select(
    var=["n_cells"],
    x=["MS4A1", "CD79A", "NKG7"],
)

print(subset.shape)
```

```{testoutput}
(2257, 3)
```

Naming the source is what prevents a gene called `n_counts` from colliding with
the `obs` column of the same name.

## Sources stay attached to AnnData

Expressions can read from `obs` and `var` metadata, selected features in `X` or
a named layer through `x=`, the immutable snapshot through `raw=`, and keyed
`obsm`/`varm` matrices. `as_frame()` additionally reaches `obsp`, `varp`, and
tabular `uns`.

```{testcode}
annotated = adata.ap.mutate(
    obs={"deep": ap.col("n_counts") >= 2000},
    raw={"ms4a1_detected": ap.col("MS4A1") > 0},
    obsm={"X_pca": {"pc1": ap.col("0")}},
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

Three derived columns from three different sources in one call. Every
matrix-like input is read-only: `mutate()` writes to metadata and leaves `X`,
`raw`, and `X_pca` untouched.

## Ownership is predictable

Axis-changing and ordering verbs return independent AnnData objects by default.
Same-shape metadata operations return a new object too, unless `inplace=True`
asks for exact input identity:

```{testcode}
copied = adata.ap.mutate(obs={"flag": ap.lit(1)})
print(copied is adata, "flag" in adata.obs)

same = adata.ap.mutate(obs={"flag": ap.lit(1)}, inplace=True)
print(same is adata, "flag" in adata.obs)
```

```{testoutput}
False False
True True
```

Validation completes before the first write, so a failing `inplace=True` call
leaves the object as it was. `copy=False` asks annplyr to avoid a guaranteed
copy, but the result may be a view or a materialized object depending on what
is safe. Do not branch on `is_view` after it.

## Positional identity protects duplicate names

Concatenated objects routinely carry repeated `obs_names`. annplyr tracks
integer positions for filtering, ordering, joins, and grouped operations, so
duplicate labels never collapse distinct cells:

```{testcode}
duplicated = sc.datasets.pbmc3k_processed()
names = duplicated.obs_names.to_numpy().copy()
names[1] = names[0]
duplicated.obs_names = names

print(duplicated.obs_names[:3].tolist())
print(duplicated.ap.filter(obs=ap.col("louvain") == "CD4 T cells").n_obs)
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
qc_rule = (ap.col("n_genes") >= 500) & (ap.col("percent_mito") < 0.03)
panel = ap.any_of(["MS4A1", "CD79A", "MISSING_GENE"])

analysis = adata.ap.filter(obs=qc_rule).ap.select(x=panel)

print(analysis.shape, analysis.var_names.tolist())
```

```{testoutput}
(2117, 2) ['MS4A1', 'CD79A']
```

`any_of()` tolerated the absent gene; `all_of()` would have raised. `where()` is
a schema selector. Its predicate receives a zero-length typed Series and should
inspect dtype, not values. Call `to_narwhals()` only when an external API needs
a raw Narwhals expression.

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
    adata.ap.to_tidy(obs=["louvain"], raw=["MS4A1", "CD79A"], max_matrix_values=100)
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
planned source(s) (to_tidy) would materialize 5276 matrix values, which exceeds max_matrix_values=100
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
