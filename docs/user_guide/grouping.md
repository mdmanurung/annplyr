# Persistent grouping

Use `group_by()` when several steps should share the same sample, condition,
cell type, or feature groups. It resolves existing metadata keys for one AnnData
axis and returns a `GroupedAnnData` wrapper.

For a single summary, `summarize(..., by=...)` is shorter. Persistent grouping
earns its keep for a sequence such as *rank cells within each type, keep the
top ones, then summarize what is left*.

All examples run against Scanpy's PBMC3K, grouped by its Louvain cell-type
labels.

```{testcode}
by_type = adata.ap.group_by(obs="louvain")

print(by_type.group_vars())
print(by_type.group_keys())
```

```{testoutput}
['louvain']
             louvain
0        CD4 T cells
1            B cells
2    CD14+ Monocytes
3           NK cells
4        CD8 T cells
5  FCGR3A+ Monocytes
6    Dendritic cells
7     Megakaryocytes
```

Groups follow **first-seen** order. `CD4 T cells` and `B cells` lead because
those are the first two cells in the object, not because of an alphabetical or
frequency rule. Grouping also includes NA values, preserves categorical dtype
and order, and omits unobserved categories.

`group_data()` adds a `.rows` column of zero-based integer positions, not axis
labels, so duplicate `obs_names` never collapse two distinct cells.

## Inspect groups before transforming

Iteration yields each key and an aligned AnnData subset:

```{testcode}
sizes = [(key["louvain"], group.n_obs) for key, group in by_type]

print(sizes[:3])
```

```{testoutput}
[('CD4 T cells', 1144), ('B cells', 342), ('CD14+ Monocytes', 480)]
```

Iteration is useful for handing an aligned subset to project code. Most
workflows can stay vectorized instead: AnnData-returning grouped methods return
another grouped wrapper, while summary and count methods return pandas tables.

## Rank within group, then select

This is the pattern persistent grouping exists for. The ranking and the
selection share one grouping, and the result is only ungrouped at the end.

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
deepest = (
    adata.ap.group_by(obs="louvain")
    .mutate(obs={"depth_rank": ap.min_rank("n_counts", descending=True)})
    .slice_max(ap.col("n_counts"), n=1)
    .ungroup()
)

print(deepest.obs[["louvain", "n_counts", "depth_rank"]])
```

```{testoutput}
                            louvain  n_counts  depth_rank
index                                                    
CATACTTGGGTTAC-1        CD4 T cells    7167.0         1.0
CAGGTTGAGGATCT-1            B cells    8011.0         1.0
AACCTACTGTGAGG-1    CD14+ Monocytes    5682.0         1.0
CAGTTTACACACGT-1           NK cells    5343.0         1.0
GCCTCAACTCTTTG-1        CD8 T cells    5981.0         1.0
ATTTAGGAACCATG-1  FCGR3A+ Monocytes    5677.0         1.0
GGGCCAACCTTGGA-1    Dendritic cells    8415.0         1.0
ACGAACTGGCTATG-1     Megakaryocytes    8875.0         1.0
```

:::

:::{tab-item} scanpy + pandas

```{testcode}
positions = adata.obs.groupby("louvain", sort=False, observed=True)["n_counts"].idxmax()
deepest_ref = adata[positions].copy()

print(deepest_ref.obs[["louvain", "n_counts"]])
```

```{testoutput}
                            louvain  n_counts
index                                        
CATACTTGGGTTAC-1        CD4 T cells    7167.0
CAGGTTGAGGATCT-1            B cells    8011.0
AACCTACTGTGAGG-1    CD14+ Monocytes    5682.0
CAGTTTACACACGT-1           NK cells    5343.0
GCCTCAACTCTTTG-1        CD8 T cells    5981.0
ATTTAGGAACCATG-1  FCGR3A+ Monocytes    5677.0
GGGCCAACCTTGGA-1    Dendritic cells    8415.0
ACGAACTGGCTATG-1     Megakaryocytes    8875.0
```

:::

::::

The baseline works for "one row per group" but stops there. Keep the top 10 per
type, keep everything above the group median, or rank now and filter two steps
later, and it has to leave the AnnData object behind and reindex by hand. The
grouped wrapper keeps the ranking column *and* the aligned matrices through
every step.

## Summarize the groups

```{testcode}
print(by_type.summarize(obs={"cells": ap.n(), "mean_genes": ap.mean("n_genes")}).round(1))
```

```{testoutput}
             louvain  cells  mean_genes
0        CD4 T cells   1144       810.3
1            B cells    342       725.4
2    CD14+ Monocytes    480       866.7
3           NK cells    154       905.5
4        CD8 T cells    316       837.6
5  FCGR3A+ Monocytes    150      1228.8
6    Dendritic cells     37      1466.9
7     Megakaryocytes     15       577.3
```

Group keys always travel with their own aggregates, including when a single
call mixes metadata and matrix sources. Call `ungroup()` at the point where an
ordinary AnnData object is needed by Scanpy, serialization, or another accessor
pipeline.

## Keep group keys valid

Grouped `select`, `relocate`, and `transmute` retain grouping keys
automatically, so a selection cannot silently destroy the grouping. `rename`
and `rename_with` update the group specification instead.

```{testcode}
print(by_type.select(obs=["n_genes"]).ungroup().obs.columns.tolist())
print(by_type.rename(obs={"cell_type": "louvain"}).group_vars())
```

```{testoutput}
['louvain', 'n_genes']
['cell_type']
```

`n_genes` was the only column requested; `louvain` was retained because it is a
key, and placed first.

If a grouped `mutate()` changes a key, expressions in that call use the groups
established at its start, and the returned wrapper resolves later verbs from the
final key values. All six grouped joins run once globally and update suffixed
key names when necessary.

Grouping is axis-specific. Pass `var=` for feature groups such as chromosome or
feature type. Computed, virtual-name, cross-axis, and empty grouping
specifications raise typed selection or axis errors.
