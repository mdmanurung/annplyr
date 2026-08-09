# AnnData-safe joins

Use joins to attach a sample sheet, donor covariates, clinical annotations, or
feature metadata to AnnData. annplyr enriches or subsets `obs` or `var` without
manufacturing aligned matrix rows, and tracks each axis by integer position, so
duplicate labels stay distinct.

The examples attach a lineage annotation to PBMC3K's Louvain labels — the
everyday case of mapping fine-grained cluster names onto a coarser level.

## Attach an annotation table

::::{tab-set}

:::{tab-item} annplyr

```{testcode}
import pandas as pd

lineage = pd.DataFrame(
    {
        "louvain": [
            "CD4 T cells",
            "CD8 T cells",
            "NK cells",
            "B cells",
            "CD14+ Monocytes",
            "FCGR3A+ Monocytes",
            "Dendritic cells",
            "Megakaryocytes",
        ],
        "lineage": [
            "Lymphoid",
            "Lymphoid",
            "Lymphoid",
            "Lymphoid",
            "Myeloid",
            "Myeloid",
            "Myeloid",
            "Megakaryocytic",
        ],
    }
)

annotated = adata.ap.left_join(lineage, by="louvain", axis="obs")

print(annotated.ap.count(by="lineage", sort=True))
print(annotated.n_obs == adata.n_obs, "lineage" in adata.obs)
```

```{testoutput}
          lineage     n
0        Lymphoid  1956
1         Myeloid   667
2  Megakaryocytic    15
True False
```

:::

:::{tab-item} pandas

```{testcode}
merged = adata.obs.merge(lineage, on="louvain", how="left")

print(merged["lineage"].value_counts().rename_axis("lineage").reset_index(name="n"))
```

```{testoutput}
          lineage     n
0        Lymphoid  1956
1         Myeloid   667
2  Megakaryocytic    15
```

:::

::::

The pandas merge produces the right column, but it produces a **detached
DataFrame**: its index is a fresh `RangeIndex`, so assigning it back with
`adata.obs = merged` silently breaks the link to `obs_names` unless you restore
the index yourself. That is the class of mistake the join verb removes — the
result is an AnnData whose every aligned container still matches.

Joins return an independent AnnData by default, the left axis order is stable,
and the original object is untouched. Pass `copy=False` only when either a view
or a materialized result is acceptable.

## Restrict an analysis cohort

Inner, semi, and anti joins subset an axis while preserving the original
left-side order. A semi join is the natural way to say "keep the cells an
external table declares eligible".

```{testcode}
lymphoid = lineage.loc[lineage["lineage"] == "Lymphoid", ["louvain"]]

print(adata.ap.semi_join(lymphoid, by="louvain", axis="obs").shape)
print(adata.ap.anti_join(lymphoid, by="louvain", axis="obs").ap.count(by="louvain", sort=True))
```

```{testoutput}
(1956, 1838)
             louvain    n
0    CD14+ Monocytes  480
1  FCGR3A+ Monocytes  150
2    Dendritic cells   37
3     Megakaryocytes   15
```

Neither adds a column: a semi join filters, an anti join takes the complement.

## Annotate features

The same verbs work on the feature axis. Axis names are labels rather than
columns, so promote them into `var` first — which is itself an ordinary
`mutate()`:

```{testcode}
panel = pd.DataFrame(
    {
        "gene": ["MS4A1", "CD79A", "NKG7", "LYZ"],
        "marker_of": ["B cells", "B cells", "NK cells", "Monocytes"],
    }
)

labelled = adata.ap.mutate(var={"gene": ap.var_names})
with_panel = labelled.ap.left_join(panel, by="gene", axis="var")

print(with_panel.var.loc[["MS4A1", "NKG7", "CD79A"]])
```

```{testoutput}
       n_cells   gene marker_of
index                          
MS4A1      423  MS4A1   B cells
NKG7       816   NKG7  NK cells
CD79A      426  CD79A   B cells
```

`LYZ` is not in the 1,838 highly variable genes, so only three of the four
panel rows find a match — and the unmatched one is dropped rather than
inventing a feature.

## Make cardinality assumptions executable

`relationship`, `multiple`, `unmatched`, `na_matches`, and `suffixes` turn join
assumptions into checks. The default relationship is `many-to-one`, so a
duplicated key in the annotation table is an error rather than silent cell
duplication:

```{testcode}
duplicated = pd.DataFrame(
    {"louvain": ["CD4 T cells", "CD4 T cells"], "lineage": ["Lymphoid", "T cell"]}
)

try:
    adata.ap.left_join(duplicated, by="louvain", axis="obs")
except ap.JoinRelationshipError as error:
    print(error)
```

```{testoutput}
join found multiple right-hand matches for at least one key
```

A pandas merge here would have returned 1,144 extra rows, which cannot be
written back to an AnnData with 2,638 cells.

A right or full join cannot add a right-only cell or feature, because no
aligned `X`, layer, embedding, or pairwise record exists for it:

```{testcode}
try:
    adata.ap.right_join(pd.DataFrame({"louvain": ["Basophils"]}), by="louvain", axis="obs")
except ap.JoinRelationshipError as error:
    print(error)
```

```{testoutput}
right_join has unmatched axis records
```

Duplicate right keys, unmatched policies, invalid suffixes, and cross-axis
requests fail the same way, with typed package errors rather than a partly
fabricated AnnData.

Grouped joins execute globally and return a grouped wrapper when the resulting
metadata still contains valid grouping keys; suffixing a key updates the stored
group specification deterministically.
