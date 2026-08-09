# Single-cell utilities

`annplyr` ships a few narrow AnnData helpers for metadata workflows. They are
deliberately small — this is not a biological QC or annotation toolkit.

## Validate marker panels before you project

A configurable marker panel should fail loudly, not produce an empty column.
`feature_present()` reports what is present, what is missing, and what is
present under different capitalization:

```{testcode}
presence = ap.feature_present(adata, ["MS4A1", "CD79A", "ms4a1", "CD4"])

print(presence.to_frame())
```

```{testoutput}
  feature  found wrong_case_found
0   MS4A1   True              NaN
1   CD79A   True              NaN
2   ms4a1  False            MS4A1
3     CD4  False              NaN
```

Both failure modes are real here: `ms4a1` is a casing mistake, while `CD4` is
genuinely absent from PBMC3K's 1,838 highly variable genes even though CD4 T
cells are the largest cluster. Checking the panel first turns that into a
decision instead of a silently empty plot.

```{testcode}
print(presence.found_features, presence.missing_features)
```

```{testoutput}
['MS4A1', 'CD79A'] ['ms4a1', 'CD4']
```

## Summarize per sample

`sample_summary()` is `summarize()` with a sample key, for the common
per-sample table:

```{testcode}
print(
    ap.sample_summary(
        adata,
        sample="louvain",
        obs={"cells": ap.n(), "mean_genes": ap.mean("n_genes")},
    ).round(1)
)
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

`sample_meta()` is stricter: it extracts one row per sample and refuses columns
that are not constant within a sample, because those are cell-level
measurements wearing a sample-level label:

```{testcode}
try:
    ap.sample_meta(adata, sample="louvain", include=["n_genes"])
except ap.AnnplyrError as error:
    print(error)
```

```{testoutput}
Column(s) are not sample-constant for 'louvain': n_genes
```

`add_sample_meta(..., inplace=False)` returns an independent object; an
explicit in-place call validates the sample relationship before writing and
returns the identical AnnData.

## Edit and audit axis names

Name writers return independent objects unless `inplace=True` is explicit:

```{testcode}
prefixed = ap.add_name_prefix(adata, "donorA", axis="obs")

print(prefixed.obs_names[:3].tolist())
print(ap.name_duplicates(adata, axis="obs").empty)
```

```{testoutput}
['donorA_AAACATACAACCAC-1', 'donorA_AAACATTGAGCTAC-1', 'donorA_AAACATTGATCAGC-1']
True
```

Prefixing is what you do before concatenating donors; `name_duplicates()` is
what you run afterwards, and it reports counts and integer positions rather
than just a boolean.

## Store plotting palettes

Palettes are stored under Scanpy-compatible `uns` keys, so Scanpy's plotting
functions pick them up without further wiring:

```{testcode}
coloured = ap.store_palette(
    adata,
    "louvain",
    ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"],
)

print("louvain_colors" in coloured.uns)
print(ap.get_palette(coloured, "louvain")["B cells"])
```

```{testoutput}
True
#2ca02c
```

`get_palette()` returns a category-to-colour mapping, so a figure legend can be
built without indexing a bare list in category order.

Out of scope for core `annplyr`: species gene registries, mitochondrial or
ribosomal scoring, cell-cycle scoring, QC plotting wrappers, marker discovery,
and cluster annotation.
