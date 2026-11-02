---
title: "annplyr: tidy dataframe verbs for AnnData"
tags:
  - Python
  - AnnData
  - scverse
  - single-cell
  - data wrangling
  - tidy data
authors:
  - name: Mikhael Manurung
    orcid: "[gap: ORCID]"
    affiliation: 1
affiliations:
  - name: "[gap: affiliation]"
    index: 1
date: 27 September 2026
bibliography: paper.bib
---

# Summary

`annplyr` adds dplyr-style verbs to AnnData [@virshup2024anndata]. One call to `ap.tidy(adata)` makes the verbs methods of the object. Users then filter, select, mutate, join, group and summarize cells and features. Verbs that return AnnData keep `obs` aligned with `X`, layers, `.raw`, embeddings, loadings and pairwise matrices. Verbs that return tables give pandas DataFrames. The boundary between the two is explicit.

# Statement of need

AnnData is the standard container for single-cell data in Python, and it underlies the scverse ecosystem [@virshup2023scverse]. Cell metadata live in `obs`, a pandas DataFrame. Most analysis questions mix that metadata with expression values, for example the mean of one marker gene per cluster.

Without a helper, users pull `obs` out, join values from `X` or `.raw` by hand, and write the result back. Each step is a chance to misalign rows. A mask built on `obs` must be applied separately to `X`, layers, `.raw` and `obsm`, and one slip leaves them out of step. A merge can reorder or duplicate cells. Duplicate cell names make label-based joins ambiguous. Large matrices are easy to densify by accident.

`annplyr` is for analysts who already work with Scanpy and want pipelines that a colleague can read in one pass. It removes the manual steps above. A single expression can read metadata and matrix features. The verb applies the result to every aligned container. Mistakes surface as typed errors rather than as a plausible wrong table.

# State of the field

The idea of tidy verbs on annotated matrices is not new. The R community has used it for several years.

R users already have this style of work. dplyr and the tidyverse [@wickham2019tidyverse] provide the grammar of verbs. `tidySingleCellExperiment` [@mangiola2025tidysce] and `tidyseurat` [@mangiola2021tidyseurat] let users apply those verbs to `SingleCellExperiment` and Seurat objects. The tidyomics project collects such packages [@hutchison2024tidyomics].

Python has no equivalent for AnnData. The usual options are listed below.

- **Plain pandas or polars on `.obs`.** This works for metadata alone. It does not carry a filter over to `X`, layers or `obsm`, and it does not look inside matrices.
- **Scanpy getters.** `sc.get.obs_df` and related functions extract tables from AnnData [@wolf2018scanpy]. They extract. They do not filter, group or join, and they do not return AnnData.
- **dplyr-style libraries for dataframes.** `siuba` [@siuba] ports dplyr to pandas and SQL. Pyjanitor offers similar cleaning helpers. Both work on dataframes, so neither knows that `obs` must stay aligned with a matrix.

`annplyr` fills this gap. It builds on Narwhals [@narwhals] to evaluate expressions against dataframe backends. It does not replace Scanpy or AnnData, and it is not a quality control or annotation workflow.

# Functionality

The main verbs are listed below. Each is available as a method on a tidied object, as an `adata.ap` accessor method, and as a free function that takes the AnnData first. Every verb returns a new object.

- `filter`, `select`, `mutate`, `transmute`, `rename`, `relocate`, `arrange`, `distinct` and the `slice_*` family.
- `left_join`, `inner_join`, `semi_join` and `anti_join` for sample sheets and feature annotations. A `relationship` argument states the expected key cardinality. The default is many-to-one, so a duplicated key raises `JoinRelationshipError` instead of duplicating cells.
- `group_by` returns a persistent grouped object for ranking, filtering, mutation and summaries within groups.
- `summarize` and `count` return DataFrames.
- `to_wide`, `join_features`, `pivot_longer` and `pull` extract plot-ready tables.
- `aggregate_cells` returns a pseudobulk AnnData, and `nest` returns one object per group.

Every name carries its container, as in tidyseurat and tidySingleCellExperiment, whose verbs these follow. `from annplyr.data import obs, var, X, raw` gives the fixed containers. Any other name is a layer or an `obsm` key, so `obs.percent_mito`, `raw.MS4A1`, `counts.MS4A1` and `X_umap[0]` each read one column. Verbs take names positionally, and new columns are keywords, as in `mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)`. One expression cannot mix containers. It fails with `IncompatibleAxisError` and names both.

Several design choices target the misalignment risk. Axis operations use integer position, so duplicate cell or feature names stay distinct. Sparse matrices stay sparse during projection. Only the requested features are read. An optional `max_matrix_values` argument sets a cumulative limit, and the planner rejects a request that would exceed it before any read. AnnData-returning verbs return `TidyAnnData`, a subclass of `AnnData`, so `isinstance` checks and Scanpy calls keep working.

The package also ships helpers to check marker panels against `var_names` (`feature_present`), extract sample metadata, store plotting palettes, and reshape the returned DataFrames (`annplyr.tidyr`). It bundles an Agent Skill that users install with `annplyr-install-skills`. The documentation examples run on the PBMC3k dataset as doctests, so a stale example fails the build.

# Example

This example summarizes the processed PBMC3k dataset by cluster. Cell counts and genes detected per cell come from `obs`. Marker means come from `.raw`.

```python
import scanpy as sc

import annplyr as ap
from annplyr.data import obs, raw

adata = ap.tidy(sc.datasets.pbmc3k_processed())

summary = adata.summarize(
    cells=ap.n(),
    median_genes=ap.median(obs.n_genes),
    mean_MS4A1=ap.mean(raw.MS4A1),
    mean_LYZ=ap.mean(raw.LYZ),
    by=obs.louvain,
)

print(summary.round(2))
```

```text
             louvain  cells  median_genes  mean_MS4A1  mean_LYZ
0        CD4 T cells   1144         809.0        0.03      0.43
1            B cells    342         677.0        0.99      0.39
2    CD14+ Monocytes    480         859.0        0.04      3.55
3           NK cells    154         890.0        0.04      0.35
4        CD8 T cells    316         824.5        0.04      0.36
5  FCGR3A+ Monocytes    150        1272.0        0.06      2.24
6    Dendritic cells     37        1544.0        0.06      3.91
7     Megakaryocytes     15         364.0        0.05      0.60
```

*MS4A1* marks B cells and *LYZ* marks monocytes and dendritic cells, as expected. The call reads only the two named genes.

The result stays attached to the data it came from. Users can pass the filtered AnnData to Scanpy, or hand the summary table to a plotting library, without converting objects in between.

# Acknowledgements

[gap: funding statement and acknowledgements]

# References
