# annplyr

`annplyr` provides tidy, dataframe-style wrangling for AnnData through the
`adata.ap` accessor. Use it to express cell and feature selection, metadata
integration, group-wise calculations, and table extraction as readable Python
pipelines while AnnData remains the aligned source of truth.

Every example in this documentation runs against Scanpy's PBMC3K during the
build, so the numbers you read are the numbers the code produces:

```{testcode}
summary = adata.ap.summarize(
    obs={"cells": ap.n(), "median_genes": ap.median("n_genes")},
    raw={"mean_MS4A1": ap.mean("MS4A1"), "mean_LYZ": ap.mean("LYZ")},
    by="louvain",
)

print(summary.round(2))
```

```{testoutput}
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

Cell counts and depth come from `obs`, marker means from the 13,714 genes in
`.raw`, one grouping, one call — and only the two requested genes are read.

AnnData-returning verbs preserve axis alignment. Table-producing verbs make the
transition to pandas explicit, so the same workflow can end in a report, plot,
model matrix, or ordinary Python function.

## Start with a task

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Analyze a small cohort
:link: quickstart
:link-type: doc

Filter cells, attach sample metadata, derive a marker score, summarize groups,
and extract a plot-ready table.
:::

:::{grid-item-card} Follow a complete workflow
:link: tutorials
:link-type: doc

Work through self-contained notebooks for cohort wrangling, the broader verb
surface, and plotting-table preparation.
:::

:::{grid-item-card} Work safely with AnnData
:link: user_guide/scverse_safety
:link-type: doc

Understand ownership, sparse and backed sources, alignment, and bounded
materialization.
:::

:::{grid-item-card} Look up an operation
:link: api
:link-type: doc

Find accessor methods, expressions, selectors, return types, and typed errors.
:::

::::

```{toctree}
:caption: Get Started
:maxdepth: 1

installation
quickstart
tutorials
```

```{toctree}
:caption: User Guide
:maxdepth: 1

user_guide/concepts
user_guide/verbs
user_guide/grouping
user_guide/joins
user_guide/tidy_tables
user_guide/scverse_safety
user_guide/utilities
```

```{toctree}
:caption: Reference
:maxdepth: 1

api
```

```{toctree}
:caption: Contributing and Maintenance
:maxdepth: 1

contributing
development/skills
development/ci-security
template_usage
```

```{toctree}
:caption: Project
:maxdepth: 1

roadmap
changelog
references
```

```{toctree}
:hidden:

development/api-contract-v0.3
development/public-typing-contract
development/performance-v0.3
development/performance-issue-10
```
