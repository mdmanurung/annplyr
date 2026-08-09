# annplyr

`annplyr` provides tidy, dataframe-style wrangling for AnnData through the
`adata.ap` accessor. Use it to express cell and feature selection, metadata
integration, group-wise calculations, and table extraction as readable Python
pipelines while AnnData remains the aligned source of truth.

```python
import annplyr as ap

analysis = (
    adata.ap.filter(obs=ap.col("qc_pass"))
    .ap.left_join(sample_sheet, by="sample_id", relationship="many-to-one")
    .ap.mutate(
        x={
            "B_markers_detected": (ap.col("MS4A1") > 0)
            & (ap.col("CD79A") > 0),
        }
    )
)

summary = analysis.ap.summarize(
    obs={
        "cells": ap.n(),
        "fraction_B_markers": ap.mean("B_markers_detected"),
    },
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by=["condition", "cell_type"],
)
```

AnnData-returning verbs preserve axis alignment. Table-producing verbs make
the transition to pandas explicit, so the same workflow can end in a report,
plot, model matrix, or ordinary Python function.

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
