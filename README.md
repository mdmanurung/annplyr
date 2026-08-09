# annplyr

[![Test](https://github.com/mdmanurung/annplyr/actions/workflows/test.yaml/badge.svg)](https://github.com/mdmanurung/annplyr/actions/workflows/test.yaml)
[![Docs](https://github.com/mdmanurung/annplyr/actions/workflows/docs.yaml/badge.svg)](https://github.com/mdmanurung/annplyr/actions/workflows/docs.yaml)
[![codecov](https://codecov.io/gh/mdmanurung/annplyr/branch/main/graph/badge.svg)](https://codecov.io/gh/mdmanurung/annplyr)

`annplyr` brings tidy, dataframe-style verbs to AnnData. It lets you filter
cells, select genes, join sample metadata, derive annotations, summarize
groups, and build plot-ready tables without separating `obs` from the matrices
and aligned containers that give it meaning.

```bash
pip install annplyr
```

Importing `annplyr` registers the `adata.ap` accessor:

```python
import annplyr as ap
```

## A cohort workflow

Suppose `adata` contains cells from several samples and `sample_sheet` contains
one row per sample. A complete wrangling pipeline can stay close to the
scientific question:

```python
analysis = (
    adata.ap.filter(
        obs=(
            ap.col("total_counts") >= 1_000,
            ap.col("pct_counts_mt") < 10,
        )
    )
    .ap.left_join(sample_sheet, by="sample_id", relationship="many-to-one")
    .ap.mutate(
        x={
            "B_markers_detected": (ap.col("MS4A1") > 0)
            & (ap.col("CD79A") > 0),
        },
    )
)

cohort_summary = analysis.ap.summarize(
    obs={
        "cells": ap.n(),
        "mean_counts": ap.mean("total_counts"),
        "fraction_B_markers": ap.mean("B_markers_detected"),
    },
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by=["condition", "cell_type"],
)

plot_data = analysis.ap.to_tidy(
    obs=["sample_id", "condition", "cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=2 * analysis.n_obs,
)
```

The AnnData-returning steps preserve aligned `X`, layers, `raw`, embeddings,
loadings, and pairwise matrices. The summary and plotting boundary is explicit:
those methods return pandas objects.

## Why use annplyr?

- **Work across AnnData sources.** Expressions can read metadata, selected
  features in `X` or a layer, `raw`, embeddings, and loadings.
- **Keep grouping in the pipeline.** `group_by()` returns a persistent grouped
  wrapper for within-sample ranking, filtering, mutation, and summaries.
- **Preserve positional identity.** Duplicate cell or feature names remain
  distinct because axis operations are position-based.
- **Make materialization deliberate.** Matrix projection happens before pandas
  conversion, and `max_matrix_values=` places a cumulative bound on reads.
- **Compose with Python.** AnnData results keep the `.ap` accessor; extracted
  tables work directly with pandas and plotting libraries.

`annplyr` is a wrangling layer, not a replacement for Scanpy, AnnData, or a
biological QC and annotation workflow.

## Learn more

- [Quickstart](https://mdmanurung.github.io/annplyr/quickstart.html)
- [End-to-end tutorials](https://mdmanurung.github.io/annplyr/tutorials.html)
- [User guide](https://mdmanurung.github.io/annplyr/user_guide/concepts.html)
- [API reference](https://mdmanurung.github.io/annplyr/api.html)

The package also bundles an annplyr Agent Skill. Install or refresh it with
`annplyr-install-skills --agent codex` or
`annplyr-install-skills --agent claude --force`.

For development setup and checks, see
[CONTRIBUTING.md](https://github.com/mdmanurung/annplyr/blob/main/CONTRIBUTING.md).
Citation metadata is available in `CITATION.cff`.
