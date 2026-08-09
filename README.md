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

## Try it on real data

Importing `annplyr` registers the `adata.ap` accessor. This runs as written on
Scanpy's PBMC3K:

```python
import scanpy as sc

import annplyr as ap

adata = sc.datasets.pbmc3k_processed()

summary = adata.ap.summarize(
    obs={"cells": ap.n(), "median_genes": ap.median("n_genes")},
    raw={"mean_MS4A1": ap.mean("MS4A1"), "mean_LYZ": ap.mean("LYZ")},
    by="louvain",
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

One call, one grouping, two sources: cell counts and depth come from `obs`,
marker means come from the 13,714 genes in `.raw`, and only the two requested
genes are ever read. `MS4A1` (CD20) marks the B cells and `LYZ` the monocytes
and dendritic cells, so the table reports biology rather than plumbing.

QC and extraction chain the same way:

```python
plot_data = (
    adata.ap.filter(obs=ap.col("percent_mito") < 0.03)
    .ap.to_tidy(obs=["louvain"], raw=["MS4A1", "CD79A"], max_matrix_values=2 * adata.n_obs)
)
print(plot_data.shape)   # (4514, 4)
```

The AnnData-returning steps preserve aligned `X`, layers, `raw`, embeddings,
loadings, and pairwise matrices. The boundary to pandas is explicit: summary and
extraction methods return DataFrames, everything else returns AnnData.

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
