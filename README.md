# annplyr

[![Test](https://github.com/mdmanurung/annplyr/actions/workflows/test.yaml/badge.svg)](https://github.com/mdmanurung/annplyr/actions/workflows/test.yaml)
[![Docs](https://github.com/mdmanurung/annplyr/actions/workflows/docs.yaml/badge.svg)](https://github.com/mdmanurung/annplyr/actions/workflows/docs.yaml)
[![codecov](https://codecov.io/gh/mdmanurung/annplyr/branch/main/graph/badge.svg)](https://codecov.io/gh/mdmanurung/annplyr)
[![Documentation](https://readthedocs.org/projects/annplyr/badge/?version=latest)](https://annplyr.readthedocs.io/)

`annplyr` provides tidy, dataframe-style wrangling for AnnData objects through
an `adata.ap` accessor. It is designed for single-cell workflows where metadata,
expression matrices, layers, embeddings, and loadings should be queried together
without losing AnnData alignment.

## Installation

```bash
pip install annplyr
```

For development:

```bash
git clone https://github.com/mdmanurung/annplyr
cd annplyr
uv sync --all-extras
```

## Quickstart

```python
import annplyr as ap

filtered = adata.ap.filter(obs=ap.col("batch") == "A", x=ap.col("GeneA") > 0)
plot_data = filtered.ap.to_df(obs=["batch"], x=["GeneA"])
```

The current V1 API includes:

- row/feature verbs: `filter`, `select`, `arrange`, `slice`, `slice_head`,
  `slice_tail`, `slice_min`, `slice_max`, and `slice_sample`;
- mutation and summaries: `mutate`, `group_by`, `summarize`, `summarise`, and
  `count`;
- extraction helpers: `pull`, `to_df`, `to_tidy`, and `pipe`;
- expression helpers such as `col`, `lit`, `desc`, `between`, `if_else`,
  `case_when`, `row_number`, and compact aggregation helpers.

See `docs/roadmap.md` for the tidyverse/scverse-grade development plan.

## Development

This repository follows the scverse cookiecutter style for local tooling:

```bash
pytest -q
uvx hatch run docs:build
prek run --all-files
uv build
```

The package metadata is managed through `pyproject.toml`, with Hatch
environments for tests and documentation.

## Citation

Citation metadata is available in `CITATION.cff`.
