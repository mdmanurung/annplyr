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

The current API includes:

- row/feature verbs: `filter`, `select`, `rename`, `rename_with`, `relocate`,
  `arrange`, `distinct`, `slice`, `slice_head`, `slice_tail`, `slice_min`,
  `slice_max`, and `slice_sample`;
- mutation and summaries: `mutate`, `transmute`, `group_by`, `summarize`,
  `summarise`, weighted `count`/`tally`, `add_count`, and `add_tally`;
- AnnData-safe metadata joins: `left_join`, `inner_join`, `right_join`,
  `full_join`, `semi_join`, and `anti_join`;
- extraction helpers: `pull`, `to_df`, `to_tidy`, `pivot_longer`,
  `pivot_wider`, `nest_by`, `nest`, `unnest`, `unnest_longer`,
  `unnest_wider`, `chop`, `unchop`, `pack`, `unpack`, `hoist`,
  `separate`, `separate_rows`, `extract`, `unite`, `drop_na`, `fill`,
  and `pipe`;
- expression helpers such as `col`, `lit`, `desc`, `between`, `if_else`,
  `case_when`, `case_match`, `recode`, `near`, `row_number`, `lead`, `lag`,
  `across`, `pick`, `if_any`, `if_all`, rank helpers including `ntile`,
  cumulative helpers, and compact aggregation helpers.

AnnData-returning verbs preserve AnnData axis alignment. Joins can enrich or
subset `obs`/`var` metadata, but they raise `JoinRelationshipError` instead of
silently adding or duplicating cells or features. Matrix-long exports require an
explicit feature selection by default; pass `allow_all_features=True` when a
whole-matrix materialization is intentional. Mutating verbs raise an
`AnnplyrError` for backed AnnData objects unless you first load them into memory.

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
