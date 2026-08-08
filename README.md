# annplyr

[![Test](https://github.com/mdmanurung/annplyr/actions/workflows/test.yaml/badge.svg)](https://github.com/mdmanurung/annplyr/actions/workflows/test.yaml)
[![Docs](https://github.com/mdmanurung/annplyr/actions/workflows/docs.yaml/badge.svg)](https://github.com/mdmanurung/annplyr/actions/workflows/docs.yaml)
[![codecov](https://codecov.io/gh/mdmanurung/annplyr/branch/main/graph/badge.svg)](https://codecov.io/gh/mdmanurung/annplyr)

`annplyr` provides tidy, dataframe-style wrangling for AnnData through the
`adata.ap` accessor. Version 0.3 focuses on safe axis ownership, positional
alignment, projected matrix reads, persistent grouping, and explicit
materialization budgets.

## Installation

```bash
pip install annplyr
```

Importing the package registers the accessor:

```python
import annplyr as ap
```

## Persistent grouping

AnnData-returning grouped verbs keep their grouping until `ungroup()`:

```python
result = (
    adata.ap.group_by(obs="batch")
    .filter(obs=ap.col("n_counts") > 1_000)
    .mutate(obs={"within_batch": ap.row_number()})
    .select(obs=["batch", "cell_type", "within_batch"])
    .ungroup()
)
```

Axis-changing verbs return independent AnnData objects by default. Pass
`copy=False` only when a view or materialized non-copy result is acceptable;
the contract does not guarantee `is_view`. Same-shape metadata operations use
`inplace=True` for exact in-place identity.

## Sources and materialization

Expressions can address `obs`, `var`, selected `X` or layer features, `raw`,
`obsm`, and `varm`. Extraction also supports `obsp`, `varp`, and tabular `uns`
through `as_frame()`. Select matrix columns explicitly and set a cumulative
read budget when materialization size matters:

```python
plot_data = adata.ap.to_tidy(
    obs=["cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=2 * adata.n_obs,
)
```

The planner rejects an over-budget multi-source request before its first matrix
read. Sparse wide exports preserve pandas sparse columns; backed axis
operations with `copy=True` return independent in-memory selections.

## Documentation

- [Quickstart](https://mdmanurung.github.io/annplyr/quickstart.html)
- [User guide](https://mdmanurung.github.io/annplyr/user_guide/concepts.html)
- [API reference](https://mdmanurung.github.io/annplyr/api.html)
- [v0.2 to v0.3 migration guide](https://mdmanurung.github.io/annplyr/migration-v0.3.html)

The package also bundles an annplyr Agent Skill. Install or refresh it with
`annplyr-install-skills --agent codex` or
`annplyr-install-skills --agent claude --force`.

## Development

```bash
pytest -q
uvx hatch run type:check
uvx hatch run docs:build
uvx hatch run docs:doctest
```

Citation metadata is available in `CITATION.cff`.
