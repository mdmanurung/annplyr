---
name: annplyr
description: "Use when writing or reviewing Python AnnData wrangling with annplyr/adata.ap for single-cell workflows: tidyverse-style filter/select/mutate/summarize/count/group_by/join/pivot over obs, var, X/layers, raw, obsm, varm, obsp, varp, and tabular uns; plotnine-ready to_df/to_tidy/as_frame extraction; and alignment, sparse/backed, or matrix-materialization safety."
license: BSD-3-Clause
---

# annplyr

Use this skill when the task involves `annplyr`, AnnData wrangling, tidyverse-style verbs for single-cell data, or plot-ready extraction from AnnData.

## First Moves

1. Import `annplyr as ap`; this registers the `adata.ap` accessor.
2. Inspect `adata.n_obs`, `adata.n_vars`, `adata.obs.columns`, `adata.var.columns`, relevant `layers`, `raw`, `obsm`, `varm`, `obsp`, `varp`, and `uns` keys before writing expressions.
3. Prefer accessor verbs for AnnData-preserving operations and pandas extraction helpers only when tabular materialization is intended.
4. Keep axis semantics explicit: `obs` operations act on cells/observations, `var` operations act on features/genes.
5. For mypy, cast an incoming `AnnData` once to `annplyr.typing.AnnDataWithAnnplyr`; accessor results then retain typed `.ap` chaining.

## Core Rules

- Use `ap.col("name")`, selectors such as `ap.starts_with(...)`, and helpers such as `ap.mean(...)`, `ap.if_else(...)`, and `ap.case_when(...)` instead of pandas string queries.
- `mutate()` and `transmute()` write only `obs` or `var` metadata. Matrix-like sources (`x`, layers, `raw`, `obsm`, `varm`) are read-only expression sources.
- AnnData-returning verbs must preserve alignment across `obs`, `var`, `X`, layers, `obsm`, `varm`, `obsp`, and `varp`.
- Axis-changing verbs default to independent results. Treat `copy=False` as permission for either a view or a materialized result, never as an `is_view` guarantee.
- Same-shape writers use `inplace=False`; grouped AnnData verbs preserve grouping until `ungroup()`.
- Joins are metadata joins for `obs` or `var`; they should not silently create duplicated cells/features.
- Whole-matrix long exports are expensive. Select features explicitly unless the user clearly asks for all features, and use cumulative `max_matrix_values=` budgets.
- Backed axis selection with `copy=True` returns an in-memory result. Load backed objects before same-shape mutation.

## Reference Files

- Read `references/quickstart.md` for minimal working examples.
- Read `references/api-patterns.md` for verb patterns, grouping, joins, pivots, and extraction.
- Read `references/safety.md` before changing AnnData axes, materializing matrices, or handling backed/sparse data.

## Validation

For package changes, run the most focused relevant tests first, then at least:

```bash
pytest -q
uvx hatch run type:check
uvx hatch run type:consumer
```

For docs or skill changes, also check the Sphinx build when practical:

```bash
uvx hatch run docs:build
uvx hatch run docs:doctest
```
