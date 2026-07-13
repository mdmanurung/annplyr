---
name: annplyr
description: Use when working with AnnData wrangling in Python using annplyr, especially tidyverse-style filter/select/mutate/summarize/join/pivot workflows over obs, var, X/layers, obsm, and varm.
---

# annplyr

Use this skill when the task involves `annplyr`, AnnData wrangling, tidyverse-style verbs for single-cell data, or plot-ready extraction from AnnData.

## First Moves

1. Import `annplyr as ap`; this registers the `adata.ap` accessor.
2. Inspect `adata.n_obs`, `adata.n_vars`, `adata.obs.columns`, `adata.var.columns`, relevant `layers`, `obsm`, and `varm` keys before writing expressions.
3. Prefer accessor verbs for AnnData-preserving operations and pandas extraction helpers only when tabular materialization is intended.
4. Keep axis semantics explicit: `obs` operations act on cells/observations, `var` operations act on features/genes.

## Core Rules

- Use `ap.col("name")`, selectors such as `ap.starts_with(...)`, and helpers such as `ap.mean(...)`, `ap.if_else(...)`, and `ap.case_when(...)` instead of pandas string queries.
- `mutate()` and `transmute()` write only `obs` or `var` metadata. Matrix-like sources (`x`, `layers`, `obsm`, `varm`) are read-only expression sources.
- AnnData-returning verbs must preserve alignment across `obs`, `var`, `X`, layers, `obsm`, `varm`, `obsp`, and `varp`.
- Joins are metadata joins for `obs` or `var`; they should not silently create duplicated cells/features.
- Whole-matrix long exports are expensive. Select features explicitly unless the user clearly asks for all features.
- Backed AnnData objects are read-only for mutating verbs unless the object is intentionally loaded into memory first.

## Reference Files

- Read `references/quickstart.md` for minimal working examples.
- Read `references/api-patterns.md` for verb patterns, grouping, joins, pivots, and extraction.
- Read `references/safety.md` before changing AnnData axes, materializing matrices, or handling backed/sparse data.

## Validation

For package changes, run the most focused relevant tests first, then at least:

```bash
pytest -q
uvx hatch run type:check
```

For docs or skill changes, also check the Sphinx build when practical:

```bash
uvx hatch run docs:build
```

