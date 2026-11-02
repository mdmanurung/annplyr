---
name: annplyr
description: "Use when writing or reviewing Python AnnData wrangling with annplyr/adata.ap for single-cell workflows: tidyverse-style filter/select/mutate/summarize/count/group_by/join/pivot over obs, var, X/layers, raw, obsm, varm, obsp, varp, and tabular uns; plotnine-ready to_wide/join_features/pull extraction; pseudobulk with aggregate_cells; and alignment, sparse/backed, or matrix-materialization safety."
license: BSD-3-Clause
---

# annplyr

Use this skill when the task involves `annplyr`, AnnData wrangling, tidyverse-style verbs for single-cell data, or plot-ready extraction from AnnData.

## First Moves

1. Import `annplyr as ap` and the containers you read, such as `from annplyr.data import obs, var, X, raw`; wrap the object with `adata = ap.tidy(adata)` so the verbs are methods. `adata.ap.<verb>` also works on a plain AnnData.
2. Inspect `adata.n_obs`, `adata.n_vars`, `adata.obs.columns`, `adata.var.columns`, relevant `layers`, `raw`, `obsm`, `varm`, `obsp`, `varp`, and `uns` keys before writing expressions.
3. Prefer accessor verbs for AnnData-preserving operations and pandas extraction helpers only when tabular materialization is intended.
4. Keep axis semantics explicit: `obs` operations act on cells/observations, `var` operations act on features/genes.
5. For mypy, cast an incoming `AnnData` once to `annplyr.typing.AnnDataWithAnnplyr`; accessor results then retain typed `.ap` chaining.

## Core Rules

- Name every column with its container: `obs.percent_mito`, `raw.MS4A1`, `X.MS4A1`, `counts.MS4A1` (any layer imported from `annplyr.data`), `X_umap[0]` (an `obsm` key), `raw["HLA-DRA"]`. Use `ap.layer("name")` for a layer whose name is reserved. Strings such as `"raw:MS4A1"` and `ap.col(...)` are for names computed in a loop. Use selectors such as `ap.starts_with(...)` (bare, they pick `obs` columns; `X[ap.starts_with(...)]` picks features) and helpers such as `ap.mean(...)`, `ap.if_else(...)` and `ap.case_when(...)` instead of pandas string queries.
- Verbs take names positionally: `select(obs.louvain, -obs.n_genes)`, `count(obs.sample)`, `group_by(obs.louvain)`. New columns are keywords: `mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)`, `summarize(cells=ap.n(), by=obs.louvain)`. A dict is for names that are not identifiers.
- `mutate()` and `transmute()` write only `obs` or `var` metadata. Matrix-like containers (`X`, layers, `raw`, `obsm`, `varm`) are read-only expression sources.
- AnnData-returning verbs must preserve alignment across `obs`, `var`, `X`, layers, `obsm`, `varm`, `obsp`, and `varp`.
- Every verb returns a new object, so assign the result. There is no `copy=` or `inplace=`. Same-shape verbs share the matrices with the input; a subsetting verb on a backed input returns the view. Grouped AnnData verbs preserve grouping until `ungroup()`.
- Joins are metadata joins for `obs` or `var`; they should not silently create duplicated cells/features.
- Whole-matrix exports are expensive. Name the features, pass `allow_all_features=True` only when the user clearly asks for all of them, and use cumulative `max_matrix_values=` budgets.
- Load backed objects with `.to_memory()` before a metadata-changing verb.

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
