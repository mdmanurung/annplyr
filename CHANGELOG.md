# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- A deterministic, redistributable downstream AnnData fixture with exact H5AD,
  Zarr, concatenation, positional-identity, ownership, grouping, and optional
  Scanpy integration checks across stable Python 3.12-3.14 CI jobs.
- One internal dense/sparse/backed chunk plan for canonical matrix summaries,
  with deterministic row and feature batching, exact eager-equivalent scalar
  results, and shared positional group plans.
- Exported `AnnplyrAccessor`, `GroupedAnnData`, and `AnnplyrExpr` types, plus a
  typing-only `AnnDataWithAnnplyr` façade and exact consumer aliases under
  `annplyr.typing`.
- A clean downstream mypy fixture covering accessor chaining, grouping,
  expressions, joins, extraction, and generic `pipe()` returns.
- A blocking pedantic zizmor workflow and grouped, seven-day-cooldown
  Dependabot updates scoped to GitHub Actions.

### Changed

- Matrix-backed `summarize()` now evaluates canonical scalar reducers in bounded
  batches without adding a public tuning flag. Cumulative `max_matrix_values`
  validation still covers the complete multi-source request before its first
  read.
- Group keys from categorical and nullable metadata retain their pandas dtypes
  in matrix summaries.
- `group_by()` now has overloads that distinguish its no-group typed AnnData
  return from `GroupedAnnData`; public selector, source, expression, join, and
  grouped-return annotations are frozen by a generated contract.
- GitHub Actions now use release-traced immutable SHA pins, non-persistent
  checkout credentials, least-privilege token scopes, quoted matrix inputs,
  and a cache-disabled serialized trusted-publishing job.

### Performance

- Canonical wide reductions avoid materializing the complete selected matrix
  source. The internal target is 25,165,824 logical values per batch; exact
  rank and distinct reducers retain at most one reduction vector's state.
- On the fixed runner, a 100,000 by 500 dense mean summary retained effectively
  neutral timing (0.670315 s to 0.660730 s median) while reducing median
  three-process peak RSS by 40.57% (985076 KiB to 585388 KiB), with identical
  result hashes. See `docs/development/performance-issue-10.md`.

### Documentation

- Rebuilt the public documentation as a task-led first-release experience,
  with an executable sample-aware cohort vignette, expanded self-contained
  notebooks, and realistic QC, marker, join, grouping, and extraction examples
  throughout the user guide.

## [0.3.0] - 2026-08-09

### Added

- Persistent `GroupedAnnData` pipelines with `ungroup`, key-aware
  select/rename/relocate/transmute, all six global joins, wrapper-aware
  `pipe`, and deterministic first-seen group plans.
- `AnnplyrExpr` expression metadata and `to_narwhals()` interoperation for raw
  Narwhals consumers.
- Dense, sparse, DataFrame, AnnData-view, and backed source adapters with
  positional row/column projection and requested-order restoration.
- Two-phase matrix request planning with cumulative `max_matrix_values=`
  budgets across filtering, ordering, mutation, summaries, extraction, and
  grouped counterparts.
- A mechanically checked public API contract, fixed ASV benchmark manifest,
  raw-result normalizer/evaluator, peak-RSS runner, and v0.3 performance report.
- Executable Sphinx examples for ownership, grouping, matrix projection,
  extraction, and typed errors.

### Changed

- Axis-changing and axis-ordering verbs now default to `copy=True`; explicit
  `copy=False` may return either a view or a materialized object.
- Same-shape public mutations use `inplace=False` instead of `copy=`, validate
  before writing, and return the identical object for `inplace=True`.
- `transmute()` is always independent and no longer accepts ownership flags.
- Filtering, arranging, distinct selection, slicing, and joins reconstruct
  AnnData by integer position, preserving duplicate axis-label identity.
- Grouped AnnData-returning verbs preserve grouping, grouping columns are
  retained deterministically, and key mutations regroup subsequent calls.
- `group_data()[".rows"]` now stores zero-based integer positions. NA groups
  are included, categorical order is preserved, and unobserved levels drop.
- `where()` is schema/dtype-only. Raw Narwhals expressions remain supported as
  opaque expressions and receive conservative full-source budget charges.
- Backed `copy=True` subsetting returns an independent in-memory object;
  selected sparse exports preserve sparse formats where the return type allows.
- Right/full joins reject right-only AnnData axis records instead of silently
  constructing cells or features without aligned matrix data.

### Performance

- Project selected matrix rows/columns before pandas conversion and reuse
  metadata frames for safe expression sequences.
- Batch proven-independent assignments and skip no-op metadata writes without
  changing sequential mutation semantics.
- Avoid Python bucket construction for persistent group plans, aggregate
  grouped counts directly, and restrict filtering joins to their key columns.
- Build observation-major long tables directly instead of melting and then
  reordering feature-major output.
- On the frozen runner, the manifest's best improvements were 99.92% for
  matrix projection, 68.87% for grouped execution, and 39.37% for metadata
  evaluation. Targeted post-freeze measurements resolve the documented
  high-cardinality, filtering-join, and chained-extraction regressions. See
  `docs/development/performance-v0.3.md` for exact ranges and caveats.

### Fixed

- Preserve aligned containers and requested order when axis labels are
  duplicated.
- Preserve nullable integer, Boolean, string, and categorical dtypes when
  grouped results are reconstructed positionally.
- Reject cumulative and later-source over-budget requests before any adapter
  reads occur.
- Preserve observation-major ordering in `to_tidy()` and `pivot_longer()`.
- Reject empty or multi-key `pull()` source mappings with `UnknownSourceError`
  instead of ignoring keys or leaking `StopIteration`.
- Align `nest_by()`, join, pipe, and `pull()` documentation with their actual
  return, ownership, and typed-failure contracts.

### Documentation

- Updated ownership, grouping, join, source, backed/sparse, budget, and typed
  failure guidance throughout the README, API reference, user guide, bundled
  skill, and affected notebook cells.

## [0.2.0] - 2026-07-15

### Added

- Roadmap and GitHub issue backlog for tidyverse/scverse-grade development.
- Scverse-style documentation, release, coverage, and community scaffolding.
- Advanced expression helpers for distinct counts, first/last/nth, lead/lag,
  null handling, ranks, and cumulative calculations.
- Tidyselect-aware helpers `across`, `pick`, `if_any`, and `if_all`, plus
  additional ranking and logic helpers including `ntile`, `max_rank`,
  `average_rank`, `cummean`, `cumany`, `cumall`, `near`, `case_match`, and
  `recode`.
- AnnData-safe metadata joins for `obs` and `var`, with relationship checks.
- Tidyr-style pandas extraction helpers: `pivot_longer`, `pivot_wider`,
  `nest_by`, `nest`, `unnest`, `unnest_longer`, `unnest_wider`, `chop`,
  `unchop`, `pack`, `unpack`, `hoist`, `drop_na`, `fill`, `separate`,
  `separate_rows`, `extract`, and `unite`.
- Group-local filtering, arranging, distinct selection, slicing, summaries,
  weighted counts, and add-count/tally helpers.
- Sparse-preserving matrix frames, pairwise `obsp`/`varp` source frames, and
  controlled `uns` frame extraction.
- Raw matrix read support for filtering, arranging, metadata mutation,
  summaries, and pandas exports.
- Public `adata.ap.as_frame()` extraction for `obs`, `var`, `x`, `raw`, `obsm`,
  `varm`, `obsp`, `varp`, and tabular `uns` sources.
- `max_matrix_values=` guards for matrix-materializing export helpers.
- Bundled Claude Code and Codex skill for AnnData wrangling guidance.

### Changed

- Internal expression reads no longer depend on `AnnData.to_df()`.
- Mutating verbs now raise package errors for backed AnnData objects.
- Long matrix exports require explicit feature selection unless full
  materialization is requested.

### Fixed

- Grouped `mutate` now uses positional indexing to correctly handle AnnData
  with duplicate `obs_names` or `var_names`.
- `relocate` raises `SelectionError` when the `before`/`after` anchor column
  is among the columns being moved.
- `unnest` preserves inner column names from the first non-empty nested frame,
  so an all-empty result still has the correct schema.
- `separate` always treats `sep` as a regular expression; NA source values now
  produce `[pd.NA] * len(into)` instead of splitting the string `"None"`.

### Documentation

- Added tutorials and improved all vignettes.
- Reorganized documentation navigation.
- Added design-lineage section crediting annsel as inspiration and framing
  tidyverse R users as the target audience.

## [0.1.0]

### Added

- Initial V1 package scaffold with an `adata.ap` accessor.
- Core verbs for filtering, selecting, arranging, slicing, mutating,
  summarizing, counting, pulling, and plot-friendly extraction.
