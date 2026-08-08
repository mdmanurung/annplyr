# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.3.0] - 2026-08-08

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
  raw-sample evaluator, peak-RSS runner, and v0.3 performance report.
- Dedicated v0.2-to-v0.3 migration documentation and executable Sphinx
  examples.

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
- On the frozen runner, the manifest's best improvements were 99.92% for
  matrix projection, 68.87% for grouped execution, and 39.37% for metadata
  evaluation. See `docs/development/performance-v0.3.md` for raw-sample ranges,
  caveats, and the non-primary ASV continuous regressions.

### Fixed

- Preserve aligned containers and requested order when axis labels are
  duplicated.
- Preserve nullable integer, Boolean, string, and categorical dtypes when
  grouped results are reconstructed positionally.
- Reject cumulative and later-source over-budget requests before any adapter
  reads occur.

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
