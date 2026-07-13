# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

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

### Changed

- Internal expression reads no longer depend on `AnnData.to_df()`.
- Mutating verbs now raise package errors for backed AnnData objects.
- Long matrix exports require explicit feature selection unless full
  materialization is requested.

## [0.1.0]

### Added

- Initial V1 package scaffold with an `adata.ap` accessor.
- Core verbs for filtering, selecting, arranging, slicing, mutating,
  summarizing, counting, pulling, and plot-friendly extraction.
