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
- AnnData-safe metadata joins for `obs` and `var`, with relationship checks.
- Tidyr-style pandas extraction helpers: `pivot_longer`, `pivot_wider`,
  `nest_by`, and `unnest`.
- Group-local filtering, slicing, summaries, and count helpers.
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
