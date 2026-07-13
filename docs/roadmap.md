# Roadmap

This roadmap tracks the work needed to move `annplyr` from the initial V1
package into a tidyverse- and scverse-grade AnnData wrangling library.

## Principles

- Preserve AnnData alignment. Any verb returning an `AnnData` object must keep
  `obs`, `var`, `X`, layers, `obsm`, `varm`, `obsp`, and `varp` consistent.
- Avoid silent materialization. Sparse and backed arrays must not be fully
  densified in internal verb paths unless the user explicitly requests a
  pandas export.
- Prefer tidy semantics where they are safe for AnnData. When tidyverse behavior
  would break axis alignment or duplicate rows, `annplyr` should fail with a
  clear package error.
- Keep storage policy explicit. View, copy, in-place mutation, and backed-mode
  behavior must be documented and tested.
- Make the public contract testable. Every new verb or helper needs behavior
  tests, invalid-input tests, and docs before release.

## Milestone 1: Public Contract And Errors

Define the current API precisely before broadening it.

Issue-sized tasks:

- Add package-level exceptions: `AnnplyrError`, `SelectionError`,
  `UnknownColumnError`, `UnknownSourceError`, `DuplicateNameError`,
  `NameRepairError`, `IncompatibleAxisError`, `SizeMismatchError`, and
  `JoinRelationshipError`.
- Replace public `ValueError` paths with typed errors where the failure belongs
  to the annplyr contract.
- Document source and axis rules for `obs`, `var`, `X`, layers, `obsm`, and
  `varm`.
- Document view, copy, and in-place behavior for each accessor verb.
- Expand the API reference to include the accessor, grouped wrapper, helpers,
  and extraction functions.

Acceptance gates:

- Public invalid-input paths are covered with `pytest.raises(..., match=...)`.
- The API reference builds with warnings treated as errors.
- README and docs describe the supported V1 surface without promising future
  behavior.

## Milestone 2: Selector And Expression Engine

Bring selection and expression semantics closer to dplyr/tidyselect while
remaining axis-aware.

Issue-sized tasks:

- Add strict selector resolution with source context.
- Add selectors: `all_of`, `any_of`, `where`, `last_col`, and `num_range`.
- Add selector algebra for union, intersection, and exclusion.
- Add expression helpers: `n_distinct`, `first`, `last`, `nth`, `lead`, `lag`,
  `coalesce`, `na_if`, `replace_na`, `is_na`, rank helpers, cumulative helpers,
  `across`, `pick`, `if_any`, and `if_all`.
- Define name repair and duplicate-name behavior.

Acceptance gates:

- Selection preserves requested order and fails clearly on missing strict names.
- Virtual index names never overwrite real metadata columns.
- Dense, CSR, and CSC fixtures produce equivalent selector results.

## Milestone 3: Core Dplyr Verbs

Add the missing dataframe verbs that can be implemented without corrupting
AnnData alignment.

Issue-sized tasks:

- Add `rename`, `rename_with`, and `relocate`.
- Add `distinct`, `transmute`, `tally`, and `add_count`.
- Add metadata joins for `obs` and `var`: `left_join`, `inner_join`,
  `right_join`, `full_join`, `semi_join`, and `anti_join`.
- Add join cardinality checks: `relationship`, `multiple`, `unmatched`,
  `na_matches`, and suffix handling.
- Extend `select` to supported matrix-adjacent sources where safe.

Acceptance gates:

- Joins cannot silently reorder, duplicate, or drop AnnData axes without a
  documented and tested policy.
- Schema-changing verbs preserve pandas dtypes and categorical metadata when
  possible.
- Existing V1 tests keep passing.

## Milestone 4: Grouped Semantics

Make grouping behavior deterministic and consistent across verbs.

Issue-sized tasks:

- Add `ungroup`, `group_vars`, `group_keys`, and `group_data`.
- Add `add` and `drop` arguments to `group_by`.
- Make `filter`, `mutate`, `slice_*`, `arrange`, `distinct`, `count`, and
  `add_count` group-aware where semantics are well-defined.
- Preserve grouping columns during `select` unless explicitly dropped through a
  controlled policy.
- Define NA group handling and categorical order behavior.

Acceptance gates:

- Grouped summaries are deterministic for strings, categoricals, and NA groups.
- Grouped `mutate` supports group-local expressions such as `row_number`.
- Grouped behavior is covered for both `obs` and `var` axes.

## Milestone 5: Scverse-Grade Internals

Refactor storage and namespace handling without changing the accessor name.

Issue-sized tasks:

- Register the accessor with `anndata.register_anndata_namespace("ap")`.
- Add an `Axis` abstraction and positional indexer normalization.
- Add a source registry for `obs`, `var`, `X`, layers, `obsm`, `varm`, `obsp`,
  `varp`, and controlled `uns` access.
- Add dense/sparse/backed matrix adapters with column projection, chunked
  reductions, and explicit densification budgets.
- Add copy and materialization policy helpers.
- Add dtype-preserving assignment helpers.

Acceptance gates:

- Internal matrix predicates do not call whole-object `adata.to_df()` or
  `.toarray()` unless explicitly in an export path.
- Backed read-only verbs work; backed mutating verbs either materialize
  explicitly or raise a typed error.
- AnnData-aligned containers remain shape-consistent after every tested verb.

## Milestone 6: Tidyr And Plot Extraction

Expand long/wide extraction while protecting users from accidental large
materializations.

Issue-sized tasks:

- Add general `pivot_longer` for metadata and selected matrix features.
- Add controlled `pivot_wider` or `from_tidy` for reconstructable tidy data.
- Add `nest`, `nest_by`, and `unnest` where they return pandas data frames or
  grouped wrappers without corrupting AnnData.
- Expand `to_df` and `to_tidy` docs with plotnine-ready examples.

Acceptance gates:

- Long exports have stable column names and deterministic row order.
- Large matrix reshaping requires explicit feature selection or an override.
- Plot-oriented examples run in documentation builds without executing heavy
  notebooks by default.

## Milestone 7: Scverse Infrastructure And Release

Complete the package infrastructure expected from a public scverse-style
project.

Issue-sized tasks:

- Add `.codecov.yaml`, GitHub Pages docs workflow, release workflow,
  and issue templates.
- Add `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`,
  `MAINTAINERS.md`, and `CITATION.cff`.
- Add `src/annplyr/py.typed` and type-check configuration.
- Expand README badges, install instructions, quickstart, citation, and support
  status.
- Add tutorial notebooks for getting started and plot-ready tables.
- Configure release checks for build, docs, coverage, and PyPI trusted
  publishing.

Acceptance gates:

- `prek run --all-files` passes.
- `pytest` and coverage thresholds pass locally and in CI.
- `sphinx-build -M html docs docs/_build -W` passes.
- `uv build` and `twine check --strict dist/*` pass.
- The first public release has changelog, tag, artifact, and citation metadata.

## Non-Goals Before 1.0

- Direct ggplot wrappers.
- Full lazy execution across arbitrary backends.
- Lossy joins that duplicate cells or genes without explicit user opt-in.
- Whole-matrix long exports by default on large sparse or backed objects.
- Reimplementing AnnData storage semantics outside AnnData-native subsetting.
