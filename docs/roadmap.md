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

## Current Hardening Status

Version 0.3 completes the public contract, ownership, grouping, projection,
tidy extraction, and release-infrastructure work described by Milestones 1-7:

- `adata.ap` is registered through the AnnData namespace mechanism.
- Axis verbs use positional identity and independent-object defaults so `X`,
  layers, `raw`, `obsm`, `varm`, `obsp`, and `varp` stay aligned even when axis
  labels are duplicated.
- Persistent grouping uses one first-seen positional group plan, retains keys,
  preserves extension/categorical dtypes, and exposes zero-based `.rows`.
- Expression sources cover `obs`, `var`, selected `X`/layers, `raw`, `obsm`,
  and `varm`; `AnnplyrExpr` metadata supports safe projection and batching.
- `as_frame()` provides controlled pandas extraction from `obs`, `var`, `x`,
  `raw`, `obsm`, `varm`, `obsp`, `varp`, and tabular `uns` entries.
- Dense, sparse, DataFrame, view, and backed adapters project rows/columns before
  conversion. Cumulative `max_matrix_values=` budgets reject complete plans
  before the first read.
- Backed axis operations support independent in-memory selections while backed
  same-shape mutation fails with typed errors.
- A fixed ASV manifest, raw-sample evaluator, peak-RSS runner, and performance
  report enforce the three v0.3 performance families.
- Long exports preserve observation-major tidy ordering, grouped key/count
  paths avoid per-group Python scans, and filtering joins read only their keys.
- Public accessor methods carry explicit return annotations, documentation
  examples are executable, and package metadata is ready for trusted
  publishing.
- A deterministic 8-cell by 6-feature downstream fixture covers nullable and
  categorical metadata, duplicate observation names, CSR data/layers/raw,
  dense axis embeddings, CSR/CSC pairwise containers, H5AD and Zarr round
  trips, concatenation, and a Scanpy preprocessing handoff.

## Remaining Work Before 1.0

The remaining work is deliberately narrower than the completed milestone
backlog:

1. Add bounded chunked reductions for dense, sparse, and backed adapters
   ([#10](https://github.com/mdmanurung/annplyr/issues/10)). They
   must match eager known-answer fixtures, respect cumulative budgets, and
   demonstrate bounded peak RSS on a fixed runner.
2. Stabilize consumer-facing typing for the dynamically registered
   `AnnData.ap` namespace and the grouped/expression return types. Internal
   mypy success alone is not sufficient; a downstream type-check fixture must
   pass ([#11](https://github.com/mdmanurung/annplyr/issues/11)).
3. Adopt only the vetted security hardening from the closed template update,
   with traceable action pins, minimal permissions, scoped dependency updates,
   and a clean zizmor audit
   ([#12](https://github.com/mdmanurung/annplyr/issues/12)).

Each item should be implemented through a focused issue with exact fixtures,
baselines, and acceptance commands. New verbs are not a 1.0 prerequisite
unless they close a demonstrated workflow gap.

Issue [#13](https://github.com/mdmanurung/annplyr/issues/13) completed the
downstream integration gate without adding a runtime dependency. Its generated
fixture is documented in `tests/integration/README.md`; mandatory CI runs the
five exact invariance/consumer cases on Python 3.12, 3.13, and 3.14, with a
separate advisory prerelease-consumer lane.

## Milestone 1: Public Contract And Errors

Status: complete for v0.3.

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

Status: complete for the supported v0.3 selector and expression surface.

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

Status: complete for v0.3.

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

Status: complete for persistent v0.3 grouping. Additive regrouping flags are
not part of the frozen interface; callers explicitly create a new grouping
with `group_by()`.

Make grouping behavior deterministic and consistent across verbs.

Issue-sized tasks:

- Add `ungroup`, `group_vars`, `group_keys`, and `group_data`.
- Keep regrouping explicit rather than adding shallow `add` and `drop` flag
  combinations to `group_by`.
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

Status: projection, adapters, ownership, and budgets are complete for v0.3;
bounded chunked reductions are the focused 1.0 follow-up above.

Refactor storage and namespace handling without changing the accessor name.

Issue-sized tasks:

- Register the accessor with `anndata.register_anndata_namespace("ap")`.
- Add an `Axis` abstraction and positional indexer normalization.
- Add a source registry for `obs`, `var`, `X`, layers, `obsm`, `varm`, `obsp`,
  `varp`, and controlled `uns` access.
- Add dense/sparse/backed matrix adapters with column projection and explicit
  densification budgets.
- Add copy and materialization policy helpers.
- Add dtype-preserving assignment helpers.

Acceptance gates:

- Internal matrix predicates do not call whole-object `adata.to_df()` or
  `.toarray()` unless explicitly in an export path.
- Backed read-only verbs work; backed mutating verbs either materialize
  explicitly or raise a typed error.
- AnnData-aligned containers remain shape-consistent after every tested verb.

## Milestone 6: Tidyr And Plot Extraction

Status: complete for the safe v0.3 extraction and rectangling surface.

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

Status: repository-side infrastructure is complete; the milestone closes when
the v0.3.0 tag, GitHub release, trusted PyPI publication, and clean-install
smoke are verified.

The incompatible generated template refresh was closed without merge. Its
security-oriented pieces are tracked separately in
[#12](https://github.com/mdmanurung/annplyr/issues/12).

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

## Scoped Single-Cell Utilities

`annplyr` can include small AnnData utility helpers when they make tidy metadata
workflows safer without turning the package into a biological analysis toolkit.
In scope: sample metadata extraction and joins, feature presence diagnostics,
safe obs/var name edits, duplicate-name reports, and Scanpy-compatible palette
storage in `uns`.

Out of scope for core `annplyr`: biological QC metric wrappers, species-specific
gene registries, mitochondrial/ribosomal/hemoglobin/pathway scoring, MALAT1 or
cell-cycle scoring, QC plotting wrappers, read/write convenience wrappers,
Seurat/LIGER conversion helpers, marker discovery, and cluster annotation.

## Non-Goals Before 1.0

- Direct ggplot wrappers.
- Full lazy execution across arbitrary backends.
- Lossy joins that duplicate cells or genes without explicit user opt-in.
- Whole-matrix long exports by default on large sparse or backed objects.
- Biological QC doctrine or curated species/pathway gene-set maintenance.
- Reimplementing AnnData storage semantics outside AnnData-native subsetting.
