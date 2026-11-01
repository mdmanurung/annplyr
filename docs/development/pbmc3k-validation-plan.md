# PBMC3k validation and filtering conveniences

Execution is authorized by the user's instruction to execute the combined validation and API plan. Existing working-tree changes are preserved; the starting diff, status, and commit are recorded under `/tmp/annplyr-validation/baseline/`.

1. Establish the baseline: unit/integration tests, executable documentation, available environments, scripts, and public API inventory. Record failures before changes.
2. Add `filter_obs` and `filter_exp` as validated delegates to the existing filter engine. Metadata defaults to obs; expression columns require `x:`, `raw:`, or a named layer qualifier. Both subset observations only. Preserve grouping, ownership, projected reads, and budgets. Reach accessor, generated flat functions and chained methods, grouped methods, exports, stubs, typing inventory, docs, bundled skill, and tests.
3. Add real PBMC3k comparisons with independent pandas/NumPy/SciPy/Scanpy references. Check container alignment, source identity, serialization, downstream execution, and deterministic numerical results. Keep biological interpretation out of these technical checks.
4. Execute scripts, notebook cells, documentation, coverage, type/lint checks, benchmark correctness smoke, and distribution checks. Fix reproducible failures at their shared cause and rerun affected checks.
5. Write `pbmc3k-validation.md` with environment/dataset provenance, exact commands and outcomes, public API/script coverage, external workflow mapping, and unresolved gaps.

## Review and acceptance

The existing source resolver and filter planner already provide projection and source partitioning. Scope checks belong in one shared validation helper, before any matrix read; the new entry points must not duplicate filtering. Qualified expression references avoid guessing between metadata and genes. General `filter` remains available for mixed-source or feature-axis operations. Opaque expressions need an explicit container for expression filtering. Literal-only metadata predicates remain supported. Flat functions are generated from accessor methods; grouped methods need explicit delegates to preserve grouped evaluation.

Scanpy >=1.12 is the checkout's integration contract; the 1.10.x tutorial is a workflow reference. Use explicit parameters and seeds rather than claiming historical output parity. Processed PBMC3k has scaled X and log-normalized raw. Preserve each expression source, and never use scaled X as input counts.

No unrelated utility is justified until a documented workflow demonstrates a gap. Existing mutation, table extraction, sample metadata, naming, and palette helpers cover the first candidate workflows. `pipe` retains the callable's return, including None; demonstrate Scanpy's return-value API where available instead of changing generic pipe semantics.

Acceptance requires passing scope/error and independent filtering regression tests, aligned real-data comparisons for all new entry points, current generated contracts, and explicit outcomes for every discovered executable. Missing runtimes, skipped tests, and unavailable external resources stay unverified. A local run does not verify the full CI Python matrix or the full performance protocol.
