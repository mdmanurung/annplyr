# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.5.0] - 2026-09-27

annplyr 0.5.0 reads like tidyseurat and tidySingleCellExperiment. Every column
name carries its container (`obs.percent_mito`, `raw.MS4A1`), verbs take names
positionally, and new columns are keywords. The "Coming from R" page of the user guide
has the side-by-side table with the R packages.

### Migration

| 0.4 | 0.5 |
|---|---|
| `adata.ap.filter(...)` | `adata = ap.tidy(adata)`, then `adata.filter(...)`; `adata.ap.filter(...)` still works |
| `ap.col("percent_mito")` | `obs.percent_mito`, after `from annplyr.data import obs` |
| `ap.col("raw:MS4A1")`, `ap.col("x:MS4A1")` | `raw.MS4A1`, `X.MS4A1` |
| `ap.col("counts:MS4A1")`, `ap.col("X_pca:0")` | `counts.MS4A1`, `X_pca[0]`; `ap.layer("name")` for a reserved layer name |
| `ap.col("raw:HLA-DRA")` | `raw["HLA-DRA"]` |
| `select(obs=["a", "b"], var=["c"], x=["G1"])` | `select(obs.a, obs.b, var.c, X.G1)` |
| none | `select(-obs.a)` drops a column |
| `rename(obs={"cell_type": "louvain"})` | `rename(cell_type=obs.louvain)` |
| `relocate(obs=["a"], before="b")` | `relocate(obs.a, before="b")` |
| `distinct(obs=["a"])`, `pull(x="G1")` | `distinct(obs.a)`, `pull(X.G1)` |
| `count(by=["a", "b"])`, `add_count(by="a")` | `count(obs.a, obs.b)`, `add_count(obs.a)` |
| `group_by("a", "b")` | `group_by(obs.a, obs.b)` |
| `mutate({"b_score": expr})` | `mutate(b_score=expr)`; a dict stays for names that are not identifiers |
| `summarize({"n": ap.n()}, by="louvain")` | `summarize(n=ap.n(), by=obs.louvain)` |
| `group_by(...).summarize({...}, by=...)` | `group_by(...).summarize(...)`; `by=` raises `TypeError` |
| `to_wide(obs=["a"], raw=["G1"], obsm={"X_umap": [...]})` | `to_wide(obs.a, raw.G1, X_umap)` |
| `to_tidy(obs=["a"], raw=["G1", "G2"])` | `join_features(raw.G1, raw.G2, shape="long")` with `.cell`, `.feature`, `.abundance_raw` |
| `as_frame("varm", key="PCs", select=["0"])` | `pull(PCs[0])` |
| `as_frame(...)` for a table | `to_wide(...)` |
| `pivot_longer(obs=[...], x=[...])` | `pivot_longer(obs.a, X.G1, names_to=..., values_to=...)` |
| `nest_by(by="a")` | `nest(obs.a)`, one object per group; `ap.unnest(table)` reverses it |
| `filter_obs(...)`, `filter_exp(...)` | `filter(...)` |
| `tally()`, `add_tally()` | `summarize(n=ap.n())`, `add_count()` |
| `sample_summary(...)` | `summarize(..., by=obs.sample)` |
| `add_sample_meta(...)` | `left_join(table, by="sample")` |
| `copy=True`, `copy=False`, `inplace=True` | removed; every verb returns a new object, so assign the result |
| `ap.pivot_wider`, `ap.nest`, `ap.unnest`, `ap.separate_rows`, and 12 more | `ap.tidyr.pivot_wider` and so on; `ap.nest` and `ap.unnest` now work on AnnData |
| `filter(obs=...)`, `mutate(obs=..., x=...)` and other per-container keywords (0.4 and earlier) | one flat payload; the expression names its container |
| `adata.ap.to_df(...)` | `adata.to_wide(...)` |
| requiring a plain `AnnData` downstream | `result.to_anndata()` |

### Changed

- `join_features(shape="long")` and `pivot_longer` return dense values in the
  source dtype, as tidySingleCellExperiment does, so plotnine and seaborn draw them
  directly. `to_wide` keeps pandas sparse columns for wide exports of many genes.

- **Breaking:** `to_wide` names columns as tidySingleCellExperiment, Seurat and
  ggann do. Genes keep their own names (`MS4A1`, not `raw_MS4A1`), and embedding
  columns are 1-based with a short prefix (`X_umap` gives `UMAP_1`, `UMAP_2`;
  `X_pca` gives `PC_1`). A gene named like an `obs` column raises
  `DuplicateNameError`.

- **Breaking:** every column name carries its container. `from annplyr.data
  import obs, var, X, raw` gives the fixed containers, and any other name is a
  layer or an `obsm`/`varm` key (`from annplyr.data import counts, X_umap`).
  `obs.name` and `raw["HLA-DRA"]` build the same expression as the old
  `ap.col("obs:name")` and `ap.col("raw:HLA-DRA")`. Strings and `ap.col` still
  work for names computed in a loop. A mistyped layer fails at the verb with
  `UnknownSourceError`, which lists the layers and keys that exist. Calling a
  container (`counts("x")`) raises a `TypeError` that names the fix.
- **Breaking:** a bare name is `obs` metadata and a matrix always needs its
  container, because there is no safe default. One expression cannot mix
  containers; it fails with `IncompatibleAxisError` naming both. A raw Narwhals
  expression is evaluated against `obs` unless the verb's `source=` names a
  container, which is charged for in full.
- **Breaking:** `select`, `rename`, `rename_with`, `relocate`, `distinct`,
  `pull`, `to_wide`, `pivot_longer`, `count`, `add_count` and `group_by` take
  names positionally: column expressions, selectors, or a whole container. Every
  per-container keyword (`obs=`, `var=`, `x=`, `raw=`, `layer=`, `obsm=`,
  `obsp=`) is gone. A whole `obs` or `obsm` key means all of its columns; a
  whole `X`, `raw` or layer in `to_wide` and `pivot_longer` needs
  `allow_all_features=True`.
- **Breaking:** `mutate`, `transmute` and `summarize` take new columns as
  keywords, like dplyr. A positional dict is still accepted for names that are
  not identifiers. A name given both ways raises `DuplicateNameError`. The flat
  functions take `adata` positional-only so the dict can be positional-only too.
  `mutate` and `transmute` write one destination per call, chosen by
  `into="obs"|"var"`.
- **Breaking:** `summarize` (with `by=` or on a grouped object) and `count` list
  groups in key order, as dplyr does, where they listed them in first-seen order.
  A categorical key follows its level order, missing keys come last, and
  `count(sort=True)` keeps key order within ties. `group_keys`, `nest`,
  `aggregate_cells` and the AnnData-returning verbs keep first-seen order.
- **Breaking:** `group_by`, `count`, `nest` and `aggregate_cells` accept a column
  expression with a single source column as a key. A computed key such as
  `obs.a + 1` raises `SelectionError` that names `mutate`.
- **Breaking:** every verb returns a new object, so `copy=` and `inplace=` are
  removed. Subsetting verbs return an independent object, except that a backed
  input returns the view, so a `.h5ad` can be subset without loading it.
  Same-shape verbs (`mutate`, `transmute`, `rename`, `rename_with`, `relocate`,
  `add_count` and the name utilities) copy `obs`, `var` and `uns` and share `X`,
  layers and the aligned matrices with the input. `select` no longer rewrites
  `obs` or `var` when their columns do not change.
- **Breaking:** removed `filter_obs`, `filter_exp`, `tally`, `add_tally`,
  `to_tidy`, `as_frame`, `nest_by`, `sample_summary` and `add_sample_meta`.
  `join_features`, `to_wide`, `pull`, `summarize`, `left_join` and `nest` cover
  them (see the table).
- **Breaking:** the 16 pandas-only tidyr helpers (`pivot_wider`, `nest`,
  `unnest`, `chop`, `unchop`, `pack`, `unpack`, `separate`, `separate_rows`,
  `extract`, `unite`, `drop_na`, `fill`, `hoist`, `unnest_longer`,
  `unnest_wider`) moved to `annplyr.tidyr`. They take DataFrames, not AnnData.
- **Breaking:** `to_df()` is renamed `to_wide()`. `AnnData.to_df(layer=...)`
  already exists upstream and dumps the whole of `X`; annplyr's verb takes
  selectors and was shadowing it on results.
- A grouped `summarize` takes no `by=`. It used to accept the argument and ignore
  it, which gave a plausible wrong table.
- `AnnDataWithAnnplyr` is now an alias of the concrete `TidyAnnData` class
  rather than a typing-only stub aliasing `AnnData`. It is unchanged for
  annotations; at runtime it now identifies the class results actually are.
- `GroupedAnnData.ungroup()` and grouped iteration return wrapped results, so a
  chain continues past the group boundary. `group_by()` with no arguments now
  returns a wrapped result like every other AnnData-returning verb.
- `__all__` is now 119 symbols, from 97 in 0.4.

### Added

- `slice_min`, `slice_max`, `slice_head`, `slice_tail` and `slice_sample` take
  `by=`, as in dplyr 1.1: `slice_max(obs.n_counts, n=3, by=obs.louvain)`. Their
  first argument is now called `order_by`.
- `X[ap.any_of(panel)]` and `X[ap.starts_with("CD")]` pick features with a
  selector inside a container.
- `ap.tidy(adata)` returns a `TidyAnnData`, a subclass of `AnnData` that carries
  every verb as a method, so `adata.filter(...).select(...)` chains without
  `.ap`. `isinstance(result, AnnData)` still holds, so Scanpy and AnnData are
  unaffected, and `result.to_anndata()` returns a plain `AnnData`. A
  user-defined `AnnData` subclass keeps its methods and `isinstance`
  relationships. `adata.ap.verb(...)` and `ap.verb(adata, ...)` call the same
  implementation.
- `annplyr.data`, `annplyr.Source` and `ap.layer(name)`. `layer` reads a layer
  whose name is reserved or not importable, such as one called `obs`.
- `join_features(*features, shape="long"|"wide")`, after tidySCE. `"long"`
  returns a DataFrame of every `obs` column plus `.cell`, `.feature` and one
  `.abundance_<container>` column per container. `"wide"` returns the object
  with each feature added to `obs`.
- `aggregate_cells(*by, source=)`, a pseudobulk AnnData with one row per observed
  key combination and `n_cells` in `obs`. `source` has no default because `X` is
  often scaled. `fun="mean"` averages instead of summing.
- `nest(*by)`, a DataFrame with one row per group and that group's object in a
  `data` column, and `ap.unnest(table)`, which stacks the objects back.
- `glimpse()` prints the `obs` head and the layer, `obsm` and `raw` names.
- `select(-obs.x)`, `relocate(-obs.x)` and `distinct(-obs.x)` read a negated bare
  column as "everything except". Elsewhere `-obs.x` is numeric negation.
- `src/annplyr/_tidy.pyi`, a generated stub that gives static checkers the same
  verb signatures as `inspect.signature`. Regenerate with
  `python scripts/render_tidy_stub.py --write`. `annplyr/data.pyi` lets mypy
  accept `from annplyr.data import counts`.
- The user guide page "Coming from R" maps each verb to its tidyseurat,
  tidySCE and tidytof spelling and lists where parity is impossible.

### Fixed

- The `obsm:` prefix resolved to the key `"obsm"` instead of the key written
  after it.
- `ap.layer("obs").MS4A1` read the `obs` column `MS4A1`, a wrong number with no
  error. It reads the layer, and checks that the layer exists.
- The `max_matrix_values` error from `join_features` names `join_features`, not
  the removed `to_tidy`.
- `left_join` and `inner_join` keep a categorical key's dtype and category
  order, so `uns["<key>_colors"]` stays aligned. `count` and `summarize(by=...)`
  keys keep categorical and nullable dtypes too.
- `arrange` applies keys in the order written. Bare names used to sort after
  every expression key.
- `arrange` sorts a categorical column by its category order, as dplyr does,
  not alphabetically.
- An unknown sort key in `arrange` or `slice_min`/`slice_max` raises
  `UnknownColumnError` instead of a raw Narwhals error.
- `select` with a repeated name raises `DuplicateNameError`.
- `mutate` and `transmute` given something other than a dict raise
  `SelectionError` naming the fix, instead of `AttributeError`.
- Verb results write to zarr. anndata's IO registry rejected the
  `TidyAnnData` subclass.
- `TidyAnnData` works with anndata 0.13, which removed `AnnData.concatenate`.
- Every docs page, notebook, the README, the skill references and the ASV
  benchmarks use the 0.5.0 API. Notebooks execute during the docs build, so a
  stale one fails it.

### Dependencies

- Floors raised to the lowest tested versions: `anndata>=0.12.4`,
  `narwhals>=2.22`, `numpy>=1.26`, `pandas>=2.2.2`, `scipy>=1.12`. A
  `lowest-direct` test environment checks them in CI.

## [0.4.0] - 2026-08-10

### Fixed

- `summarize()` no longer attaches group keys to the wrong aggregates. When a
  request combined a metadata source with a matrix source, the pieces were
  joined with an outer merge, which sorts its join keys, and the key column was
  then rewritten in first-seen order. Every group whose first-seen position
  differed from its sorted position reported another group's numbers. Summary
  sources are now keyed by an internal group id and keys are attached per row,
  which makes the misalignment unrepresentable.
- Expressions that combine two matrix columns now work for every sparse
  subtype. pandas ships sparse binary kernels for `float64`, `int64`, and
  `bool` only, so a two-gene score over a `float32` `X`, layer, or `raw`, the
  usual single-cell case, previously failed with an opaque
  `sparse_add_float32` error in `mutate()`, `filter()`, `arrange()`, and
  `summarize()`. Evaluation now densifies only the affected projected columns,
  so results match the dense path exactly and exported frames keep their
  sparse columns and subtype.

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
- Every executable example runs against Scanpy's PBMC3K instead of a four-cell
  synthetic fixture, and prints the values it produces. The quickstart, user
  guide, and cohort vignette contrast each annplyr call with the scanpy and
  pandas assembly it replaces, and both sides execute during the documentation
  build. The three notebooks are re-authored on the same dataset and committed
  with executed outputs.
- `scanpy` joins the `doc` dependency group, `docs/_pbmc.py` caches the shared
  fixture, and the Docs workflow caches the downloaded dataset.

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
