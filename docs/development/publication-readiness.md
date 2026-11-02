# Publication readiness of annplyr 0.5.0

Audit run on 2026-11-01 on branch `release/0.5.0-readiness`. The target is a
JOSS paper and a listing in the scverse ecosystem.

## Verdict

The code is **ready, with caveats**: every check that can run locally passes.
It cannot be cited or installed yet, because there is no release, no PyPI
upload and no DOI. Only the maintainer can do those steps (see
"Maintainer gaps" below).

## Blockers found

Each blocker is fixed on this branch unless marked open.

| Blocker | Evidence | Status |
|---|---|---|
| Docs used the pre-0.5.0 API | `docs:doctest`: 36 of 69 failed, e.g. `filter() got an unexpected keyword argument 'obs'`. README and the 3 notebooks had the same problem. | Fixed: 78 of 78 pass |
| Notebooks never ran at build | `nb_execution_mode` was `"off"` in `docs/conf.py` | Fixed: `cache` mode, and errors fail the build |
| ASV benchmarks used the old API | `bench_smoke.py:33` called `group_by(obs=...)`, so the `benchmark-smoke` CI job failed | Fixed: all 48 benchmark methods run |
| anndata 0.13 broke `TidyAnnData` | `test_chaining.py`: `'super' object has no attribute 'concatenate'` | Fixed: the override only exists when anndata has the method |
| Categorical join key demoted to `str` | After `left_join`, `uns["louvain_colors"]` no longer matched the categories | Fixed, with a regression test |
| `arrange` ignored key order | `arrange("louvain", desc("n_genes"))` sorted by `n_genes` first | Fixed, with a regression test |
| `arrange` sorted categoricals alphabetically | dplyr and pandas sort by level order | Fixed, with a regression test |
| `write_zarr` failed on every verb result | anndata `IORegistryError` for the subclass | Fixed, with a regression test |
| Declared floors were untested | `narwhals 2.14` lacks `Expr.__neg__`. `anndata<0.12.4` fails backed fancy indexing. | Fixed: floors raised, and a `lowest-direct` env runs in CI |
| Not on PyPI, no tag, no DOI | `pypi.org/pypi/annplyr/json` returns 404. `git tag` is empty. | **Open**: maintainer step |

## Results per check

Unit tests and their environments:

| Check | Command | Outcome |
|---|---|---|
| Unit tests, 5 envs | `uvx hatch test --all -n auto` | Pass: 392 passed and 2 skipped in each of py3.12/3.13/3.14 stable, 3.14 pre and 3.12 min |
| Integration | `uvx hatch run integration:run` | Pass: 64 of 64, including `aggregate_cells` against `sc.get.aggregate` |
| Coverage | `uvx hatch test -py 3.12 -i deps=stable --cover`, then `cov-combine` and `cov-report` | Pass: 87% against `fail_under = 85`; `_verbs.py` is lowest at 80% |
| Clean wheel install | Fresh py3.12 venv, wheel plus pytest and scanpy, tests run against site-packages | Pass: 371 of 371, before later tests were added |
| Lowest dependency versions | `uvx hatch test -py 3.12 -i deps=min` | Pass: anndata 0.12.4, narwhals 2.22.0, numpy 2.0.0, pandas 2.2.2, scipy 1.13.0 |
| numpy 1.26 and scipy 1.12 floors | Forced into the env by hand | Pass, but checked manually only; CI resolves numpy 2.0 |
| Test-order independence | `pytest-randomly` seeds 11 and 222 | Pass: no order-dependent failures |

Types, lint and generated files:

| Check | Command | Outcome |
|---|---|---|
| Types | `uvx hatch run type:check` and `type:consumer` | Pass, after fixing 5 mypy errors |
| Lint and format | `uvx prek run --all-files` | Pass, except two hooks below |
| `python` entry hooks | `render-typing-contract`, `render-tidy-stub` | **Unverified as hooks**: this host has no `python` on PATH. Both scripts pass `--check` when run from the hatch env. |
| Generated files | `scripts/render_typing_contract.py --check` and `scripts/render_tidy_stub.py --check` | Pass |

Docs, build and benchmarks:

| Check | Command | Outcome |
|---|---|---|
| Docs build | `uvx hatch run docs:build`, which runs with `-W` | Pass: 0 warnings, notebooks executed |
| Doctests | `uvx hatch run docs:doctest` | Pass: 78 of 78, on a fresh `-E` build |
| Build | `uvx hatch build -c`, then `uvx twine check dist/annplyr-0.5.0*` | Pass |
| Benchmark smoke, run directly | Every `time_*` and `peakmem_*` method run once in the hatch env | Pass: 48 of 48 |
| Benchmark smoke through ASV | `asv_local.py run --quick -E existing -b "^bench_smoke"`, in a scratch copy with the CI benchmark stack | Passed 5 of 5 before the API redesign. **Not rerun** since; the direct run above covers the same methods. |
| CITATION.cff | `uvx cffconvert --validate -i CITATION.cff` | Pass: schema 1.2.0 |
| JOSS PDF | `.github/workflows/draft-pdf.yml` | **Unverified**: runs in GitHub Actions only |
| Workflow audit | `uvx zizmor@1.29.0 --offline --persona pedantic .github/workflows/` | Pass: 0 findings, after adding `concurrency` to `draft-pdf.yml`. The action is now in `ci-security.md`. |
| Paper example output | Diff of the `paper.md` table against the doctested `docs/index.md` output | Pass: identical code and output |
| Docs CI in Actions | `docs.yaml` now runs the notebooks and downloads PBMC3k (about 24 MB) | **Unverified**: run locally only. It needs network, and adds a few minutes. |
| Error messages name the fix | Raise sites sampled, and the typed errors asserted in the PBMC3k tests | **Partial**: unknown prefixes, sources and a non-dict `mutate` name the fix. Some planner and selection messages state only the problem. |

## Scientific validation

`tests/integration/test_pbmc3k_reference.py` compares annplyr against an
independent computation on Scanpy's processed PBMC3k. The references are
pandas, NumPy and `sc.get.obs_df`. That data has a scaled dense `X` over 1,838
genes and a log-normalized sparse `.raw` over 13,714 genes. The tests cover:

- values for every verb family, grouped and ungrouped, at `rtol` 0 for
  selection, 1e-6 for float32 arithmetic and 1e-5 for aggregates;
- `x:` and `raw:` over dense and CSR at float32 and float64, with `raw:`
  proven not to read `X`;
- alignment of `X`, `.raw`, `obsm`, `varm`, `obsp`, `layers` and `uns` after
  each verb, and identical results on shuffled input once aligned by name;
- typed failures for many-to-many joins, unknown columns, sources and genes,
  and duplicate names; NA for unmatched keys; views never write to the parent;
- h5ad and zarr round trips, and repeat runs giving identical output.

These are technical checks of wrangling semantics. They make no biological
claim.

## Changes made

| Change | Guarded by |
|---|---|
| `TidyAnnData` supports anndata 0.13 | `tests/test_chaining.py::test_anndata_methods_that_build_new_objects_stay_chainable` |
| Typing contract renders the inherited constructor as inherited | `tests/test_consumer_typing.py::test_generated_typing_contract_is_current` |
| `filter_obs` and `filter_exp` docstrings carry the full contract | `tests/test_public_docs_v03.py::test_accessor_and_grouped_docstrings_cover_contract_sections` |
| Joins keep a categorical key's dtype and order | `tests/test_joins_tidyr_scverse.py::test_join_and_count_keep_categorical_key_dtype_and_order` |
| `arrange` applies keys in written order | `tests/test_tidyverse_parity.py::test_arrange_applies_keys_left_to_right_with_bare_names` |
| `arrange` sorts a categorical by level order | `tests/test_tidyverse_parity.py::test_arrange_sorts_a_categorical_by_category_order` |
| Unknown sort key raises `UnknownColumnError` | `tests/test_errors.py::test_unknown_sort_key_raises_unknown_column_error` |
| Repeated name in `select` raises `DuplicateNameError` | `tests/test_errors.py::test_select_with_a_repeated_name_raises_duplicate_name_error` |
| `mutate` without a dict names the fix | `tests/test_errors.py::test_mutate_without_a_dict_names_the_fix` |
| Verb results write to zarr | `tests/integration/test_downstream_fixture.py::test_verb_results_can_be_written_to_zarr` |
| Docs, README, notebooks and benchmarks use the 0.5.0 API | `docs:doctest`, the notebook-executing `docs:build` and the `benchmark-smoke` job |
| Dependency floors raised | the `hatch-test.py3.12-min` env |
| Pre-commit `check-ast` runs on Python 3.12 | `prek run --all-files` |
| Trailing-whitespace hook skips docs Markdown | `docs:doctest`. The hook had broken 14 doctests. |
| JOSS `paper.md`/`paper.bib`, `.zenodo.json`, `CITATION.cff`, README citation section | `cffconvert --validate` |
| Donor-level tutorial section | `docs:doctest` on `docs/vignettes/cohort_workflow.md` |

## API redesign for R parity

After the audit, the 0.5.0 API was reworked to match tidyseurat,
tidySingleCellExperiment and tidytof as closely as Python allows. The full
side-by-side table is in `docs/user_guide/coming_from_r.md`, and the migration
table is in the CHANGELOG.

- Container namespaces: `obs.percent_mito`, `raw.MS4A1`, `counts.MS4A1`,
  `raw["HLA-DRA"]`, from `annplyr.data`.
- Positional names in every verb, keyword assignments in `mutate` and
  `summarize`, and `by=` on the slice verbs.
- New verbs: `join_features`, `aggregate_cells`, `nest`/`unnest`, `glimpse`, and
  `ap.tidy`.
- Removed: duplicate verbs, `copy=`/`inplace=`, and per-container keywords. The
  pandas tidyr helpers moved to `annplyr.tidyr`.

The API review found and fixed these bugs, each with a regression test:
- `ap.layer("obs")` read obs metadata instead of the layer;
- the `obsm:` prefix used the wrong key;
- grouped `summarize(by=)` was ignored;
- `annplyr.data` minted fake names for Sphinx, which broke a fresh docs build.

A fresh agent that had read only the docs did ten tidySCE vignette tasks on
PBMC3k, and all ten worked on the first try.

R parity was then measured with `scripts/r_parity/`. The same PBMC3k matrix
(`.raw`, with synthetic donors) went to annplyr and to tidySingleCellExperiment
(R 4.5.2, Bioconductor 3.22). The tasks were filter, count, mutate, grouped
summarise, `slice_max(by=)`, `join_features` long, rename with `if_else`,
`left_join`, `aggregate_cells` and `arrange`. All ten match at `rtol` 1e-6:
- the same cells, in the same order, for filter, slice and arrange;
- the same counts;
- `aggregate_cells` sums within 4.4e-15 relative.

Two details:
- **Fixed during this check.** `aggregate_cells` summed float32 data in
  float32, which drifted by 1.2e-6; it now accumulates floats in float64.
- **Still different.** `count` and grouped `summarize` return groups in
  first-seen order, where dplyr sorts them by key; the values are identical.

anndataR did not build here, because its `rhdf5filters` dependency fails to
compile, so the scripts pass the matrix to R as Matrix Market + CSV.

## Maintainer gaps

The maintainer has to supply these, and the files hold `[gap: ...]` markers
where they belong.

- **ORCID and affiliation**: `CITATION.cff`, `paper.md` front matter, and
  `.zenodo.json` creators. JSON has no comments, so add the keys there.
- **Funding and acknowledgements**: `paper.md`.
- **Release date**: `CHANGELOG.md`, `CITATION.cff` and `tests/test_skills.py`
  say 2026-09-27, but nothing was released then. Set the real date when tagging.
- **Release steps, in order**:
  1. Merge this branch.
  2. Turn on the Zenodo GitHub integration.
  3. On pypi.org, register a pending trusted publisher for project `annplyr`.
     Use workflow `release.yaml` and environment `pypi`. Create the `pypi`
     environment in the GitHub repo settings.
  4. Publish GitHub release `v0.5.0`. `release.yaml` then builds and uploads
     to PyPI through trusted publishing, and Zenodo mints the DOI.
  5. Put the concept DOI into `CITATION.cff` and the README badge.
- **JOSS**: submit `paper.md`. Once JOSS accepts it, uncomment and fill
  `preferred-citation` in `CITATION.cff`. It stays commented until then, so
  GitHub cites the software rather than an article that does not exist.
- **scverse**: open a `meta.yaml` PR to `scverse/ecosystem-packages`. That
  needs the PyPI release first, plus your consent to the listing and to the
  scverse code of conduct.

## Known limits, not fixed

- `slice_sample` defaults to `random_state=None`, so it is unseeded unless the
  caller passes a seed. This matches pandas and dplyr.
- After `filter`, `semi_join` and `anti_join`, anndata drops unused categories
  from the subset. `left_join` and `inner_join` keep the full list.
- `ad.io.write_elem` on a `TidyAnnData` into zarr is untested and probably
  needs `to_anndata()` first.
- Backed-mode round trips and the full performance protocol were not run.

## Environment and provenance

- Host: Linux 4.18.0 (el8).
- Commits audited: start at `7303d49`; the readiness fixes come after it on
  `release/0.5.0-readiness`.
- Integration env: Python 3.14.5, anndata 0.13.4, narwhals 2.28.0,
  pandas 3.1.0, numpy 2.5.4, scipy 1.18.1, scanpy 1.12.5, zarr 3.4.2,
  h5py 3.16.0.
- Unit tests ran on Python 3.12, 3.13 and 3.14 in the stable envs, on 3.14 with
  prereleases, and on 3.12 at the lowest direct versions listed above.
- Dataset: `scanpy.datasets.pbmc3k_processed()`, cached at
  `data/pbmc3k_processed.h5ad`. It is the 10x Genomics PBMC3k data, processed
  as in the Scanpy tutorial. PBMC3k is a single donor, so the donor labels in
  the docs are synthetic and marked as such.
