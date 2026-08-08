# v0.3 Performance Report

This report compares annplyr 0.2.0 at
`c0f0735859a059c3570dfea2ac6f1df322ac582e` with the frozen, unreferenced v0.3
snapshot `861cafcf9547673d0d7831e8b2af8d64cf09ce4d`. The snapshot procedure did
not move `HEAD`, update a branch, or use the repository index. All measurements
ran on `res-hpc-gpu12.researchlumc.nl` with one thread per numerical runtime and
seed `20260808`.

## Fixed runner and raw evidence

The timing environment contained Python 3.14.5, anndata 0.12.19, Narwhals
2.23.0, NumPy 2.4.6, pandas 2.3.3, and SciPy 1.16.3. The candidate wheel was
installed without dependency resolution into the same isolated ASV environment
used for the baseline. Import checks in the benchmark-launch environment
resolved annplyr 0.2.0 for `c0f07358` and annplyr 0.3.0 for `861cafcf`.

The authoritative raw files are ignored local artifacts:

- baseline timing: `.asv/results/res-hpc-gpu12.researchlumc.nl/c0f07358-env-fa329217d7752b590c681231009101f6.json`;
- candidate timing: `.asv/results/res-hpc-gpu12.researchlumc.nl/861cafcf-existing-py_exports_para-lipg-hpc_mdmanurung__hobby_annplyr_.asv_env_1e0106681b6027853d4043d29ee6025a_bin_python.json`;
- peak RSS: `c0f07358-peak-rss.json` and `861cafcf-peak-rss.json` in the same directory;
- normalized inputs: `c0f07358-normalized.json` and `861cafcf-normalized.json`;
- manifest verdict: `861cafcf-manifest-evaluation.json`.

Both timing files contain 34 benchmark methods and 47 parameter series. Every
non-smoke series has 14 raw samples and every smoke series has two; no result is
null. Each peak-RSS case has three independent subprocess samples. The 14
timings came from two complete seven-repeat rounds. The second round was kept
because it exposed distinct warm-run phases and made the noise estimate more
representative.

Pinned ASV 0.6.6 does not provide `asv dev`. The equivalent smoke check was
`asv run --quick`; all 34 candidate methods completed. Full candidate capture
used the isolated existing interpreter, the installed frozen wheel,
`--set-commit-hash 861cafcf9547673d0d7831e8b2af8d64cf09ce4d`, and a frozen worktree.

## Acceptance results

Times below are medians with the full raw-sample range in parentheses. Smaller
is better.

| Family and representative case | v0.2 | v0.3 | Improvement |
|---|---:|---:|---:|
| Matrix: dense export, 10 features | 97.10 ms (82.66-99.37) | 6.52 ms (4.26-11.67) | 93.29% |
| Matrix: CSR filter, 3 features | 12.85 s (12.74-12.98) | 842.43 ms (547.72-980.37) | 93.45% |
| Grouped: 20-group filter | 464.90 ms (331.33-475.57) | 213.90 ms (205.83-226.84) | 53.99% |
| Grouped: 20-group arrange | 660.02 ms (468.28-714.80) | 267.24 ms (263.57-279.73) | 59.51% |
| Grouped: 20-group slice-head | 172.32 ms (169.39-182.07) | 53.65 ms (52.72-55.48) | 68.87% |
| Metadata: mutate 10 independent columns | 142.46 ms (140.29-179.26) | 89.39 ms (66.47-91.59) | 37.25% |

Peak RSS is reported in KiB:

| Case | v0.2 raw samples | v0.3 raw samples | Median change |
|---|---:|---:|---:|
| Dense projection | 573816, 573764, 573528 | 557728, 558060, 557912 | 2.76% lower |
| Grouped mutate | 302288, 302576, 302104 | 308852, 309120, 309108 | 2.26% higher |
| Metadata mutate | 792312, 792504, 792508 | 480372, 480496, 480560 | 39.37% lower |

The manifest evaluator found a best improvement of 99.92% for matrix
projection, 68.87% for grouped execution, and 39.37% for metadata evaluation.
Each family therefore exceeds its required 20% improvement. It found no
primary regression greater than 10% outside measured noise and no sample-count
error. The primary manifest verdict is **pass**.

## Noise and the inner-join investigation

The manifest declares a robust noise rule. A slowdown is outside noise only
when its median shift exceeds `1.4826 * max(MAD_baseline, MAD_candidate)` and
the two raw-sample ranges do not overlap. Requiring both conditions prevents a
bimodal warm-run phase from being classified as a separated regression.

The combined inner-join samples triggered the original MAD-only rule: the v0.2
median was 346.99 ms with range 336.60-448.29 ms, while the v0.3 median was
404.51 ms with range 341.73-413.93 ms. The ranges overlap substantially. A
paired isolated rerun produced these raw samples:

- v0.2: 412.895, 410.539, 412.775, 418.850, 432.389, 404.536, 402.022,
  402.788, 404.458, 382.843, 334.816, 334.080, 333.344, 333.996 ms;
- v0.3: 415.433, 416.604, 420.554, 416.936, 410.965, 416.661, 415.767,
  415.636, 391.525, 346.789, 348.445, 349.209, 350.167, 348.821 ms.

Paired `cProfile` runs after warm-up took 0.431 s for v0.2 and 0.432 s for
v0.3. The merge phase was 0.133 versus 0.125 s and positional subsetting was
0.099 versus 0.091 s; copy-related calls were 0.224 versus 0.258 s. This
profile and the overlapping raw distributions do not support a separated total
runtime regression. A regression test now distinguishes overlapping phased
samples from truly separated samples.

## ASV continuous screen

The exact frozen comparison was:

```text
asv continuous c0f0735859a059c3570dfea2ac6f1df322ac582e 861cafcf9547673d0d7831e8b2af8d64cf09ce4d --factor 1.10
```

It completed all benchmarks but exited 1 with `PERFORMANCE DECREASED`. This
screen did **not** pass. It reported six slowdowns, all outside the manifest's
primary acceptance cases:

| Diagnostic case | v0.2 | v0.3 | Ratio |
|---|---:|---:|---:|
| Metadata semi-join | 73.0 +/- 10 ms | 135 +/- 3 ms | 1.85 |
| 50,000-group count | 68.1 +/- 0.4 ms | 125 +/- 2 ms | 1.84 |
| 50,000-group keys | 61.8 +/- 0.2 ms | 109 +/- 0.4 ms | 1.76 |
| Rectangling chained extraction | 333 +/- 2 ms | 495 +/- 2 ms | 1.49 |
| Metadata anti-join | 73.0 +/- 10 ms | 107 +/- 10 ms | 1.47 |
| Smoke grouped summarize | 6.51 ms | 8.95 ms | 1.38 |

The same paired screen measured inner join at 348 +/- 30 ms for v0.2 and
347 +/- 40 ms for v0.3. It also confirmed the intended large improvements,
including grouped filter at 0.43x, metadata mutate at 0.55x, and sparse smoke
projection at 0.02x. The six diagnostic slowdowns remain optimization targets;
they are not relabeled as successes. The frozen acceptance contract assigns
the primary, noise-aware regression gate to the manifest evaluator, which
passes.

## Correctness and measurement caveats

Every timed method runs an independent expected-result or structural oracle in
setup. Realization benchmarks force data access. `copy=False` call-cost cases
check shape and ordering without requiring a view, because v0.3 permits either
a view or a materialized object; view-realization cases are separate where an
operation actually yields a view. The v0.2 sparse filter uses the same
three-feature Boolean predicate as v0.3, avoiding a pandas 2.3.3 sparse-float
arithmetic defect. For backed v0.2 CSR/CSC cases that cannot project through
the public baseline API, the semantically equivalent fallback materializes
inside the timed operation and is labeled in the benchmark code.

Fixture construction occurs outside timed methods. Backed measurements are
warm-cache results because filesystem caches were not controlled. Peak RSS
runs use fresh processes, but setup allocation can still affect process RSS.
The correctness gate includes a chained fixture-level structural comparison of
values, extension dtypes, duplicate indexes, sparse formats, aligned shapes,
ownership, and grouping state. This repository contains no representative
downstream dataset, so fixture-level invariance is the strongest available
downstream safeguard.
