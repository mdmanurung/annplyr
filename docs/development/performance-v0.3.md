# v0.3 Performance Report

This report compares annplyr 0.2.0 at
`c0f0735859a059c3570dfea2ac6f1df322ac582e` with the frozen, unreferenced v0.3
snapshot `861cafcf9547673d0d7831e8b2af8d64cf09ce4d`. The snapshot procedure did
not move `HEAD`, update a branch, or use the repository index. All
frozen-snapshot measurements ran on `res-hpc-gpu12.researchlumc.nl` with one
thread per numerical runtime and seed `20260808`.

## Fixed runner and raw evidence

The timing environment contained Python 3.14.5, anndata 0.12.19, Narwhals
2.23.0, NumPy 2.4.6, pandas 2.3.3, and SciPy 1.16.3. The candidate wheel was
installed without dependency resolution into the same isolated ASV environment
used for the baseline. Import checks in the benchmark-launch environment
resolved annplyr 0.2.0 for `c0f07358` and annplyr 0.3.0 for `861cafcf`.

The authoritative frozen-snapshot raw files are ignored local artifacts:

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
projection at 0.02x. At the time of the frozen snapshot, the six diagnostic
slowdowns remained optimization targets; they are not relabeled as successes
in this historical screen. The frozen acceptance contract assigns the primary,
noise-aware regression gate to the manifest evaluator, which passes. The
release-hardening measurements below reassess the named regressions without
rewriting this frozen evidence.

## Post-freeze release hardening

The frozen screen drove bounded fixes for high-cardinality grouping, filtering
joins, and long extraction. The pre-fix source was merged-main commit
`ef136845cc425eb157530247a626eac4d23a2219`; the release-candidate source was
commit `a6edfb2c05be3f865fa46d77b2be22db454fec2d`. Direct same-host
comparisons used the pinned
Python 3.14.5 ASV interpreter, seed, one numerical thread, warmed fixtures, and
disabled garbage collection during samples. These targeted measurements are
separate from, and do not replace, the commit-level ASV continuous screen.
Times are medians in milliseconds with full sample ranges.

| Diagnostic case | v0.2 | Pre-fix v0.3 | Release candidate | Change from pre-fix |
|---|---:|---:|---:|---:|
| 50,000-group keys | 80.457 (79.634-83.200) | 128.027 (125.325-131.691) | 2.570 (2.533-3.144) | 97.99% lower |
| 50,000-group count | 107.015 (105.639-111.486) | 150.680 (147.928-163.993) | 9.473 (8.430-10.709) | 93.71% lower |
| Metadata semi-join | 121.249 (120.171-122.949) | 189.276 (187.415-193.312) | 115.059 (104.428-118.742) | 39.21% lower |
| Metadata anti-join | 130.929 (126.960-167.317) | 190.017 (186.640-198.171) | 118.308 (113.481-131.224) | 37.74% lower |
| Dense `pivot_longer` | 124.301 (110.650-132.468) | 113.191 (106.781-122.435) | 61.753 (56.149-89.533) | 45.44% lower |
| Chained extraction | 876.682 (843.969-929.550) | 904.878 (821.129-930.365) | 697.213 (680.012-721.597) | 22.95% lower |
| CSR `to_tidy` | 20,201.242 (20,054.017-20,563.892) | 80.406 (75.215-81.589) | 70.055 (67.075-72.315) | 12.87% lower |
| CSC `to_tidy` | 20,496.622 (20,479.783-21,139.518) | 94.671 (90.327-97.267) | 86.847 (83.073-93.042) | 8.26% lower |

The release candidate has a lower median than both the pre-fix v0.3 code and
v0.2 in all eight targeted cases. The release-candidate and v0.2 anti-join
ranges overlap, so that comparison is treated as noisy rather than as a
separated improvement. Group and join measurements used 11 repeats; chained
extraction used seven. The v0.2 sparse cases used three repeats because each
sample took about 20 seconds.

Result gates accompanied every comparison. Dense long output was exactly
observation-major with shape `(1,000,000, 5)`, matching dtypes, and row hash
`7426557883411314069`. CSR and CSC outputs both had shape `(1,000,000, 3)`,
preserved `Sparse[float32, 0.0]`, and had hashes `15146737095588008042` and
`3289997898914792372`. The chained result had shape `(20,000, 53)`, 52
`float32` columns plus one object column, and hash `10322975911543133315`.
Independent grouping and join oracles additionally covered categorical and
mixed-null keys, multiple keys, duplicate indexes, weighted counts, both axes,
and both null-matching policies.

### Release-candidate commit-level rerun

A complete ASV 0.6.6 comparison on
`res-hpc-exe050.researchlumc.nl` then compared the v0.2 commit with
`a6edfb2c05be3f865fa46d77b2be22db454fec2d`:

```text
asv continuous c0f0735859a059c3570dfea2ac6f1df322ac582e \
  a6edfb2c05be3f865fa46d77b2be22db454fec2d \
  --factor 1.10 --show-stderr
```

It completed all 34 benchmark methods and 47 parameter series over two rounds,
exited zero, reported `PERFORMANCE INCREASED`, and reported no statistically
classified decreases. The same raw timing files and fresh-process RSS samples
also passed the normalized manifest: zero sample errors, zero primary
regressions outside noise, and best improvements of 99.93% for matrix
projection, 90.24% for grouped execution, and 47.02% for metadata evaluation.

Representative continuous-screen results were:

| Diagnostic case | v0.2 | Release candidate | Ratio or status |
|---|---:|---:|---:|
| 50,000-group count | 104 +/- 1 ms | 10.7 +/- 0.7 ms | 0.10x |
| 50,000-group keys | 92.1 +/- 1 ms | 2.78 +/- 0.09 ms | 0.03x |
| Metadata semi-join | 120 +/- 1 ms | 113 +/- 2 ms | lower candidate median |
| Metadata anti-join | 122 +/- 1 ms | 113 +/- 100 ms | noisy, not separated |
| Chained extraction | 885 +/- 10 ms | 785 +/- 30 ms | 0.89x |
| Smoke grouped summarize | 8.03 +/- 0.03 ms | 9.61 +/- 0.01 ms | statistically inconclusive |

The small-workload grouped-summarize median remains 19.7% higher. Its ASV 99%
confidence intervals, 4.59-11.46 ms for v0.2 and 8.96-10.27 ms for the release
candidate, overlap; the continuous screen therefore did not classify the
two-sample diagnostic as a regression. It remains an unresolved small-workload
overhead rather than a resolved improvement. In contrast, the representative
20-group summarize case improved from 238 +/- 2 ms to 165 +/- 0.9 ms.

Fresh-process peak RSS, including fixture allocation, was:

| Case | v0.2 raw samples (KiB) | Release candidate raw samples (KiB) | Median change |
|---|---:|---:|---:|
| Dense projection | 572940, 573192, 573588 | 558020, 558308, 557896 | 2.65% lower |
| Grouped mutate | 302160, 302664, 301768 | 312528, 330816, 312392 | 3.43% higher |
| Metadata mutate | 792300, 792216, 792552 | 599956, 480164, 600052 | 24.28% lower |

One metadata candidate sample is substantially lower and one grouped candidate
sample is elevated, so only the three-process medians feed the manifest gate.
The exact current-host timing, RSS, normalized, and verdict files are retained
as ignored local artifacts under
`.asv/results/res-hpc-exe050.researchlumc.nl/`. Commits after `a6edfb2c` in the
release branch change documentation or repository configuration only; the
benchmarked `src/annplyr` tree is unchanged.

The raw ASV and RSS JSON can now be converted reproducibly before evaluation:

```text
python benchmarks/normalize_results.py --timing RAW_ASV.json \
  --peak-rss RAW_RSS.json --output NORMALIZED.json
python benchmarks/evaluate_manifest.py --baseline BASELINE.json \
  --candidate NORMALIZED.json --output VERDICT.json
```

## Correctness and measurement caveats

Timed methods contain lightweight realization, cardinality, or schema checks;
the primary benchmark gates and `test_v03_structural_invariance.py` provide the
stronger independent result comparison. `copy=False` call-cost cases check
shape and ordering without requiring a view, because v0.3 permits either a
view or a materialized object; view-realization cases are separate where an
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
ownership, and grouping state. That was the strongest downstream safeguard for
the frozen v0.3 candidate; the post-v0.3 integration evidence below supersedes
that limitation without changing the historical benchmark results.

## Post-v0.3 downstream integration evidence

Issue #13 adds a wholly synthetic, redistributable 8-cell by 6-feature AnnData
fixture generated from explicit values. No generated H5AD or Zarr dataset is
committed. Five tests cover the in-memory fixture, H5AD and Zarr v2 round trips,
split-and-reconstruct concatenation, and a Scanpy normalization/log1p/PCA
handoff. Every annplyr path compares exact matrix values and dtypes, pandas
extension dtypes and categorical metadata, axis order, CSR/CSC formats, aligned
container shapes, independent ownership, duplicate observation-name positional
identity, and retained grouping state.

The fixed-stack acceptance run passed all five tests on Python 3.14.5 with
AnnData 0.12.19 and Scanpy 1.12.2. Isolated stable runs then passed 5/5 on
Python 3.12.13, 3.13.13, and 3.14.5 with AnnData 0.13.2, Scanpy 1.12.3, NumPy
2.4.6, pandas 3.0.5, SciPy 1.18.0, and Zarr 3.3.0 writing the current
AnnData-default Zarr v2 representation. The advisory Python 3.14 prerelease
run also passed 5/5 with Scanpy 1.13.0a1 and NumPy 2.5.1; its remaining stack
matched the stable run. Dedicated mandatory CI repeats the stable matrix. The
prerelease-dependency consumer job remains advisory so an identified external
prerelease incompatibility cannot mask stable-environment status.
