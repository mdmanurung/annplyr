# Issue 10 bounded-reduction performance evidence

This report compares the last merged `main` before issue 10
(`981204e482ec45aa5356e555c4ec812d819b61be`) with the implementation and
benchmark commit (`6a3937395e8392254bd2757b9659cbda347d256a`). Both sides ran on
`res-hpc-exe050.researchlumc.nl` with the Python 3.14.5 executable and exact
package stack frozen in `benchmarks/manifest.json`. BLIS, OpenBLAS, OpenMP,
MKL, NumExpr, and vecLib thread counts were fixed at one.

The new case summarizes every feature in a deterministic 100,000 by 500
`float32` dense matrix. Fixture construction is excluded from ASV timing. The
oracle independently computes one NumPy mean per column and requires the exact
SHA-256 result hash
`cb89489dbc9e8f4485803d3311cd662a0879f8178c5328df9ed47f8823bece5a`.
The baseline and candidate produced that same hash in every fresh-process RSS
run.

## Target result

The exact dense-reduction timing samples in seconds were:

- baseline: `0.657457531, 0.669028748, 1.699717453, 1.023851365,
  0.811600719, 0.773238842, 0.745438653, 0.626717608, 0.616092381,
  0.626487509, 0.635981191, 0.655441234, 0.671601439, 0.695199905`;
- candidate: `0.624950536, 0.636253546, 0.649653471, 0.672218355,
  0.696109598, 0.692794364, 0.731372488, 0.617569094, 0.620925128,
  0.641435086, 0.656598042, 0.664862120, 0.673265718, 0.715967264`.

The timing median changed from 0.670315 s (range 0.616092-1.699717) to
0.660730 s (range 0.617569-0.731372), a 1.43% improvement. This is effectively
timing-neutral evidence, not a basis for a broad speed claim.

Fresh-process peak RSS, including fixture allocation, was:

| Revision | Exact samples (KiB) | Median (KiB) | Range (KiB) |
|---|---:|---:|---:|
| Baseline | 985076, 984976, 985728 | 985076 | 984976-985728 |
| Candidate | 585388, 591388, 584020 | 585388 | 584020-591388 |

The median peak RSS reduction is 40.57% (399,688 KiB). It clears the frozen
20% family-improvement threshold without changing the exact result hash.

## Existing-family regression screen

Negative changes are improvements. Every primary timing case retained at least
14 raw samples, and every RSS case used exactly three new processes.

| Case | Baseline median [range] | Candidate median [range] | Change |
|---|---:|---:|---:|
| Grouped arrange | 1.14244 s [0.963959, 1.86764] | 1.02162 s [0.969404, 1.10660] | -10.58% |
| Grouped filter | 0.678819 s [0.629967, 0.746237] | 0.696111 s [0.663401, 0.728234] | +2.55% |
| Grouped mutate | 0.404620 s [0.358644, 0.458352] | 0.394924 s [0.367164, 0.429833] | -2.40% |
| Grouped slice head | 0.267356 s [0.252179, 0.282546] | 0.268884 s [0.250606, 0.292969] | +0.57% |
| Grouped summarize | 0.223175 s [0.209021, 0.355743] | 0.223541 s [0.208093, 0.237615] | +0.16% |
| Metadata inner join | 0.650787 s [0.531341, 1.88487] | 0.604299 s [0.525197, 1.56471] | -7.14% |
| Metadata left join | 1.12371 s [0.801540, 1.93691] | 0.914412 s [0.826235, 1.06446] | -18.63% |
| Metadata mutate ten | 0.405872 s [0.216468, 1.07234] | 0.236686 s [0.224682, 0.271961] | -41.68% |
| Metadata select twenty | 0.252839 s [0.238682, 0.884309] | 0.236800 s [0.228688, 0.267672] | -6.34% |
| Backed CSC selected read | 20.4925 s [19.8898, 22.2156] | 20.0675 s [18.5106, 21.5856] | -2.07% |
| Backed CSR selected read | 20.0097 s [19.8177, 20.7361] | 20.0251 s [18.6824, 21.3697] | +0.08% |
| Backed dense selected read | 0.207387 s [0.202109, 0.217188] | 0.206278 s [0.203895, 0.211589] | -0.53% |
| Dense export ten | 0.157894 s [0.151268, 0.869221] | 0.155706 s [0.149167, 0.169813] | -1.39% |
| Dense filter three | 0.429102 s [0.415937, 0.547365] | 0.431514 s [0.422875, 0.461707] | +0.56% |
| Sparse CSC filter three | 21.5657 s [21.1184, 22.6879] | 20.6429 s [20.2136, 21.8876] | -4.28% |
| Sparse CSR filter three | 21.6942 s [21.4962, 22.8987] | 21.4310 s [20.2839, 23.0443] | -1.21% |
| Sparse CSC export ten | 20.1098 s [19.7755, 20.5229] | 19.6144 s [19.0910, 21.1464] | -2.46% |
| Sparse CSR export ten | 20.0496 s [19.5585, 21.1734] | 19.9813 s [19.4153, 21.1910] | -0.34% |
| Dense reduction | 0.670315 s [0.616092, 1.69972] | 0.660730 s [0.617569, 0.731372] | -1.43% |
| RSS dense projection | 624496 KiB [616460, 626304] | 618500 KiB [616728, 622408] | -0.96% |
| RSS grouped mutate | 378216 KiB [378024, 379280] | 378460 KiB [370828, 379060] | +0.06% |
| RSS metadata mutate | 546252 KiB [538260, 658128] | 547180 KiB [547024, 665996] | +0.17% |

The manifest evaluator found no sample errors and no primary regressions outside
its median-absolute-deviation plus raw-range-overlap noise rule. Incremental
issue comparisons can name the family required to clear the 20% improvement
gate with `--improvement-family`; omitting that option preserves the historical
release-wide requirement that every family improve. Issue 10 focused
`bounded_reductions`, while projection, grouped execution, and metadata
evaluation remained subject to the unchanged global regression gate. The
focused verdict passed.

The complete normalized timing/RSS samples and machine-readable verdict are
committed as
`benchmarks/evidence/issue-10/baseline-normalized.json`,
`candidate-normalized.json`, and `manifest-evaluation.json`.

## Limitations

- This is one fixed HPC runner, one pinned stack, and one run period. It does
  not establish performance on other hardware or dependency versions.
- Backed timings use warm filesystem caches. Setup is excluded from ASV timing,
  while fresh-process peak RSS deliberately includes fixture allocation.
- Several baseline timing ranges contain high outliers. The report therefore
  keeps every sample and applies the manifest noise policy rather than selecting
  favorable observations.
- The 100,000 by 500 case measures canonical dense means. Correctness tests,
  rather than this benchmark, cover every reducer and dense, CSR, CSC, backed
  dense, backed CSR, and backed CSC storage.
- Exact rank and distinct reducers may retain one reduction vector of state.
  Opaque or cross-row custom expressions remain on the documented conservative
  eager path because they cannot be safely decomposed without changing results.
