from __future__ import annotations

import json
from pathlib import Path
from runpy import run_path

_outside_noise = run_path(str(Path(__file__).parents[1] / "benchmarks" / "evaluate_manifest.py"))["_outside_noise"]
_normalize = run_path(str(Path(__file__).parents[1] / "benchmarks" / "normalize_results.py"))["_normalize"]


def test_noise_gate_distinguishes_overlapping_phases_from_separated_regression() -> None:
    baseline = [0.340, 0.345, 0.350, 0.420]
    overlapping_candidate = [0.345, 0.400, 0.410, 0.415]
    separated_candidate = [0.500, 0.510, 0.520, 0.530]

    assert not _outside_noise(baseline, overlapping_candidate)
    assert _outside_noise(baseline, separated_candidate)


def test_raw_benchmark_results_normalize_parameter_samples_and_skip_smoke(tmp_path: Path) -> None:
    timing_path = tmp_path / "timing.json"
    rss_path = tmp_path / "rss.json"
    timing_path.write_text(
        json.dumps(
            {
                "commit_hash": "abc123",
                "result_columns": ["result", "params", "samples"],
                "results": {
                    "bench_family.Case.time_plain": [0.2, [], [[0.1, 0.2]]],
                    "bench_family.Case.time_parameterized": [
                        [0.2, 0.3],
                        [["'csr'", "'csc'"]],
                        [[0.2, 0.3], [0.4, 0.5]],
                    ],
                    "bench_smoke.Smoke.time_case": [0.1, [], [[0.1]]],
                },
            }
        )
    )
    rss_path.write_text(json.dumps({"label": "abc123", "samples": {"dense_projection": [1, 2, 3]}}))

    normalized = _normalize(timing_path, rss_path)

    assert normalized["timing_seconds"] == {
        "bench_family.Case.time_plain": [0.1, 0.2],
        "bench_family.Case.time_parameterized('csr')": [0.2, 0.3],
        "bench_family.Case.time_parameterized('csc')": [0.4, 0.5],
    }
    assert normalized["peak_rss_kib"] == {"dense_projection": [1, 2, 3]}
