from __future__ import annotations

import json
from pathlib import Path
from runpy import run_path

_outside_noise = run_path(str(Path(__file__).parents[1] / "benchmarks" / "evaluate_manifest.py"))["_outside_noise"]
_evaluate = run_path(str(Path(__file__).parents[1] / "benchmarks" / "evaluate_manifest.py"))["_evaluate"]
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
    rss_path.write_text(
        json.dumps(
            {
                "label": "abc123",
                "samples": {"dense_projection": [1, 2, 3]},
                "result_hashes": {"dense_reduction": ["deadbeef"]},
            }
        )
    )

    normalized = _normalize(timing_path, rss_path)

    assert normalized["timing_seconds"] == {
        "bench_family.Case.time_plain": [0.1, 0.2],
        "bench_family.Case.time_parameterized('csr')": [0.2, 0.3],
        "bench_family.Case.time_parameterized('csc')": [0.4, 0.5],
    }
    assert normalized["peak_rss_kib"] == {"dense_projection": [1, 2, 3]}
    assert normalized["result_hashes"] == {"dense_reduction": ["deadbeef"]}


def test_manifest_can_focus_improvement_without_weakening_regression_gate() -> None:
    manifest = {
        "measurement": {"timing_repeats_minimum": 7},
        "gates": {"family_improvement_fraction": 0.20, "primary_regression_fraction": 0.10},
        "families": {
            "changed": {"primary": ["bench.changed"], "rss": "changed"},
            "unchanged": {"primary": ["bench.unchanged"], "rss": "unchanged"},
        },
    }
    baseline = {
        "timing_seconds": {"bench.changed": [10.0] * 7, "bench.unchanged": [10.0] * 7},
        "peak_rss_kib": {"changed": [100, 100, 100], "unchanged": [100, 100, 100]},
    }
    candidate = {
        "timing_seconds": {"bench.changed": [7.0] * 7, "bench.unchanged": [10.0] * 7},
        "peak_rss_kib": {"changed": [100, 100, 100], "unchanged": [100, 100, 100]},
    }

    release_wide = _evaluate(manifest, baseline, candidate)
    focused = _evaluate(manifest, baseline, candidate, improvement_families={"changed"})

    assert not release_wide["pass"]
    assert focused["pass"]
    assert focused["families"]["changed"]["improvement_required"]
    assert not focused["families"]["unchanged"]["improvement_required"]


def test_focused_manifest_still_rejects_unrelated_regression() -> None:
    manifest = {
        "measurement": {"timing_repeats_minimum": 7},
        "gates": {"family_improvement_fraction": 0.20, "primary_regression_fraction": 0.10},
        "families": {
            "changed": {"primary": ["bench.changed"], "rss": "changed"},
            "unchanged": {"primary": ["bench.unchanged"], "rss": "unchanged"},
        },
    }
    baseline = {
        "timing_seconds": {"bench.changed": [10.0] * 7, "bench.unchanged": [10.0] * 7},
        "peak_rss_kib": {"changed": [100, 100, 100], "unchanged": [100, 100, 100]},
    }
    candidate = {
        "timing_seconds": {"bench.changed": [7.0] * 7, "bench.unchanged": [20.0] * 7},
        "peak_rss_kib": {"changed": [100, 100, 100], "unchanged": [100, 100, 100]},
    }

    report = _evaluate(manifest, baseline, candidate, improvement_families={"changed"})

    assert not report["pass"]
    assert report["primary_regressions_outside_noise"]
