from __future__ import annotations

from pathlib import Path
from runpy import run_path

_outside_noise = run_path(str(Path(__file__).parents[1] / "benchmarks" / "evaluate_manifest.py"))["_outside_noise"]


def test_noise_gate_distinguishes_overlapping_phases_from_separated_regression() -> None:
    baseline = [0.340, 0.345, 0.350, 0.420]
    overlapping_candidate = [0.345, 0.400, 0.410, 0.415]
    separated_candidate = [0.500, 0.510, 0.520, 0.530]

    assert not _outside_noise(baseline, overlapping_candidate)
    assert _outside_noise(baseline, separated_candidate)
