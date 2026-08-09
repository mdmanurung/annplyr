"""Normalize raw ASV timing and peak-RSS output for the manifest gate."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _case_names(name: str, parameters: list[list[str]]) -> list[str]:
    if not parameters:
        return [name]
    return [f"{name}({','.join(values)})" for values in itertools.product(*parameters)]


def _timing_samples(payload: dict[str, Any]) -> dict[str, list[float]]:
    columns = payload["result_columns"]
    parameter_index = columns.index("params")
    sample_index = columns.index("samples")
    normalized: dict[str, list[float]] = {}
    for name, result in payload["results"].items():
        if name.startswith("bench_smoke."):
            continue
        case_names = _case_names(name, result[parameter_index])
        sample_groups = result[sample_index]
        if len(case_names) != len(sample_groups):
            msg = f"ASV parameter/sample mismatch for {name}: {len(case_names)} cases, {len(sample_groups)} samples"
            raise ValueError(msg)
        for case_name, samples in zip(case_names, sample_groups, strict=True):
            if samples is None:
                raise ValueError(f"ASV benchmark has no raw samples: {case_name}")
            normalized[case_name] = [float(value) for value in samples]
    return normalized


def _normalize(timing_path: Path, peak_rss_path: Path) -> dict[str, Any]:
    timing = _load(timing_path)
    peak_rss = _load(peak_rss_path)
    timing_label = str(timing["commit_hash"])
    rss_label = str(peak_rss.get("label", timing_label))
    if not (timing_label.startswith(rss_label) or rss_label.startswith(timing_label)):
        raise ValueError(f"timing/RSS labels do not match: {timing_label!r} != {rss_label!r}")
    return {
        "label": timing_label,
        "peak_rss_kib": peak_rss["samples"],
        "result_hashes": peak_rss.get("result_hashes", {}),
        "source_peak_rss": str(peak_rss_path.resolve()),
        "source_timing": str(timing_path.resolve()),
        "timing_seconds": _timing_samples(timing),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timing", type=Path, required=True, help="Raw ASV result JSON")
    parser.add_argument("--peak-rss", type=Path, required=True, help="Raw run_peak_rss.py JSON")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    normalized = _normalize(args.timing, args.peak_rss)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
