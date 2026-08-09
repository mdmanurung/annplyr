"""Evaluate normalized timing/RSS samples against the frozen v0.3 gates."""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from typing import Any


def _median(values: list[float]) -> float:
    if not values:
        raise ValueError("sample list is empty")
    return float(statistics.median(values))


def _mad(values: list[float]) -> float:
    center = _median(values)
    return float(statistics.median(abs(value - center) for value in values))


def _outside_noise(baseline: list[float], candidate: list[float]) -> bool:
    baseline_median = _median(baseline)
    candidate_median = _median(candidate)
    scale = max(1.4826 * _mad(baseline), 1.4826 * _mad(candidate))
    raw_ranges_overlap = max(min(baseline), min(candidate)) <= min(max(baseline), max(candidate))
    return candidate_median - baseline_median > scale and not raw_ranges_overlap


def _matches(samples: dict[str, list[float]], prefix: str) -> list[str]:
    return sorted(
        name for name in samples if name == prefix or name.startswith(f"{prefix}(") or name.startswith(f"{prefix}-")
    )


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _evaluate(
    manifest: dict[str, Any],
    baseline: dict[str, Any],
    candidate: dict[str, Any],
    *,
    improvement_families: set[str] | None = None,
) -> dict[str, Any]:
    """Evaluate regressions globally and improvement gates for selected families."""

    configured_families = set(manifest["families"])
    required_families = configured_families if improvement_families is None else improvement_families
    if unknown := required_families - configured_families:
        names = ", ".join(sorted(unknown))
        raise ValueError(f"unknown improvement families: {names}")

    minimum = int(manifest["measurement"]["timing_repeats_minimum"])
    improvement_gate = float(manifest["gates"]["family_improvement_fraction"])
    regression_gate = float(manifest["gates"]["primary_regression_fraction"])

    timing_base: dict[str, list[float]] = baseline["timing_seconds"]
    timing_new: dict[str, list[float]] = candidate["timing_seconds"]
    rss_base: dict[str, list[float]] = baseline["peak_rss_kib"]
    rss_new: dict[str, list[float]] = candidate["peak_rss_kib"]
    hashes_base: dict[str, list[str]] = baseline.get("result_hashes", {})
    hashes_new: dict[str, list[str]] = candidate.get("result_hashes", {})

    sample_errors: list[str] = []
    for label, samples in (("baseline", timing_base), ("candidate", timing_new)):
        for name, values in samples.items():
            if len(values) < minimum:
                sample_errors.append(f"{label} {name} has {len(values)} timing samples; need {minimum}")
    for label, samples in (("baseline", rss_base), ("candidate", rss_new)):
        for name, values in samples.items():
            if len(values) != 3:
                sample_errors.append(f"{label} {name} has {len(values)} RSS samples; need exactly 3")
    for name in sorted(set(hashes_base) | set(hashes_new)):
        if hashes_base.get(name) != hashes_new.get(name):
            sample_errors.append(f"result hash mismatch for {name}: {hashes_base.get(name)} != {hashes_new.get(name)}")

    families: dict[str, Any] = {}
    regressions: list[dict[str, Any]] = []
    for family, config in manifest["families"].items():
        improvements: list[dict[str, Any]] = []
        for prefix in config["primary"]:
            for name in _matches(timing_base, prefix):
                if name not in timing_new:
                    sample_errors.append(f"candidate missing timing benchmark {name}")
                    continue
                base_values = timing_base[name]
                new_values = timing_new[name]
                base_median = _median(base_values)
                new_median = _median(new_values)
                fraction = (base_median - new_median) / base_median
                record = {
                    "metric": "time",
                    "name": name,
                    "baseline_median": base_median,
                    "candidate_median": new_median,
                    "improvement_fraction": fraction,
                }
                improvements.append(record)
                regression = (new_median - base_median) / base_median
                if regression > regression_gate and _outside_noise(base_values, new_values):
                    regressions.append({**record, "regression_fraction": regression})

        rss_name = config["rss"]
        if rss_name not in rss_base or rss_name not in rss_new:
            sample_errors.append(f"missing RSS case {rss_name}")
        else:
            base_median = _median(rss_base[rss_name])
            new_median = _median(rss_new[rss_name])
            improvements.append(
                {
                    "metric": "peak_rss",
                    "name": rss_name,
                    "baseline_median": base_median,
                    "candidate_median": new_median,
                    "improvement_fraction": (base_median - new_median) / base_median,
                }
            )
        best = max((item["improvement_fraction"] for item in improvements), default=float("-inf"))
        improvement_required = family in required_families
        families[family] = {
            "pass": not improvement_required or best >= improvement_gate,
            "improvement_required": improvement_required,
            "best_improvement_fraction": best,
            "cases": improvements,
        }

    verdict = not sample_errors and not regressions and all(item["pass"] for item in families.values())
    return {
        "pass": verdict,
        "improvement_families": sorted(required_families),
        "sample_errors": sample_errors,
        "primary_regressions_outside_noise": regressions,
        "families": families,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("manifest.json"))
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--improvement-family",
        action="append",
        dest="improvement_families",
        help="Require the >=20%% improvement gate only for this family; repeat as needed (default: all)",
    )
    args = parser.parse_args()

    manifest = _load(args.manifest)
    baseline = _load(args.baseline)
    candidate = _load(args.candidate)
    try:
        report = _evaluate(
            manifest,
            baseline,
            candidate,
            improvement_families=set(args.improvement_families) if args.improvement_families else None,
        )
    except ValueError as error:
        parser.error(str(error))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
