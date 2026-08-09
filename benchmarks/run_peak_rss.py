"""Run three fresh subprocess peak-RSS samples per performance family."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

CASES = ("dense_projection", "dense_reduction", "grouped_mutate", "metadata_mutate")
THREAD_VARIABLES = (
    "BLIS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--source-root", type=Path)
    args = parser.parse_args()
    if args.processes != 3:
        raise SystemExit("the frozen manifest requires exactly three fresh processes")

    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env.update(dict.fromkeys(THREAD_VARIABLES, "1"))
    source_root = (args.source_root or repo_root).resolve()
    source_path = source_root / "src"
    python_paths = [str(source_path), str(repo_root)]
    if existing_pythonpath := env.get("PYTHONPATH"):
        python_paths.append(existing_pythonpath)
    env["PYTHONPATH"] = os.pathsep.join(python_paths)
    samples: dict[str, list[int]] = {case: [] for case in CASES}
    result_hashes: dict[str, set[str]] = {case: set() for case in CASES}
    for case in CASES:
        for _ in range(args.processes):
            completed = subprocess.run(
                [args.python, "-m", "benchmarks.peak_rss_worker", case],
                check=True,
                capture_output=True,
                text=True,
                env=env,
                cwd=repo_root,
            )
            payload = json.loads(completed.stdout)
            samples[case].append(int(payload["peak_rss_kib"]))
            if payload["result_hash"]:
                result_hashes[case].add(str(payload["result_hash"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "label": args.label,
                "fresh_processes": args.processes,
                "unit": "KiB",
                "setup_allocation_included": True,
                "samples": samples,
                "result_hashes": {case: sorted(values) for case, values in result_hashes.items() if values},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
