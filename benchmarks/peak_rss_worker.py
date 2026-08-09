"""Execute one full benchmark case and report process peak RSS as JSON."""

from __future__ import annotations

import argparse
import json
import resource

try:
    from . import bench_grouped as _bench_grouped
    from . import bench_metadata as _bench_metadata
    from . import bench_projection as _bench_projection
    from . import bench_reductions as _bench_reductions
except ImportError:  # pragma: no cover - direct worker-script execution
    import bench_grouped as _bench_grouped
    import bench_metadata as _bench_metadata
    import bench_projection as _bench_projection
    import bench_reductions as _bench_reductions


def _run(case: str) -> str:
    if case == "dense_projection":
        benchmark = _bench_projection.DenseProjection()
        benchmark.setup()
        benchmark.time_filter_three_features()
        return ""
    elif case == "grouped_mutate":
        benchmark = _bench_grouped.GroupedTwenty()
        benchmark.setup()
        benchmark.time_mutate()
        return ""
    elif case == "metadata_mutate":
        benchmark = _bench_metadata.WideMetadata()
        benchmark.setup()
        benchmark.time_mutate_ten_independent()
        return ""
    elif case == "dense_reduction":
        benchmark = _bench_reductions.DenseReductions()
        benchmark.setup()
        benchmark.time_mean_all_features()
        return benchmark.result_hash
    else:
        raise ValueError(f"unknown peak-RSS case: {case}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "case",
        choices=["dense_projection", "dense_reduction", "grouped_mutate", "metadata_mutate"],
    )
    args = parser.parse_args()
    result_hash = _run(args.case)
    # Linux reports ru_maxrss in KiB. Fixture setup is intentionally included
    # and the performance report must retain that caveat.
    print(
        json.dumps(
            {
                "case": args.case,
                "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                "result_hash": result_hash,
            }
        )
    )


if __name__ == "__main__":
    main()
