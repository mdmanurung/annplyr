"""Execute one full benchmark case and report process peak RSS as JSON."""

from __future__ import annotations

import argparse
import json
import resource

try:
    from . import bench_grouped as _bench_grouped
    from . import bench_metadata as _bench_metadata
    from . import bench_projection as _bench_projection
except ImportError:  # pragma: no cover - direct worker-script execution
    import bench_grouped as _bench_grouped
    import bench_metadata as _bench_metadata
    import bench_projection as _bench_projection


def _run(case: str) -> None:
    if case == "dense_projection":
        benchmark = _bench_projection.DenseProjection()
        benchmark.setup()
        benchmark.time_filter_three_features()
    elif case == "grouped_mutate":
        benchmark = _bench_grouped.GroupedTwenty()
        benchmark.setup()
        benchmark.time_mutate()
    elif case == "metadata_mutate":
        benchmark = _bench_metadata.WideMetadata()
        benchmark.setup()
        benchmark.time_mutate_ten_independent()
    else:
        raise ValueError(f"unknown peak-RSS case: {case}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=["dense_projection", "grouped_mutate", "metadata_mutate"])
    args = parser.parse_args()
    _run(args.case)
    # Linux reports ru_maxrss in KiB. Fixture setup is intentionally included
    # and the performance report must retain that caveat.
    print(json.dumps({"case": args.case, "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


if __name__ == "__main__":
    main()
