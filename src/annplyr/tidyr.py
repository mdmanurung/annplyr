"""tidyr helpers for the pandas DataFrames annplyr returns.

These work on plain DataFrames, not AnnData, so they live apart from the
AnnData verbs at the package root::

    import annplyr as ap

    wide = ap.tidyr.pivot_wider(long, names_from="feature", values_from="value")
"""

from annplyr._verbs import (
    chop,
    drop_na,
    extract,
    fill,
    hoist,
    nest,
    pack,
    pivot_wider,
    separate,
    separate_rows,
    unchop,
    unite,
    unnest,
    unnest_longer,
    unnest_wider,
    unpack,
)

__all__ = [
    "chop",
    "drop_na",
    "extract",
    "fill",
    "hoist",
    "nest",
    "pack",
    "pivot_wider",
    "separate",
    "separate_rows",
    "unchop",
    "unite",
    "unnest",
    "unnest_longer",
    "unnest_wider",
    "unpack",
]
