"""Nest an AnnData into one subset per group, and stack the subsets back.

This is tidySingleCellExperiment's ``nest(data = -cell_class)``: a table with
one row per group whose ``data`` column holds that group's object, so a
per-group analysis is a loop or a list comprehension over ``data``.
"""

from __future__ import annotations

from typing import Any

import anndata as ad
import pandas as pd
from anndata import AnnData

from annplyr._errors import SelectionError
from annplyr._groups import GroupPlan, GroupSpec
from annplyr._tidy import TidyAnnData


def nest_adata(adata: AnnData, *by: Any, name: str = "data") -> pd.DataFrame:
    """Return the observed keys, in first-seen order, plus one subset per group."""
    keys = by[0] if len(by) == 1 else by
    plan = GroupPlan.build(adata, GroupSpec.resolve(adata, obs=keys))
    if name in plan.keys.columns:
        raise SelectionError(f"nest name={name!r} is also a grouping key; pick another name")
    table = plan.keys.copy()
    table[name] = [TidyAnnData.wrap(adata[positions].copy()) for positions in plan.positions]
    return table


def unnest_adata(table: pd.DataFrame, name: str = "data") -> AnnData:
    """Stack the objects in ``table[name]`` back into one AnnData, in row order.

    ``anndata.concat`` keeps the obs and var tables and the aligned matrices,
    with ``join="outer"`` so a per-group step that added an obs column keeps it.
    """
    if name not in table.columns:
        raise SelectionError(f"unnest needs a {name!r} column of AnnData objects")
    return TidyAnnData.wrap(ad.concat(list(table[name]), join="outer", merge="same", uns_merge="same"))
