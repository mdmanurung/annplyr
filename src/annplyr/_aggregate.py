"""Pseudobulk: sum or average a matrix over groups of cells.

This is the Python form of tidySingleCellExperiment's ``aggregate_cells``. The
result is an AnnData with one row per observed key combination, so the next
step (edgeR, DESeq2, pydeseq2, decoupler) receives an object, not a table.
"""

from __future__ import annotations

from typing import Any, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from annplyr._errors import AnnplyrError, UnknownSourceError
from annplyr._groups import GroupPlan, GroupSpec


def _source_matrix(adata: AnnData, source: Any) -> tuple[Any, Any]:
    """Return the matrix and its feature table for ``"x"``, ``"raw"`` or a layer."""
    name = getattr(source, "prefix", source)  # an annplyr.data Source carries its prefix
    if name in {"x", "X"}:
        return adata.X, adata.var
    if name == "raw":
        if adata.raw is None:
            raise UnknownSourceError("aggregate_cells source='raw' needs adata.raw, which is None")
        return adata.raw.X, adata.raw.var
    if isinstance(name, str) and name in adata.layers:
        return adata.layers[name], adata.var
    layers = ", ".join(repr(key) for key in adata.layers.keys() if key is not None)
    msg = f"aggregate_cells source={name!r} is not 'x', 'raw' or a layer; layers are: {layers or 'none'}"
    raise UnknownSourceError(msg)


def aggregate_cells_adata(
    adata: AnnData,
    *by: Any,
    source: Any,
    fun: Literal["sum", "mean"] = "sum",
) -> AnnData:
    """Sum (or average) ``source`` over each observed combination of ``by``.

    ``source`` has no default on purpose: ``X`` often holds scaled or
    log-normalized values, and a pseudobulk for differential expression must sum
    raw counts, for example ``source="counts"``. ``obs`` of the result holds the
    keys plus ``n_cells``; ``var`` is copied from the source matrix.
    """
    if adata.isbacked:
        raise AnnplyrError("aggregate_cells cannot read a backed AnnData; call .to_memory() first")
    if fun not in {"sum", "mean"}:
        raise AnnplyrError(f"aggregate_cells fun= must be 'sum' or 'mean', not {fun!r}")
    keys = by[0] if len(by) == 1 else by
    plan = GroupPlan.build(adata, GroupSpec.resolve(adata, obs=keys))
    matrix, features = _source_matrix(adata, source)

    n_groups = len(plan.keys)
    n_cells = np.bincount(plan.group_ids, minlength=n_groups)
    # a sum keeps the source dtype, so integer counts stay integers
    weights = np.ones(adata.n_obs, dtype=matrix.dtype) if fun == "sum" else 1.0 / n_cells[plan.group_ids]
    indicator = sp.csr_matrix((weights, (plan.group_ids, np.arange(adata.n_obs))), shape=(n_groups, adata.n_obs))
    values = indicator @ matrix  # sparse stays sparse, dense stays dense

    obs = plan.keys.copy()
    obs["n_cells"] = n_cells
    obs.index = pd.Index(obs[list(plan.spec.columns)].astype(str).agg("_".join, axis=1), name=None)
    return AnnData(X=values, obs=obs, var=features.copy())
