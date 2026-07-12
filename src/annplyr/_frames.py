from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any, cast

import narwhals as nw
import numpy as np
import pandas as pd
from anndata import AnnData
from narwhals.exceptions import ColumnNotFoundError
from scipy import sparse

from annplyr._errors import SelectionError, UnknownColumnError, UnknownSourceError

ROW_NUMBER = "__annplyr_row_number__"
VIRTUAL_COLUMNS = {"obs_names", "var_names", ROW_NUMBER}


def obs_frame(adata: AnnData) -> pd.DataFrame:
    frame = cast(pd.DataFrame, adata.obs).copy()
    frame["obs_names"] = adata.obs_names.to_numpy()
    return frame


def var_frame(adata: AnnData) -> pd.DataFrame:
    frame = cast(pd.DataFrame, adata.var).copy()
    frame["var_names"] = adata.var_names.to_numpy()
    return frame


def x_frame(adata: AnnData, layer: str | None = None) -> pd.DataFrame:
    return adata.to_df(layer=layer)


def obsm_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        matrix = adata.obsm[key]
    except KeyError as exc:
        msg = f"Unknown obsm key: {key!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, adata.obs_names)


def varm_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        matrix = adata.varm[key]
    except KeyError as exc:
        msg = f"Unknown varm key: {key!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, adata.var_names)


def matrix_frame(matrix: Any, index: pd.Index) -> pd.DataFrame:
    if isinstance(matrix, pd.DataFrame):
        frame = matrix.copy()
        frame.index = index
        frame.columns = [str(column) for column in frame.columns]
        return frame

    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    return pd.DataFrame(values, index=index, columns=[str(i) for i in range(values.shape[1])])


def with_row_number(frame: pd.DataFrame) -> pd.DataFrame:
    work = frame.copy()
    work[ROW_NUMBER] = np.arange(1, len(work) + 1)
    return work


def drop_internal_columns(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.drop(columns=[ROW_NUMBER], errors="ignore")


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [value]
    if isinstance(value, Iterable):
        return list(value)
    return [value]


def selector_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, str):
        return [nw.col(value)]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if all(isinstance(item, str) for item in value):
            return [nw.col(*value)]
        return [nw.col(item) if isinstance(item, str) else item for item in value]
    return [value]


def expr_for(value: Any) -> Any:
    if isinstance(value, str):
        return nw.col(value)
    return value


def evaluate_select(frame: pd.DataFrame, selectors: Any) -> pd.DataFrame:
    if selectors is None:
        return drop_internal_columns(frame.copy())
    exprs = selector_list(selectors)
    if not exprs:
        return pd.DataFrame(index=frame.index)
    work = with_row_number(frame)
    columns = [str(column) for column in frame.columns]
    public_columns = [column for column in columns if column not in VIRTUAL_COLUMNS]
    resolved_exprs: list[Any] = []
    for expr in exprs:
        if hasattr(expr, "resolve"):
            names = expr.resolve(frame, columns, public_columns)
            if names:
                resolved_exprs.append(nw.col(*names))
        else:
            resolved_exprs.append(expr)
    if not resolved_exprs:
        return pd.DataFrame(index=frame.index)
    try:
        selected = nw.from_native(work).select(*resolved_exprs).to_native()
    except ColumnNotFoundError as exc:
        msg = f"Unknown column in selector: {exc}"
        raise UnknownColumnError(msg) from exc
    except Exception as exc:
        msg = f"Selection failed: {exc}"
        raise SelectionError(msg) from exc
    selected = drop_internal_columns(selected)
    if len(selected) == len(frame):
        selected.index = frame.index
    return selected


def evaluate_filter(frame: pd.DataFrame, predicates: Any) -> pd.Index:
    exprs = [expr_for(predicate) for predicate in as_list(predicates)]
    if not exprs:
        return frame.index
    work = with_row_number(frame)
    try:
        filtered = nw.from_native(work).filter(*exprs).to_native()
    except ColumnNotFoundError as exc:
        msg = f"Unknown column in filter predicate: {exc}"
        raise UnknownColumnError(msg) from exc
    except Exception as exc:
        msg = f"Filter failed: {exc}"
        raise SelectionError(msg) from exc
    return filtered.index


def evaluate_assignments(frame: pd.DataFrame, assignments: Mapping[str, Any] | None) -> pd.DataFrame:
    if not assignments:
        return pd.DataFrame(index=frame.index)
    exprs = [expr_for(expr).alias(name) for name, expr in assignments.items()]
    try:
        assigned = nw.from_native(with_row_number(frame)).select(*exprs).to_native()
    except ColumnNotFoundError as exc:
        msg = f"Unknown column in assignment: {exc}"
        raise UnknownColumnError(msg) from exc
    except Exception as exc:
        msg = f"Assignment failed: {exc}"
        raise SelectionError(msg) from exc
    if len(assigned) == len(frame):
        assigned.index = frame.index
    return assigned


def intersect_ordered(base: pd.Index, *indices: pd.Index) -> pd.Index:
    if not indices:
        return base
    keep = pd.Series(True, index=base)
    for index in indices:
        keep &= base.isin(index)
    return base[keep.to_numpy()]


def source_frame(adata: AnnData, source: str, key: str | None = None, layer: str | None = None) -> pd.DataFrame:
    if source == "obs":
        return obs_frame(adata)
    if source == "var":
        return var_frame(adata)
    if source == "x":
        return x_frame(adata, layer=layer)
    if source == "obsm":
        if key is None:
            msg = "obsm source requires a key"
            raise UnknownSourceError(msg)
        return obsm_frame(adata, key)
    if source == "varm":
        if key is None:
            msg = "varm source requires a key"
            raise UnknownSourceError(msg)
        return varm_frame(adata, key)
    msg = f"Unknown AnnData source: {source!r}"
    raise UnknownSourceError(msg)
