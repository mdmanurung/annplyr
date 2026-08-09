from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any, cast

import narwhals as nw
import numpy as np
import pandas as pd
from anndata import AnnData
from narwhals.exceptions import ColumnNotFoundError
from scipy import sparse

from annplyr._errors import (
    DuplicateNameError,
    SelectionError,
    SizeMismatchError,
    UnknownColumnError,
    UnknownSourceError,
)
from annplyr._expr import AnnplyrExpr, to_narwhals

ROW_NUMBER = "__annplyr_row_number__"
OBS_NAMES = "__annplyr_obs_names__"
VAR_NAMES = "__annplyr_var_names__"
VIRTUAL_COLUMNS = {OBS_NAMES, VAR_NAMES, ROW_NUMBER}
VIRTUAL_ATTR = "annplyr_virtual_columns"


def obs_frame(adata: AnnData) -> pd.DataFrame:
    frame = cast(pd.DataFrame, adata.obs).copy(deep=False)
    frame.index = pd.RangeIndex(len(frame))
    frame[OBS_NAMES] = adata.obs_names.to_numpy()
    frame.attrs[VIRTUAL_ATTR] = {OBS_NAMES}
    return frame


def var_frame(adata: AnnData) -> pd.DataFrame:
    frame = cast(pd.DataFrame, adata.var).copy(deep=False)
    frame.index = pd.RangeIndex(len(frame))
    frame[VAR_NAMES] = adata.var_names.to_numpy()
    frame.attrs[VIRTUAL_ATTR] = {VAR_NAMES}
    return frame


def x_frame(adata: AnnData, layer: str | None = None) -> pd.DataFrame:
    try:
        matrix = adata.layers[layer] if layer is not None else adata.X
    except KeyError as exc:
        msg = f"Unknown layer: {layer!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, pd.RangeIndex(adata.n_obs), columns=adata.var_names)


def raw_frame(adata: AnnData) -> pd.DataFrame:
    if adata.raw is None:
        msg = "AnnData object has no raw matrix"
        raise UnknownSourceError(msg)
    return matrix_frame(adata.raw.X, pd.RangeIndex(adata.n_obs), columns=adata.raw.var_names)


def obsm_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        matrix = adata.obsm[key]
    except KeyError as exc:
        msg = f"Unknown obsm key: {key!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, pd.RangeIndex(adata.n_obs))


def varm_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        matrix = adata.varm[key]
    except KeyError as exc:
        msg = f"Unknown varm key: {key!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, pd.RangeIndex(adata.n_vars))


def obsp_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        matrix = adata.obsp[key]
    except KeyError as exc:
        msg = f"Unknown obsp key: {key!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, adata.obs_names, columns=adata.obs_names)


def varp_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        matrix = adata.varp[key]
    except KeyError as exc:
        msg = f"Unknown varp key: {key!r}"
        raise UnknownSourceError(msg) from exc
    return matrix_frame(matrix, adata.var_names, columns=adata.var_names)


def uns_frame(adata: AnnData, key: str) -> pd.DataFrame:
    try:
        value = adata.uns[key]
    except KeyError as exc:
        msg = f"Unknown uns key: {key!r}"
        raise UnknownSourceError(msg) from exc
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if isinstance(value, pd.Series):
        return value.to_frame().copy()
    if isinstance(value, Mapping):
        try:
            return pd.DataFrame(value)
        except Exception as exc:
            msg = f"uns key {key!r} cannot be represented as a pandas DataFrame"
            raise UnknownSourceError(msg) from exc
    msg = f"uns key {key!r} cannot be represented as a pandas DataFrame"
    raise UnknownSourceError(msg)


def matrix_frame(matrix: Any, index: pd.Index, *, columns: pd.Index | Sequence[str] | None = None) -> pd.DataFrame:
    if isinstance(matrix, pd.DataFrame):
        frame = matrix.copy()
        frame.index = index
        frame.columns = [str(column) for column in frame.columns]
        return frame

    if sparse.issparse(matrix):
        spmatrix = matrix.tocsc()
        sparse_array = cast(Any, pd.arrays.SparseArray)
        names = (
            [str(column) for column in columns] if columns is not None else [str(i) for i in range(spmatrix.shape[1])]
        )
        series = [
            pd.Series(sparse_array.from_spmatrix(spmatrix[:, i]), index=index, name=names[i])
            for i in range(spmatrix.shape[1])
        ]
        return pd.concat(series, axis=1) if series else pd.DataFrame(index=index, columns=names)

    values = np.asarray(matrix)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    if columns is None:
        columns = [str(i) for i in range(values.shape[1])]
    return pd.DataFrame(values, index=index, columns=[str(column) for column in columns])


#: pandas only ships sparse binary kernels for these subtypes.
SPARSE_ARITHMETIC_SUBTYPES = frozenset({np.dtype(np.float64), np.dtype(np.int64), np.dtype(bool)})


def _needs_dense_evaluation(dtype: Any) -> bool:
    return isinstance(dtype, pd.SparseDtype) and np.dtype(dtype.subtype) not in SPARSE_ARITHMETIC_SUBTYPES


def prepare_sparse_for_arithmetic(frame: pd.DataFrame) -> pd.DataFrame:
    """Densify sparse columns pandas cannot compute with.

    Single-cell matrices are almost always `float32`, but pandas implements
    sparse binary operations for `float64`, `int64`, and `bool` only. Without
    this step, combining two projected matrix columns in one expression fails
    with an opaque ``sparse_add_float32`` attribute error.

    Densifying preserves the subtype, so a sparse source yields exactly the
    dtype its dense equivalent would. Only evaluation scratch frames are
    converted, so exported frames keep their sparse columns, and only the
    already-projected columns are materialized.
    """
    positions = [position for position, dtype in enumerate(frame.dtypes) if _needs_dense_evaluation(dtype)]
    if not positions:
        return frame
    # Positional assignment keeps duplicate axis names addressable.
    prepared = frame.copy(deep=False)
    for position in positions:
        column = cast(Any, frame.iloc[:, position])
        prepared.isetitem(position, column.sparse.to_dense())
    return prepared


def with_row_number(frame: pd.DataFrame) -> pd.DataFrame:
    # Evaluation only adds or replaces whole columns, so a shallow scratch
    # frame avoids copying every metadata block without exposing source writes.
    work = frame.copy(deep=False)
    work[ROW_NUMBER] = np.arange(1, len(work) + 1)
    virtual = set(work.attrs.get(VIRTUAL_ATTR, set()))
    virtual.add(ROW_NUMBER)
    work.attrs[VIRTUAL_ATTR] = virtual
    return work


def drop_internal_columns(frame: pd.DataFrame) -> pd.DataFrame:
    virtual = set(frame.attrs.get(VIRTUAL_ATTR, set())) | {ROW_NUMBER}
    return frame.drop(columns=[column for column in virtual if column in frame.columns], errors="ignore")


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
        value = nw.col(value)
    return to_narwhals(value)


def evaluate_select(frame: pd.DataFrame, selectors: Any) -> pd.DataFrame:
    if selectors is None:
        return drop_internal_columns(frame.copy())
    exprs = selector_list(selectors)
    if not exprs:
        return pd.DataFrame(index=frame.index)
    work = with_row_number(frame)
    columns = [str(column) for column in frame.columns]
    virtual_columns = set(frame.attrs.get(VIRTUAL_ATTR, set())) | VIRTUAL_COLUMNS
    public_columns = [column for column in columns if column not in virtual_columns]
    resolved_exprs: list[Any] = []
    for expr in exprs:
        try:
            if hasattr(expr, "resolve"):
                names = expr.resolve(frame, columns, public_columns)
                if names:
                    resolved_exprs.append(nw.col(*names))
            else:
                resolved_exprs.append(to_narwhals(expr))
        except UnknownColumnError:
            raise
        except Exception as exc:
            msg = f"Selection failed while resolving selector: {exc}"
            raise SelectionError(msg) from exc
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
    exprs = [_predicate_expr(frame, predicate) for predicate in as_list(predicates)]
    if not exprs:
        return frame.index
    work = with_row_number(prepare_sparse_for_arithmetic(frame))
    try:
        filtered = nw.from_native(work).filter(*exprs).to_native()
    except ColumnNotFoundError as exc:
        msg = f"Unknown column in filter predicate: {exc}"
        raise UnknownColumnError(msg) from exc
    except Exception as exc:
        msg = f"Filter failed: {exc}"
        raise SelectionError(msg) from exc
    return filtered.index


def _predicate_expr(frame: pd.DataFrame, predicate: Any) -> Any:
    if hasattr(predicate, "to_expr"):
        return to_narwhals(predicate.to_expr(frame))
    return expr_for(predicate)


def expand_assignments(frame: pd.DataFrame, assignments: Mapping[str, Any] | Any | None) -> dict[str, Any]:
    if not assignments:
        return {}
    if hasattr(assignments, "expand"):
        return dict(assignments.expand(frame))
    expanded: dict[str, Any] = {}
    for name, expr in assignments.items():
        if hasattr(expr, "expand"):
            for expanded_name, expanded_expr in expr.expand(frame).items():
                _add_expanded_assignment(expanded, expanded_name, expanded_expr)
        else:
            _add_expanded_assignment(expanded, str(name), expr)
    return expanded


def _add_expanded_assignment(expanded: dict[str, Any], name: str, expr: Any) -> None:
    if name in expanded:
        msg = f"Duplicate assignment output name: {name!r}"
        raise DuplicateNameError(msg)
    expanded[name] = expr


def evaluate_assignments(frame: pd.DataFrame, assignments: Mapping[str, Any] | Any | None) -> pd.DataFrame:
    assignments = expand_assignments(frame, assignments)
    if not assignments:
        return pd.DataFrame(index=frame.index)
    work = with_row_number(prepare_sparse_for_arithmetic(frame))
    assigned = pd.DataFrame(index=frame.index)
    for batch in _assignment_batches(assignments):
        try:
            expressions = [expr_for(expr).alias(name) for name, expr in batch]
            result = nw.from_native(work).select(*expressions).to_native()
        except ColumnNotFoundError as exc:
            msg = f"Unknown column in assignment: {exc}"
            raise UnknownColumnError(msg) from exc
        except Exception as exc:
            msg = f"Assignment failed: {exc}"
            raise SelectionError(msg) from exc
        for name, _ in batch:
            series = _aligned_assignment_series(result[name], frame.index, name=name)
            assigned[name] = series
            work[name] = series
    return assigned


def _assignment_batches(assignments: Mapping[str, Any]) -> list[list[tuple[str, Any]]]:
    """Batch only contiguous wrapped expressions proven mutually independent."""
    batches: list[list[tuple[str, Any]]] = []
    current: list[tuple[str, Any]] = []
    current_outputs: set[str] = set()
    current_dependencies: set[str] = set()
    current_cardinality: str | None = None

    def flush() -> None:
        nonlocal current, current_outputs, current_dependencies, current_cardinality
        if current:
            batches.append(current)
        current = []
        current_outputs = set()
        current_dependencies = set()
        current_cardinality = None

    for name, expression in assignments.items():
        batchable = (
            isinstance(expression, AnnplyrExpr)
            and expression.dependencies is not None
            and expression.output_width == 1
            and expression.cardinality in {"row", "scalar"}
        )
        dependencies = set(expression.dependencies or ()) if isinstance(expression, AnnplyrExpr) else set()
        compatible = (
            batchable
            and (current_cardinality is None or current_cardinality == expression.cardinality)
            and not dependencies.intersection(current_outputs)
            and name not in current_dependencies
        )
        if not compatible:
            flush()
        if batchable:
            current.append((name, expression))
            current_outputs.add(name)
            current_dependencies.update(dependencies)
            current_cardinality = expression.cardinality
        else:
            batches.append([(name, expression)])
    flush()
    return batches


def _aligned_assignment_series(series: pd.Series, index: pd.Index, *, name: str) -> pd.Series:
    if len(series) == len(index):
        out = series.copy()
        out.index = index
        return out
    if len(series) == 1:
        return pd.Series([series.iloc[0]] * len(index), index=index, name=name)
    msg = f"Assignment {name!r} returned {len(series)} rows for an axis of length {len(index)}"
    raise SizeMismatchError(msg)


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
    if source == "raw":
        if layer is not None:
            msg = "raw source does not support layer"
            raise UnknownSourceError(msg)
        return raw_frame(adata)
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
    if source == "obsp":
        if key is None:
            msg = "obsp source requires a key"
            raise UnknownSourceError(msg)
        return obsp_frame(adata, key)
    if source == "varp":
        if key is None:
            msg = "varp source requires a key"
            raise UnknownSourceError(msg)
        return varp_frame(adata, key)
    if source == "uns":
        if key is None:
            msg = "uns source requires a key"
            raise UnknownSourceError(msg)
        return uns_frame(adata, key)
    msg = f"Unknown AnnData source: {source!r}"
    raise UnknownSourceError(msg)
