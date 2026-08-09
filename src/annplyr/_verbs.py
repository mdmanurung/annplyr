from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from typing import Any, cast

import narwhals as nw
import numpy as np
import pandas as pd
from anndata import AnnData

from annplyr._errors import (
    AnnplyrError,
    DuplicateNameError,
    IncompatibleAxisError,
    JoinRelationshipError,
    NameRepairError,
    SelectionError,
    UnknownColumnError,
    UnknownSourceError,
)
from annplyr._expr import Desc, to_narwhals
from annplyr._frames import (
    OBS_NAMES,
    VAR_NAMES,
    VIRTUAL_ATTR,
    evaluate_assignments,
    evaluate_filter,
    evaluate_select,
    expand_assignments,
    intersect_ordered,
    obs_frame,
    obsm_frame,
    prepare_sparse_for_arithmetic,
    raw_frame,
    source_frame,
    var_frame,
    varm_frame,
    with_row_number,
    x_frame,
)
from annplyr._groups import GroupPlan
from annplyr._reductions import ReductionPlan, SummaryGroupPlan, summarize_chunked
from annplyr._sources import (
    DEFAULT_REDUCTION_CHUNK_VALUES,
    RequestPlanner,
    request_dependencies,
    source_adapter,
)


def _decorate_projected_frame(
    frame: pd.DataFrame,
    adata: AnnData,
    *,
    axis: str,
    request: Any,
    row_positions: np.ndarray | None = None,
) -> pd.DataFrame:
    """Attach an axis-name virtual only when dependency metadata requires it."""
    dependencies = request_dependencies(request, frame.iloc[:0, :])
    virtual_name = OBS_NAMES if axis == "obs" else VAR_NAMES
    if dependencies is not None and virtual_name not in dependencies:
        return frame
    frame = frame.copy()
    if axis == "obs":
        values = adata.obs_names.to_numpy()
        frame[OBS_NAMES] = values if row_positions is None else values[row_positions]
        frame.attrs[VIRTUAL_ATTR] = {OBS_NAMES}
    else:
        values = adata.var_names.to_numpy()
        frame[VAR_NAMES] = values if row_positions is None else values[row_positions]
        frame.attrs[VIRTUAL_ATTR] = {VAR_NAMES}
    return frame


def _matrix_virtual_values(
    adata: AnnData,
    *,
    axis: str,
    request: Any,
    schema: pd.DataFrame,
) -> tuple[str | None, np.ndarray | None]:
    dependencies = request_dependencies(request, schema)
    virtual_name = OBS_NAMES if axis == "obs" else VAR_NAMES
    if dependencies is not None and virtual_name not in dependencies:
        return None, None
    values = adata.obs_names.to_numpy() if axis == "obs" else adata.var_names.to_numpy()
    return virtual_name, values


def _add_matrix_request(
    planner: RequestPlanner,
    adata: AnnData,
    source: str,
    request: Any,
    *,
    mode: str = "expression",
    key: str | None = None,
    layer: str | None = None,
    axis: str,
    context: str,
    row_positions: Any = None,
) -> tuple[int, str]:
    adapter = source_adapter(adata, source, key=key, layer=layer)
    token = planner.add(
        adapter,
        request=request,
        mode=mode,
        row_positions=row_positions,
        context=context,
    )
    return token, axis


def _axis_positions(indexer: Any, size: int, *, axis: str) -> np.ndarray:
    """Normalize an axis indexer to zero-based integer positions."""
    base = np.arange(size, dtype=np.intp)
    if indexer is None:
        return base
    if isinstance(indexer, slice):
        return base[indexer]
    if isinstance(indexer, (int, np.integer)):
        values = np.asarray([indexer], dtype=np.intp)
    else:
        values = np.asarray(indexer)
        if values.ndim != 1:
            msg = f"{axis} positions must be a one-dimensional integer or boolean indexer"
            raise SelectionError(msg)
        if np.issubdtype(values.dtype, np.bool_):
            if len(values) != size:
                msg = f"{axis} boolean indexer has length {len(values)} for an axis of length {size}"
                raise SelectionError(msg)
            return base[values]
        if not np.issubdtype(values.dtype, np.integer):
            msg = f"{axis} positions must be integers, not axis labels"
            raise SelectionError(msg)
        values = values.astype(np.intp, copy=False)
    values = values.copy()
    values[values < 0] += size
    if ((values < 0) | (values >= size)).any():
        msg = f"{axis} position is outside an axis of length {size}"
        raise SelectionError(msg)
    return values


def _subset_positions(
    adata: AnnData,
    obs_positions: Any = None,
    var_positions: Any = None,
    *,
    copy: bool = True,
) -> AnnData:
    """Subset every aligned AnnData container using integer axis positions.

    A copied subset of a backed object is explicitly materialized in memory.
    ``copy=False`` is non-mutating and may return either an AnnData view or a
    materialized object, as allowed by the v0.3 ownership contract.
    """
    obs_idx = _axis_positions(obs_positions, adata.n_obs, axis="obs")
    var_idx = _axis_positions(var_positions, adata.n_vars, axis="var")
    selected = adata[obs_idx, var_idx]
    if not copy:
        return selected
    if adata.isbacked:
        return selected.to_memory(copy=True)
    return selected.copy()


def _subset(adata: AnnData, obs_idx: Any, var_idx: Any, *, copy: bool = True) -> AnnData:
    """Compatibility alias for internal callers; indexers are positions only."""
    return _subset_positions(adata, obs_idx, var_idx, copy=copy)


def _same_shape_target(adata: AnnData, *, inplace: bool, verb: str) -> AnnData:
    _ensure_not_backed(adata, verb)
    return adata if inplace else adata.copy()


def _inplace_reorder_axis(adata: AnnData, positions: np.ndarray, *, axis: str) -> None:
    """Reorder one complete axis while preserving the AnnData object's identity."""
    expected = adata.n_obs if axis == "obs" else adata.n_vars
    if len(positions) != expected or not np.array_equal(np.sort(positions), np.arange(expected)):
        msg = f"in-place {axis} reorder requires a complete positional permutation"
        raise SelectionError(msg)
    if np.array_equal(positions, np.arange(expected)):
        return
    if axis == "obs":
        adata._inplace_subset_obs(positions)
    else:
        adata._inplace_subset_var(positions)


def _ensure_not_backed(adata: AnnData, verb: str) -> None:
    if adata.isbacked:
        msg = f"{verb} cannot modify an AnnData object in backed mode; call .to_memory() first"
        raise AnnplyrError(msg)


def _obs_table(adata: AnnData) -> pd.DataFrame:
    return cast(pd.DataFrame, adata.obs)


def _var_table(adata: AnnData) -> pd.DataFrame:
    return cast(pd.DataFrame, adata.var)


def _axis(axis: str | int) -> str:
    if axis in {"obs", "observation", "observations", "0", 0}:
        return "obs"
    if axis in {"var", "variable", "variables", "1", 1}:
        return "var"
    msg = "axis must be 'obs' or 'var'"
    raise IncompatibleAxisError(msg)


def filter_adata(
    adata: AnnData,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    raw: Any = None,
    obs_names: Any = None,
    var_names: Any = None,
    obsm: Mapping[str, Any] | None = None,
    varm: Mapping[str, Any] | None = None,
    layer: str | None = None,
    copy: bool = True,
    max_matrix_values: int | None = None,
    _group_positions: tuple[np.ndarray, ...] | None = None,
    _group_axis: str | None = None,
) -> AnnData:
    obs_indices: list[pd.Index] = []
    var_indices: list[pd.Index] = []
    planner = RequestPlanner(max_matrix_values)
    planned: list[tuple[int, str, Any]] = []

    if x is not None:
        token, axis = _add_matrix_request(planner, adata, "x", x, layer=layer, axis="obs", context="filter x")
        planned.append((token, axis, x))
    if raw is not None:
        token, axis = _add_matrix_request(planner, adata, "raw", raw, axis="obs", context="filter raw")
        planned.append((token, axis, raw))
    for key, predicates in (obsm or {}).items():
        token, axis = _add_matrix_request(
            planner, adata, "obsm", predicates, key=key, axis="obs", context=f"filter obsm {key!r}"
        )
        planned.append((token, axis, predicates))
    for key, predicates in (varm or {}).items():
        token, axis = _add_matrix_request(
            planner, adata, "varm", predicates, key=key, axis="var", context=f"filter varm {key!r}"
        )
        planned.append((token, axis, predicates))
    projected = planner.execute()

    if obs is not None:
        obs_indices.append(
            _evaluate_filter_positions(obs_frame(adata), obs, _group_positions if _group_axis == "obs" else None)
        )
    if obs_names is not None:
        obs_indices.append(
            _evaluate_filter_positions(
                pd.DataFrame(
                    {"obs_names": adata.obs_names, OBS_NAMES: adata.obs_names},
                    index=pd.RangeIndex(adata.n_obs),
                ),
                obs_names,
                _group_positions if _group_axis == "obs" else None,
            )
        )
    for token, axis, predicates in planned:
        frame = _decorate_projected_frame(projected[token], adata, axis=axis, request=predicates)
        if axis == "obs":
            obs_indices.append(
                _evaluate_filter_positions(frame, predicates, _group_positions if _group_axis == "obs" else None)
            )
        else:
            var_indices.append(
                _evaluate_filter_positions(frame, predicates, _group_positions if _group_axis == "var" else None)
            )

    if var is not None:
        var_indices.append(
            _evaluate_filter_positions(var_frame(adata), var, _group_positions if _group_axis == "var" else None)
        )
    if var_names is not None:
        var_indices.append(
            _evaluate_filter_positions(
                pd.DataFrame(
                    {"var_names": adata.var_names, VAR_NAMES: adata.var_names},
                    index=pd.RangeIndex(adata.n_vars),
                ),
                var_names,
                _group_positions if _group_axis == "var" else None,
            )
        )
    obs_base = (
        pd.Index(np.concatenate(_group_positions))
        if _group_axis == "obs" and _group_positions
        else pd.RangeIndex(adata.n_obs)
    )
    var_base = (
        pd.Index(np.concatenate(_group_positions))
        if _group_axis == "var" and _group_positions
        else pd.RangeIndex(adata.n_vars)
    )
    obs_idx = intersect_ordered(obs_base, *obs_indices)
    var_idx = intersect_ordered(var_base, *var_indices)
    return _subset(adata, obs_idx, var_idx, copy=copy)


def _evaluate_filter_positions(
    frame: pd.DataFrame,
    predicates: Any,
    groups: tuple[np.ndarray, ...] | None,
) -> pd.Index:
    if groups is None:
        return evaluate_filter(frame, predicates)
    pieces: list[np.ndarray] = []
    for positions in groups:
        local = frame.iloc[positions, :].reset_index(drop=True)
        selected = evaluate_filter(local, predicates).to_numpy(dtype=np.intp)
        pieces.append(positions[selected])
    return pd.Index(np.concatenate(pieces) if pieces else np.empty(0, dtype=np.intp))


def select_adata(adata: AnnData, obs: Any = None, var: Any = None, x: Any = None, copy: bool = True) -> AnnData:
    obs_columns = (
        _selected_real_columns(evaluate_select(obs_frame(adata), obs).columns, _obs_table(adata).columns, source="obs")
        if obs is not None
        else _obs_table(adata).columns
    )
    var_columns = (
        _selected_real_columns(evaluate_select(var_frame(adata), var).columns, _var_table(adata).columns, source="var")
        if var is not None
        else _var_table(adata).columns
    )
    selected_var_names = (
        _selected_real_columns(evaluate_select(x_frame(adata), x).columns, adata.var_names, source="x")
        if x is not None
        else list(adata.var_names)
    )
    var_positions = _name_occurrence_positions(adata.var_names, selected_var_names, source="x")
    out = _subset_positions(adata, None, var_positions, copy=copy)
    out.obs = _obs_table(out).loc[:, list(obs_columns)].copy()
    out.var = _var_table(out).loc[:, list(var_columns)].copy()
    return out


def _name_occurrence_positions(available: pd.Index, selected: Sequence[str], *, source: str) -> np.ndarray:
    """Translate selected names to positions without indexing AnnData by labels."""
    positions_by_name: dict[str, list[int]] = {}
    for position, name in enumerate(available):
        positions_by_name.setdefault(str(name), []).append(position)
    offsets: dict[str, int] = {}
    positions: list[int] = []
    for name_value in selected:
        name = str(name_value)
        matches = positions_by_name.get(name, [])
        offset = offsets.get(name, 0)
        if offset >= len(matches):
            msg = f"Unknown or over-selected {source} column: {name!r}"
            raise UnknownColumnError(msg)
        positions.append(matches[offset])
        offsets[name] = offset + 1
    return np.asarray(positions, dtype=np.intp)


def _real_columns(selected: pd.Index, available: pd.Index) -> list[str]:
    return [column for column in selected if column in available]


def _selected_real_columns(selected: pd.Index, available: pd.Index, *, source: str) -> list[str]:
    columns = [str(column) for column in selected]
    real = [column for column in columns if column in available]
    computed = [column for column in columns if column not in available]
    if computed:
        msg = (
            f"{source} select can only keep existing AnnData-aligned columns; "
            f"computed or renamed column(s) are not supported: {', '.join(computed)}"
        )
        raise SelectionError(msg)
    return real


def rename_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, str] | None = None,
    var: Mapping[str, str] | None = None,
    x: Mapping[str, str] | None = None,
    inplace: bool = False,
) -> AnnData:
    _ensure_not_backed(adata, "rename")
    obs_mapping = _rename_mapping(_obs_table(adata).columns, obs, source="obs") if obs else None
    var_mapping = _rename_mapping(_var_table(adata).columns, var, source="var") if var else None
    x_names = _renamed_names(adata.var_names, x, source="x") if x else None
    out = _same_shape_target(adata, inplace=inplace, verb="rename")
    if obs_mapping:
        out.obs = _obs_table(out).rename(columns=obs_mapping).copy()
    if var_mapping:
        out.var = _var_table(out).rename(columns=var_mapping).copy()
    if x_names is not None:
        out.var_names = x_names
    return out


def rename_with_adata(
    adata: AnnData,
    func: Callable[[str], str],
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    inplace: bool = False,
) -> AnnData:
    _ensure_not_backed(adata, "rename_with")
    obs_mapping = (
        _rename_with_mapping(obs_frame(adata), obs, func, _obs_table(adata).columns, source="obs")
        if obs is not None
        else None
    )
    var_mapping = (
        _rename_with_mapping(var_frame(adata), var, func, _var_table(adata).columns, source="var")
        if var is not None
        else None
    )
    x_mapping = _rename_with_mapping(x_frame(adata), x, func, adata.var_names, source="x") if x is not None else None
    return rename_adata(adata, obs=obs_mapping, var=var_mapping, x=x_mapping, inplace=inplace)


def _rename_with_mapping(
    frame: pd.DataFrame,
    selector: Any,
    func: Callable[[str], str],
    available: pd.Index,
    *,
    source: str,
) -> dict[str, str]:
    selected = _real_columns(evaluate_select(frame, selector).columns, available)
    names = [func(old) for old in selected]
    _ensure_unique(names, source=source)
    return dict(zip(names, selected, strict=True))


def _rename_mapping(available: pd.Index, mapping: Mapping[str, str], *, source: str) -> dict[str, str]:
    duplicated_sources = (
        pd.Index(list(mapping.values()))[pd.Index(list(mapping.values())).duplicated()].unique().tolist()
    )
    if duplicated_sources:
        msg = f"Duplicate {source} source name(s): {', '.join(duplicated_sources)}"
        raise DuplicateNameError(msg)
    missing = [old for old in mapping.values() if old not in available]
    if missing:
        msg = f"Unknown {source} column(s): {', '.join(missing)}"
        raise UnknownColumnError(msg)
    proposed = [mapping.get(column, column) for column in available]
    _ensure_unique(proposed, source=source)
    return {old: new for new, old in mapping.items()}


def _renamed_names(available: pd.Index, mapping: Mapping[str, str], *, source: str) -> list[str]:
    missing = [old for old in mapping.values() if old not in available]
    if missing:
        msg = f"Unknown {source} name(s): {', '.join(missing)}"
        raise UnknownColumnError(msg)
    old_to_new = {old: new for new, old in mapping.items()}
    renamed = [old_to_new.get(name, name) for name in available]
    _ensure_unique(renamed, source=source)
    return renamed


def _ensure_unique(names: Sequence[str], *, source: str) -> None:
    duplicated = pd.Index(names)[pd.Index(names).duplicated()].unique().tolist()
    if duplicated:
        msg = f"Duplicate {source} name(s) after operation: {', '.join(duplicated)}"
        raise DuplicateNameError(msg)


def relocate_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    before: str | None = None,
    after: str | None = None,
    inplace: bool = False,
) -> AnnData:
    _ensure_not_backed(adata, "relocate")
    obs_order = (
        _relocated_order(
            _obs_table(adata).columns,
            _selected_columns(obs_frame(adata), obs, _obs_table(adata).columns),
            before,
            after,
        )
        if obs is not None
        else None
    )
    var_order = (
        _relocated_order(
            _var_table(adata).columns,
            _selected_columns(var_frame(adata), var, _var_table(adata).columns),
            before,
            after,
        )
        if var is not None
        else None
    )
    x_positions = _relocated_feature_positions(adata, x, before=before, after=after) if x is not None else None
    out = _same_shape_target(adata, inplace=inplace, verb="relocate")
    if obs_order is not None:
        out.obs = _obs_table(out).loc[:, obs_order].copy()
    if var_order is not None:
        out.var = _var_table(out).loc[:, var_order].copy()
    if x_positions is not None:
        _inplace_reorder_axis(out, x_positions, axis="var")
    return out


def _relocated_feature_positions(
    adata: AnnData,
    selector: Any,
    *,
    before: str | None,
    after: str | None,
) -> np.ndarray:
    selected_names = _selected_columns(x_frame(adata), selector, adata.var_names)
    selected = _name_occurrence_positions(adata.var_names, selected_names, source="x")
    selected_set = set(selected.tolist())
    remaining = [position for position in range(adata.n_vars) if position not in selected_set]

    def _anchor_position(name: str, *, kind: str) -> int:
        matches = [position for position in remaining if str(adata.var_names[position]) == name]
        if len(matches) > 1:
            msg = f"Relocate {kind} anchor {name!r} is ambiguous because the variable name is duplicated"
            raise DuplicateNameError(msg)
        if matches:
            return remaining.index(matches[0])
        if name in {str(adata.var_names[position]) for position in selected}:
            msg = f"Relocate anchor {name!r} is among the columns being moved; use a stationary column as anchor"
            raise SelectionError(msg)
        msg = f"Unknown relocate anchor: {name!r}"
        raise UnknownColumnError(msg)

    if before is not None and after is not None:
        msg = "relocate received both before and after anchors for the same source"
        raise SelectionError(msg)
    if before is not None:
        index = _anchor_position(before, kind="before")
    elif after is not None:
        index = _anchor_position(after, kind="after") + 1
    else:
        index = 0
    return np.asarray([*remaining[:index], *selected.tolist(), *remaining[index:]], dtype=np.intp)


def _selected_columns(frame: pd.DataFrame, selector: Any, available: pd.Index) -> list[str]:
    return _real_columns(evaluate_select(frame, selector).columns, available)


def _relocated_order(
    columns: pd.Index,
    selected: Sequence[str],
    before: str | None,
    after: str | None,
) -> list[str]:
    selected = [column for column in selected if column in columns]
    remaining = [column for column in columns if column not in selected]
    before_valid = before in remaining if before is not None else False
    after_valid = after in remaining if after is not None else False
    if before_valid and after_valid:
        msg = "relocate received both before and after anchors for the same source"
        raise SelectionError(msg)
    if before_valid:
        index = remaining.index(cast(str, before))
    elif after_valid:
        index = remaining.index(cast(str, after)) + 1
    elif before is None and after is None:
        index = 0
    elif before is not None:
        if before in columns:
            msg = f"Relocate anchor {before!r} is among the columns being moved; use a stationary column as anchor"
            raise SelectionError(msg)
        msg = f"Unknown relocate anchor: {before!r}"
        raise UnknownColumnError(msg)
    elif after is not None:
        if after in columns:
            msg = f"Relocate anchor {after!r} is among the columns being moved; use a stationary column as anchor"
            raise SelectionError(msg)
        msg = f"Unknown relocate anchor: {after!r}"
        raise UnknownColumnError(msg)
    else:
        index = 0
    return [*remaining[:index], *selected, *remaining[index:]]


def distinct_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    axis: str = "obs",
    keep_all: bool = False,
    copy: bool = True,
    max_matrix_values: int | None = None,
    _group_positions: tuple[np.ndarray, ...] | None = None,
    _group_axis: str | None = None,
) -> AnnData:
    axis = _axis(axis)
    if x is not None:
        if obs is not None or var is not None:
            msg = "distinct accepts one source at a time"
            raise UnknownSourceError(msg)
        if axis != "obs":
            msg = "x distinct is only defined on the obs axis"
            raise IncompatibleAxisError(msg)
        planner = RequestPlanner(max_matrix_values)
        token, _ = _add_matrix_request(planner, adata, "x", x, axis="obs", context="distinct x")
        frame = _decorate_projected_frame(planner.execute()[token], adata, axis="obs", request=x)
        selected = evaluate_select(frame, x)
        obs_idx = _distinct_positions(
            selected,
            _group_positions if _group_axis == "obs" else None,
        )
        return _subset(adata, obs_idx, slice(None), copy=copy)
    RequestPlanner(max_matrix_values).validate()
    if axis == "obs":
        frame, selector, available = _distinct_source(adata, obs=obs, x=x, axis=axis)
        selected = evaluate_select(frame, selector)
        obs_idx = _distinct_positions(
            selected,
            _group_positions if _group_axis == "obs" else None,
        )
        out = _subset(adata, obs_idx, slice(None), copy=copy)
        if not keep_all and obs is not None:
            out.obs = _obs_table(out).loc[:, _real_columns(selected.columns, available)].copy()
        return out

    frame, selector, available = _distinct_source(adata, var=var, axis=axis)
    selected = evaluate_select(frame, selector)
    var_idx = _distinct_positions(
        selected,
        _group_positions if _group_axis == "var" else None,
    )
    out = _subset(adata, slice(None), var_idx, copy=copy)
    if not keep_all and var is not None:
        out.var = _var_table(out).loc[:, _real_columns(selected.columns, available)].copy()
    return out


def _distinct_positions(frame: pd.DataFrame, groups: tuple[np.ndarray, ...] | None) -> np.ndarray:
    if groups is None:
        return frame.drop_duplicates(keep="first").index.to_numpy(dtype=np.intp)
    pieces = [
        frame.iloc[positions, :].drop_duplicates(keep="first").index.to_numpy(dtype=np.intp) for positions in groups
    ]
    return np.concatenate(pieces) if pieces else np.empty(0, dtype=np.intp)


def _distinct_source(
    adata: AnnData,
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    axis: str,
) -> tuple[pd.DataFrame, Any, pd.Index]:
    provided = [value is not None for value in (obs, var, x)]
    if sum(provided) > 1:
        msg = "distinct accepts one source at a time"
        raise UnknownSourceError(msg)
    if x is not None:
        if axis != "obs":
            msg = "x distinct is only defined on the obs axis"
            raise IncompatibleAxisError(msg)
        return x_frame(adata), x, adata.var_names
    if axis == "obs":
        return obs_frame(adata), obs, _obs_table(adata).columns
    return var_frame(adata), var, _var_table(adata).columns


def _sort_values_for_frame(frame: pd.DataFrame, by: Any) -> pd.Index:
    keys = _sort_keys(by)
    if not keys:
        return frame.index

    work = with_row_number(prepare_sparse_for_arithmetic(frame))
    sort_columns: list[str] = []
    ascending: list[bool] = []
    for i, (expr, is_descending) in enumerate(keys):
        name = f"__annplyr_sort_{i}__"
        selected = nw.from_native(work).select(_sort_expr(expr).alias(name)).to_native()
        work[name] = selected[name].to_numpy()
        sort_columns.append(name)
        ascending.append(not is_descending)

    sorted_frame = work.sort_values(sort_columns, ascending=ascending, kind="mergesort")
    return sorted_frame.index


def _sort_expr(expr: Any) -> Any:
    if isinstance(expr, str):
        return nw.col(expr)
    return to_narwhals(expr)


def _sort_keys(by: Any) -> list[tuple[Any, bool]]:
    if by is None:
        return []
    if isinstance(by, Desc):
        return [(by.expr, True)]
    if isinstance(by, (str, bytes)):
        return [(by, False)]
    if isinstance(by, Sequence):
        keys: list[tuple[Any, bool]] = []
        for item in by:
            if isinstance(item, Desc):
                keys.append((item.expr, True))
            else:
                keys.append((item, False))
        return keys
    return [(by, False)]


def arrange_adata(
    adata: AnnData,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    raw: Any = None,
    obsm: Mapping[str, Any] | None = None,
    varm: Mapping[str, Any] | None = None,
    layer: str | None = None,
    copy: bool = True,
    max_matrix_values: int | None = None,
    _group_positions: tuple[np.ndarray, ...] | None = None,
    _group_axis: str | None = None,
) -> AnnData:
    planner = RequestPlanner(max_matrix_values)
    obs_frames: list[tuple[pd.DataFrame | int, Any, str]] = []
    var_frames: list[tuple[pd.DataFrame | int, Any, str]] = []
    if obs is not None:
        obs_frames.append((obs_frame(adata), obs, "obs"))
    if x is not None:
        token, _ = _add_matrix_request(planner, adata, "x", x, layer=layer, axis="obs", context="arrange x")
        obs_frames.append((token, x, "obs"))
    if raw is not None:
        token, _ = _add_matrix_request(planner, adata, "raw", raw, axis="obs", context="arrange raw")
        obs_frames.append((token, raw, "obs"))
    for key, by in (obsm or {}).items():
        token, _ = _add_matrix_request(planner, adata, "obsm", by, key=key, axis="obs", context=f"arrange obsm {key!r}")
        obs_frames.append((token, by, "obs"))
    if var is not None:
        var_frames.append((var_frame(adata), var, "var"))
    for key, by in (varm or {}).items():
        token, _ = _add_matrix_request(planner, adata, "varm", by, key=key, axis="var", context=f"arrange varm {key!r}")
        var_frames.append((token, by, "var"))
    projected = planner.execute()

    def _resolved(items: list[tuple[pd.DataFrame | int, Any, str]]) -> list[tuple[pd.DataFrame, Any]]:
        resolved: list[tuple[pd.DataFrame, Any]] = []
        for frame_or_token, by, axis in items:
            frame = (
                _decorate_projected_frame(projected[frame_or_token], adata, axis=axis, request=by)
                if isinstance(frame_or_token, int)
                else frame_or_token
            )
            resolved.append((frame, by))
        return resolved

    resolved_obs = _resolved(obs_frames)
    resolved_var = _resolved(var_frames)
    obs_idx = _sort_values_for_group_frames(
        adata.n_obs,
        resolved_obs,
        _group_positions if _group_axis == "obs" else None,
    )
    var_idx = _sort_values_for_group_frames(
        adata.n_vars,
        resolved_var,
        _group_positions if _group_axis == "var" else None,
    )

    return _subset(adata, obs_idx, var_idx, copy=copy)


def _sort_values_for_group_frames(
    size: int,
    frame_by: list[tuple[pd.DataFrame, Any]],
    groups: tuple[np.ndarray, ...] | None,
) -> pd.Index:
    if groups is None:
        return _sort_values_for_frames(pd.RangeIndex(size), frame_by)
    pieces = [_sort_values_for_frames(pd.Index(positions), frame_by).to_numpy(dtype=np.intp) for positions in groups]
    return pd.Index(np.concatenate(pieces) if pieces else np.empty(0, dtype=np.intp))


def _sort_values_for_frames(base_index: pd.Index, frame_by_iter: Any) -> pd.Index:
    work = pd.DataFrame(index=base_index)
    sort_columns: list[str] = []
    ascending: list[bool] = []
    for frame, by in frame_by_iter:
        frame = frame.loc[base_index]
        for expr, is_descending in _sort_keys(by):
            name = f"__annplyr_sort_{len(sort_columns)}__"
            work_frame = with_row_number(prepare_sparse_for_arithmetic(frame))
            selected = nw.from_native(work_frame).select(_sort_expr(expr).alias(name)).to_native()
            work[name] = selected[name].to_numpy()
            sort_columns.append(name)
            ascending.append(not is_descending)
    if not sort_columns:
        return base_index
    return work.sort_values(sort_columns, ascending=ascending, kind="mergesort").index


def _obs_sort_frames(
    adata: AnnData,
    *,
    obs: Any = None,
    x: Any = None,
    raw: Any = None,
    obsm: Mapping[str, Any] | None = None,
    layer: str | None = None,
):
    if obs is not None:
        yield obs_frame(adata), obs
    if x is not None:
        yield x_frame(adata, layer=layer), x
    if raw is not None:
        yield raw_frame(adata), raw
    for key, by in (obsm or {}).items():
        yield obsm_frame(adata, key), by


def _var_sort_frames(adata: AnnData, *, var: Any = None, varm: Mapping[str, Any] | None = None):
    if var is not None:
        yield var_frame(adata), var
    for key, by in (varm or {}).items():
        yield varm_frame(adata, key), by


def slice_adata(adata: AnnData, *indices: Any, axis: str = "obs", copy: bool = True) -> AnnData:
    axis = _axis(axis)
    selector = _slice_selector(indices)
    if axis == "obs":
        return _subset(adata, selector, slice(None), copy=copy)
    return _subset(adata, slice(None), selector, copy=copy)


def _slice_selector(indices: tuple[Any, ...]) -> Any:
    if not indices:
        return slice(None)
    if len(indices) == 1:
        only = indices[0]
        if isinstance(only, slice):
            return only
        if isinstance(only, Sequence) and not isinstance(only, (str, bytes)):
            return list(only)
    return list(indices)


def slice_head_adata(adata: AnnData, n: int = 5, *, axis: str = "obs", copy: bool = True) -> AnnData:
    _validate_slice_n(n)
    return slice_adata(adata, slice(0, n), axis=axis, copy=copy)


def slice_tail_adata(adata: AnnData, n: int = 5, *, axis: str = "obs", copy: bool = True) -> AnnData:
    _validate_slice_n(n)
    if n == 0:
        return slice_adata(adata, slice(0, 0), axis=axis, copy=copy)
    return slice_adata(adata, slice(-n, None), axis=axis, copy=copy)


def _validate_slice_n(n: int) -> None:
    if n < 0:
        msg = "slice n must be non-negative"
        raise SelectionError(msg)


def slice_min_adata(adata: AnnData, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = True) -> AnnData:
    axis = _axis(axis)
    _validate_slice_n(n)
    frame = obs_frame(adata) if axis == "obs" else var_frame(adata)
    positions = _sort_values_for_frame(frame, by).to_numpy(dtype=np.intp)[:n]
    return (
        _subset_positions(adata, positions, None, copy=copy)
        if axis == "obs"
        else _subset_positions(adata, None, positions, copy=copy)
    )


def slice_max_adata(adata: AnnData, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = True) -> AnnData:
    axis = _axis(axis)
    _validate_slice_n(n)
    by_desc = _desc_order_by(by)
    frame = obs_frame(adata) if axis == "obs" else var_frame(adata)
    positions = _sort_values_for_frame(frame, by_desc).to_numpy(dtype=np.intp)[:n]
    return (
        _subset_positions(adata, positions, None, copy=copy)
        if axis == "obs"
        else _subset_positions(adata, None, positions, copy=copy)
    )


def _desc_order_by(by: Any) -> Any:
    if isinstance(by, Sequence) and not isinstance(by, (str, bytes)):
        return [item if isinstance(item, Desc) else Desc(item) for item in by]
    return by if isinstance(by, Desc) else Desc(by)


def slice_sample_adata(
    adata: AnnData,
    n: int | None = None,
    *,
    prop: float | None = None,
    replace: bool = False,
    random_state: int | None = None,
    axis: str = "obs",
    copy: bool = True,
) -> AnnData:
    axis = _axis(axis)
    size = adata.n_obs if axis == "obs" else adata.n_vars
    if n is not None and prop is not None:
        msg = "slice_sample accepts n and prop as mutually exclusive arguments"
        raise SelectionError(msg)
    if n is not None and n < 0:
        msg = "slice_sample n must be non-negative"
        raise SelectionError(msg)
    if prop is not None and prop < 0:
        msg = "slice_sample prop must be non-negative"
        raise SelectionError(msg)
    if n is None:
        n = int(round(size * prop)) if prop is not None else min(size, 1)
    if not replace and n > size:
        msg = "slice_sample n cannot be larger than the axis size unless replace=True"
        raise SelectionError(msg)
    rng = np.random.default_rng(random_state)
    selected = rng.choice(size, size=n, replace=replace)
    return slice_adata(adata, selected.tolist(), axis=axis, copy=copy)


def mutate_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    var: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    raw: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
    layer: str | None = None,
    inplace: bool = False,
    max_matrix_values: int | None = None,
    _group_positions: tuple[np.ndarray, ...] | None = None,
    _group_axis: str | None = None,
) -> AnnData:
    _ensure_not_backed(adata, "mutate")
    obs_values, var_values = _evaluate_mutation_sources(
        adata,
        obs=obs,
        var=var,
        x=x,
        raw=raw,
        obsm=obsm,
        varm=varm,
        layer=layer,
        max_matrix_values=max_matrix_values,
        group_positions=_group_positions,
        group_axis=_group_axis,
    )
    out = _same_shape_target(adata, inplace=inplace, verb="mutate")
    for values in obs_values:
        for column in values.columns:
            _assign_positional(_obs_table(out), str(column), values[column])
    for values in var_values:
        for column in values.columns:
            _assign_positional(_var_table(out), str(column), values[column])
    return out


def _assign_positional(table: pd.DataFrame, name: str, values: pd.Series) -> None:
    """Assign a Series by position while retaining its pandas extension dtype."""
    series = values.reset_index(drop=True)
    if name in table and table[name].reset_index(drop=True).equals(series):
        return
    table[name] = pd.Series(series.array, index=table.index, name=name)


def transmute_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    var: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    raw: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
    layer: str | None = None,
    max_matrix_values: int | None = None,
    _group_positions: tuple[np.ndarray, ...] | None = None,
    _group_axis: str | None = None,
) -> AnnData:
    obs_values, var_values = _evaluate_mutation_sources(
        adata,
        obs=obs,
        var=var,
        x=x,
        raw=raw,
        obsm=obsm,
        varm=varm,
        layer=layer,
        max_matrix_values=max_matrix_values,
        group_positions=_group_positions,
        group_axis=_group_axis,
    )
    obs_columns = [str(column) for values in obs_values for column in values.columns]
    var_columns = [str(column) for values in var_values for column in values.columns]
    out = adata.to_memory(copy=True) if adata.isbacked else adata.copy()
    for values in obs_values:
        for column in values.columns:
            _assign_positional(_obs_table(out), str(column), values[column])
    for values in var_values:
        for column in values.columns:
            _assign_positional(_var_table(out), str(column), values[column])
    if obs_columns:
        out.obs = _obs_table(out).loc[:, obs_columns].copy()
    if var_columns:
        out.var = _var_table(out).loc[:, var_columns].copy()
    return out


def _evaluate_mutation_sources(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None,
    var: Mapping[str, Any] | None,
    x: Mapping[str, Any] | None,
    raw: Mapping[str, Any] | None,
    obsm: Mapping[str, Mapping[str, Any]] | None,
    varm: Mapping[str, Mapping[str, Any]] | None,
    layer: str | None,
    max_matrix_values: int | None,
    group_positions: tuple[np.ndarray, ...] | None = None,
    group_axis: str | None = None,
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    planner = RequestPlanner(max_matrix_values)
    obs_sources: list[tuple[pd.DataFrame | int, Mapping[str, Any], str]] = []
    var_sources: list[tuple[pd.DataFrame | int, Mapping[str, Any], str]] = []
    if obs:
        obs_sources.append((obs_frame(adata), obs, "obs"))
    if x:
        token, _ = _add_matrix_request(planner, adata, "x", x, layer=layer, axis="obs", context="mutate x")
        obs_sources.append((token, x, "obs"))
    if raw:
        token, _ = _add_matrix_request(planner, adata, "raw", raw, axis="obs", context="mutate raw")
        obs_sources.append((token, raw, "obs"))
    for key, assignments in (obsm or {}).items():
        token, _ = _add_matrix_request(
            planner, adata, "obsm", assignments, key=key, axis="obs", context=f"mutate obsm {key!r}"
        )
        obs_sources.append((token, assignments, "obs"))
    if var:
        var_sources.append((var_frame(adata), var, "var"))
    for key, assignments in (varm or {}).items():
        token, _ = _add_matrix_request(
            planner, adata, "varm", assignments, key=key, axis="var", context=f"mutate varm {key!r}"
        )
        var_sources.append((token, assignments, "var"))
    projected = planner.execute()

    def _evaluate(sources: list[tuple[pd.DataFrame | int, Mapping[str, Any], str]]) -> list[pd.DataFrame]:
        values: list[pd.DataFrame] = []
        for frame_or_token, assignments, axis in sources:
            frame = (
                _decorate_projected_frame(projected[frame_or_token], adata, axis=axis, request=assignments)
                if isinstance(frame_or_token, int)
                else frame_or_token
            )
            if group_positions is not None and group_axis == axis:
                pieces = [
                    (positions, evaluate_assignments(frame.iloc[positions, :].reset_index(drop=True), assignments))
                    for positions in group_positions
                ]
                values.append(_scatter_group_assignments(pieces, len(frame)))
            else:
                values.append(evaluate_assignments(frame, assignments))
        return values

    return _evaluate(obs_sources), _evaluate(var_sources)


def _scatter_group_assignments(
    pieces: list[tuple[np.ndarray, pd.DataFrame]],
    size: int,
) -> pd.DataFrame:
    if not pieces:
        return pd.DataFrame(index=pd.RangeIndex(size))
    columns = [str(column) for column in pieces[0][1].columns]
    output = pd.DataFrame(index=pd.RangeIndex(size))
    all_positions = np.concatenate([positions for positions, _ in pieces])
    for column in columns:
        concatenated = pd.concat([values[column].reset_index(drop=True) for _, values in pieces], ignore_index=True)
        concatenated.index = all_positions
        ordered = concatenated.sort_index(kind="stable")
        output[column] = pd.Series(ordered.array, index=pd.RangeIndex(size), name=column)
    return output


def summarize_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    var: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    raw: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
    by: Any = None,
    layer: str | None = None,
    max_matrix_values: int | None = None,
    _group_plan: GroupPlan | None = None,
) -> pd.DataFrame:
    obs_axis_requested = any(source is not None for source in (obs, x, raw, obsm))
    var_axis_requested = any(source is not None for source in (var, varm))
    if obs_axis_requested and var_axis_requested:
        msg = "summarize accepts obs-axis or var-axis sources, not both at once"
        raise IncompatibleAxisError(msg)

    planner = RequestPlanner(max_matrix_values)
    descriptors: list[tuple[pd.DataFrame | int, Mapping[str, Any], str]] = []
    if var_axis_requested:
        if var:
            descriptors.append((var_frame(adata), var, "var"))
        for key, assignments in (varm or {}).items():
            token, _ = _add_matrix_request(
                planner,
                adata,
                "varm",
                assignments,
                key=key,
                axis="var",
                context=f"summarize varm {key!r}",
            )
            descriptors.append((token, assignments, "var"))
        by_source = var_frame(adata)
    else:
        if obs:
            descriptors.append((obs_frame(adata), obs, "obs"))
        if x:
            token, _ = _add_matrix_request(planner, adata, "x", x, layer=layer, axis="obs", context="summarize x")
            descriptors.append((token, x, "obs"))
        if raw:
            token, _ = _add_matrix_request(planner, adata, "raw", raw, axis="obs", context="summarize raw")
            descriptors.append((token, raw, "obs"))
        for key, assignments in (obsm or {}).items():
            token, _ = _add_matrix_request(
                planner,
                adata,
                "obsm",
                assignments,
                key=key,
                axis="obs",
                context=f"summarize obsm {key!r}",
            )
            descriptors.append((token, assignments, "obs"))
        by_source = obs_frame(adata)
    planner.validate()
    groups = SummaryGroupPlan.build(by_source, by, grouped=_group_plan)
    pieces: list[pd.DataFrame] = []
    for frame_or_token, assignments, axis in descriptors:
        if not isinstance(frame_or_token, int):
            positions = np.arange(len(frame_or_token), dtype=np.intp)
            work = groups.add_to_frame(frame_or_token, positions)
            pieces.append(summarize_frame(work, assignments=assignments, by=groups.internal_columns))
            continue

        request = planner.requests[frame_or_token]
        schema = request.adapter.schema
        virtual_name, virtual_values = _matrix_virtual_values(
            adata,
            axis=axis,
            request=assignments,
            schema=schema,
        )
        planning_schema = schema.copy()
        if virtual_name is not None and virtual_values is not None:
            planning_schema[virtual_name] = pd.Series([], dtype=virtual_values.dtype)
        reductions = ReductionPlan.resolve(planning_schema, assignments)
        if reductions is not None:
            chunks = planner.chunk_plan(
                frame_or_token,
                target_values=DEFAULT_REDUCTION_CHUNK_VALUES,
            )
            pieces.append(
                summarize_chunked(
                    chunks,
                    reductions,
                    groups,
                    virtual_name=virtual_name,
                    virtual_values=virtual_values,
                )
            )
            continue

        frame = _decorate_projected_frame(request.read(), adata, axis=axis, request=assignments)
        positions = request.row_positions
        work = groups.add_to_frame(frame, positions)
        pieces.append(summarize_frame(work, assignments=assignments, by=groups.internal_columns))

    if not pieces:
        return count_frame(by_source, by=by)
    result = _merge_summary_pieces(pieces, list(groups.internal_columns))
    return groups.restore_keys(result)


def _merge_summary_pieces(pieces: list[pd.DataFrame], by_columns: list[str]) -> pd.DataFrame:
    result = pieces[0].reset_index(drop=True)
    for piece in pieces[1:]:
        piece = piece.reset_index(drop=True)
        if by_columns:
            result = result.merge(piece, on=by_columns, how="outer", sort=False)
        else:
            result = pd.concat([result, piece], axis=1)
    return result


def summarize_frame(frame: pd.DataFrame, assignments: Mapping[str, Any] | None, by: Any = None) -> pd.DataFrame:
    work, by_columns = prepare_by_frame(frame, by)
    work = prepare_sparse_for_arithmetic(work)
    assignments = expand_assignments(work, assignments)
    exprs = [
        to_narwhals(expr.alias(name)) if hasattr(expr, "alias") else nw.col(expr).alias(name)
        for name, expr in assignments.items()
    ]
    if by_columns:
        return nw.from_native(work).group_by(*by_columns).agg(*exprs).to_native()
    if exprs:
        return nw.from_native(with_row_number(work)).select(*exprs).to_native()
    return pd.DataFrame(index=[0])


def prepare_by_frame(frame: pd.DataFrame, by: Any) -> tuple[pd.DataFrame, list[str]]:
    if by is None:
        return frame.copy(), []
    by_frame = evaluate_select(frame, by)
    by_columns = list(by_frame.columns)
    work = frame.copy()
    for column in by_columns:
        if column in work.columns and not work[column].reset_index(drop=True).equals(
            by_frame[column].reset_index(drop=True)
        ):
            msg = f"Grouping column {column!r} would overwrite an existing source column"
            raise DuplicateNameError(msg)
        work[column] = by_frame[column].to_numpy()
    return work, by_columns


def count_adata(
    adata: AnnData,
    by: Any = None,
    *,
    wt: Any = None,
    sort: bool = False,
    axis: str = "obs",
    name: str = "n",
) -> pd.DataFrame:
    frame = obs_frame(adata) if _axis(axis) == "obs" else var_frame(adata)
    return count_frame(frame, by=by, wt=wt, sort=sort, name=name)


def tally_adata(
    adata: AnnData,
    by: Any = None,
    *,
    wt: Any = None,
    sort: bool = False,
    axis: str = "obs",
    name: str = "n",
) -> pd.DataFrame:
    return count_adata(adata, by=by, wt=wt, sort=sort, axis=axis, name=name)


def add_count_adata(
    adata: AnnData,
    by: Any = None,
    *,
    wt: Any = None,
    sort: bool = False,
    axis: str = "obs",
    name: str = "n",
    inplace: bool = False,
) -> AnnData:
    _ensure_not_backed(adata, "add_count")
    axis = _axis(axis)
    frame = obs_frame(adata) if axis == "obs" else var_frame(adata)
    values = _count_values(frame, by=by, wt=wt)
    positions = (
        values.sort_values(ascending=False, kind="mergesort").index.to_numpy(dtype=np.intp)
        if sort
        else np.arange(len(values), dtype=np.intp)
    )
    out = _same_shape_target(adata, inplace=inplace, verb="add_count")
    if axis == "obs":
        _assign_positional(_obs_table(out), name, values)
    else:
        _assign_positional(_var_table(out), name, values)
    if sort:
        _inplace_reorder_axis(out, positions, axis=axis)
    return out


def add_tally_adata(
    adata: AnnData,
    *,
    wt: Any = None,
    sort: bool = False,
    axis: str = "obs",
    name: str = "n",
    inplace: bool = False,
) -> AnnData:
    return add_count_adata(adata, wt=wt, sort=sort, axis=axis, name=name, inplace=inplace)


def _count_values(frame: pd.DataFrame, by: Any = None, wt: Any = None) -> pd.Series:
    by_frame = evaluate_select(frame, by) if by is not None else pd.DataFrame(index=frame.index)
    by_columns = list(by_frame.columns)
    if wt is not None:
        weights = evaluate_assignments(frame, {"__annplyr_wt__": wt})["__annplyr_wt__"]
    else:
        weights = pd.Series(1, index=frame.index)
    if not by_columns:
        return pd.Series(weights.sum(), index=frame.index)
    work = by_frame.copy()
    work["__annplyr_count_row__"] = np.arange(len(work))
    work["__annplyr_wt__"] = weights.to_numpy()
    if wt is None:
        return work.groupby(by_columns, sort=False, dropna=False)["__annplyr_count_row__"].transform("size")
    return work.groupby(by_columns, sort=False, dropna=False)["__annplyr_wt__"].transform("sum")


def count_frame(
    frame: pd.DataFrame,
    by: Any = None,
    *,
    wt: Any = None,
    sort: bool = False,
    name: str = "n",
) -> pd.DataFrame:
    work, by_columns = prepare_by_frame(frame, by)
    if wt is not None:
        work["__annplyr_wt__"] = evaluate_assignments(frame, {"__annplyr_wt__": wt})["__annplyr_wt__"].to_numpy()
    if not by_columns:
        value = work["__annplyr_wt__"].sum() if wt is not None else len(frame)
        return pd.DataFrame({name: [value]})
    grouped = work.groupby(by_columns, sort=False, observed=True, dropna=False)
    result = (
        grouped["__annplyr_wt__"].sum().reset_index(name=name)
        if wt is not None
        else grouped.size().reset_index(name=name)
    )
    if sort:
        result = result.sort_values(name, ascending=False, kind="mergesort").reset_index(drop=True)
    return result


def left_join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None = None,
    axis: str = "obs",
    relationship: str = "many-to-one",
    multiple: str = "error",
    unmatched: str = "drop",
    na_matches: str = "na",
    suffixes: tuple[str, str] = ("", "_right"),
    copy: bool = True,
) -> AnnData:
    return _join_adata(
        adata,
        other,
        by=by,
        axis=axis,
        how="left",
        relationship=relationship,
        multiple=multiple,
        unmatched=unmatched,
        na_matches=na_matches,
        suffixes=suffixes,
        copy=copy,
    )


def inner_join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None = None,
    axis: str = "obs",
    relationship: str = "many-to-one",
    multiple: str = "error",
    unmatched: str = "drop",
    na_matches: str = "na",
    suffixes: tuple[str, str] = ("", "_right"),
    copy: bool = True,
) -> AnnData:
    return _join_adata(
        adata,
        other,
        by=by,
        axis=axis,
        how="inner",
        relationship=relationship,
        multiple=multiple,
        unmatched=unmatched,
        na_matches=na_matches,
        suffixes=suffixes,
        copy=copy,
    )


def right_join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None = None,
    axis: str = "obs",
    relationship: str = "many-to-one",
    multiple: str = "error",
    unmatched: str = "error",
    na_matches: str = "na",
    suffixes: tuple[str, str] = ("", "_right"),
    copy: bool = True,
) -> AnnData:
    return _join_adata(
        adata,
        other,
        by=by,
        axis=axis,
        how="right",
        relationship=relationship,
        multiple=multiple,
        unmatched=unmatched,
        na_matches=na_matches,
        suffixes=suffixes,
        copy=copy,
    )


def full_join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None = None,
    axis: str = "obs",
    relationship: str = "many-to-one",
    multiple: str = "error",
    unmatched: str = "error",
    na_matches: str = "na",
    suffixes: tuple[str, str] = ("", "_right"),
    copy: bool = True,
) -> AnnData:
    return _join_adata(
        adata,
        other,
        by=by,
        axis=axis,
        how="outer",
        relationship=relationship,
        multiple=multiple,
        unmatched=unmatched,
        na_matches=na_matches,
        suffixes=suffixes,
        copy=copy,
    )


def semi_join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None = None,
    axis: str = "obs",
    na_matches: str = "na",
    copy: bool = True,
) -> AnnData:
    mask = _join_filter_mask(adata, other, by=by, axis=axis, keep_matches=True, na_matches=na_matches)
    return (
        _subset(adata, mask.to_numpy(), slice(None), copy=copy)
        if _axis(axis) == "obs"
        else _subset(adata, slice(None), mask.to_numpy(), copy=copy)
    )


def anti_join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None = None,
    axis: str = "obs",
    na_matches: str = "na",
    copy: bool = True,
) -> AnnData:
    mask = _join_filter_mask(adata, other, by=by, axis=axis, keep_matches=False, na_matches=na_matches)
    return (
        _subset(adata, mask.to_numpy(), slice(None), copy=copy)
        if _axis(axis) == "obs"
        else _subset(adata, slice(None), mask.to_numpy(), copy=copy)
    )


def _join_adata(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None,
    axis: str,
    how: str,
    relationship: str,
    multiple: str,
    unmatched: str,
    na_matches: str,
    suffixes: tuple[str, str],
    copy: bool,
) -> AnnData:
    axis = _axis(axis)
    _validate_join_unmatched(unmatched)
    _validate_join_na_matches(na_matches)
    left = _axis_table(adata, axis).reset_index(drop=True).copy()
    right = _coerce_join_frame(other)
    by_columns = _join_by_columns(left, right, by)
    _validate_left_join_relationship(left, by_columns, relationship=relationship)
    right = _prepare_join_right(right, by_columns, relationship=relationship, multiple=multiple)
    if (
        unmatched == "error"
        and not _join_filter_mask(
            adata,
            right,
            by=by_columns,
            axis=axis,
            keep_matches=True,
            na_matches=na_matches,
        ).all()
    ):
        msg = f"{how}_join has unmatched axis records"
        raise JoinRelationshipError(msg)
    left_key = "__annplyr_axis_position__"
    left[left_key] = np.arange(len(left), dtype=np.intp)
    right_order = "__annplyr_right_order__"
    if how in {"right", "outer"}:
        right = right.copy()
        right[right_order] = np.arange(len(right), dtype=np.intp)
    merge_how = "outer" if how == "outer" else how
    joined = _merge_join_frames(
        left,
        right,
        by_columns=by_columns,
        how=merge_how,
        suffixes=suffixes,
        na_matches=na_matches,
    )
    if how in {"right", "outer"} and joined[left_key].isna().any():
        msg = f"{how}_join would add axis records, which AnnData-aligned joins do not support"
        raise JoinRelationshipError(msg)
    if unmatched == "error" and joined["_merge"].eq("left_only").any():
        msg = f"{how}_join has unmatched axis records"
        raise JoinRelationshipError(msg)
    joined = joined[joined[left_key].notna()].copy()
    if joined[left_key].duplicated().any():
        msg = f"{how}_join would duplicate AnnData axis records"
        raise JoinRelationshipError(msg)
    if how == "right":
        joined = joined.sort_values([right_order, left_key], kind="mergesort")
    elif how == "outer":
        joined = joined.sort_values(left_key, kind="mergesort")
    positions = joined[left_key].astype(np.intp).to_numpy()
    table = joined.drop(columns=[left_key, right_order, "_merge"], errors="ignore")
    labels = adata.obs_names.take(positions) if axis == "obs" else adata.var_names.take(positions)
    table.index = labels
    _ensure_unique([str(column) for column in table.columns], source=f"{axis} join")
    out = (
        _subset(adata, positions, slice(None), copy=copy)
        if axis == "obs"
        else _subset(adata, slice(None), positions, copy=copy)
    )
    if axis == "obs":
        out.obs = table.copy()
    else:
        out.var = table.copy()
    return out


def _axis_table(adata: AnnData, axis: str) -> pd.DataFrame:
    return _obs_table(adata) if axis == "obs" else _var_table(adata)


def _validate_join_unmatched(unmatched: str) -> None:
    if unmatched not in {"drop", "error"}:
        msg = "unmatched must be 'drop' or 'error'"
        raise JoinRelationshipError(msg)


def _validate_join_na_matches(na_matches: str) -> None:
    if na_matches not in {"na", "never"}:
        msg = "na_matches must be 'na' or 'never'"
        raise JoinRelationshipError(msg)


def _coerce_join_frame(other: pd.DataFrame | Mapping[str, Any]) -> pd.DataFrame:
    return other.copy() if isinstance(other, pd.DataFrame) else pd.DataFrame(other)


def _join_by_columns(left: pd.DataFrame, right: pd.DataFrame, by: str | Sequence[str] | None) -> list[str]:
    if by is None:
        columns = [column for column in left.columns if column in right.columns]
    elif isinstance(by, str):
        columns = [by]
    else:
        columns = list(by)
    if not columns:
        msg = "join requires at least one shared key column"
        raise JoinRelationshipError(msg)
    missing_left = [column for column in columns if column not in left.columns]
    missing_right = [column for column in columns if column not in right.columns]
    if missing_left or missing_right:
        missing = [*(f"left.{column}" for column in missing_left), *(f"right.{column}" for column in missing_right)]
        msg = f"Unknown join key column(s): {', '.join(missing)}"
        raise UnknownColumnError(msg)
    return [str(column) for column in columns]


def _prepare_join_right(
    right: pd.DataFrame,
    by_columns: Sequence[str],
    *,
    relationship: str,
    multiple: str,
) -> pd.DataFrame:
    if multiple not in {"error", "first", "all"}:
        msg = "multiple must be 'error', 'first', or 'all'"
        raise JoinRelationshipError(msg)
    duplicated_right = right.duplicated(list(by_columns), keep=False)
    if multiple == "first":
        right = right.drop_duplicates(list(by_columns), keep="first")
    elif multiple == "error" and duplicated_right.any():
        msg = "join found multiple right-hand matches for at least one key"
        raise JoinRelationshipError(msg)

    if relationship in {"many-to-one", "one-to-one"} and right.duplicated(list(by_columns), keep=False).any():
        msg = f"join relationship {relationship!r} requires unique right-hand keys"
        raise JoinRelationshipError(msg)
    if relationship not in {"many-to-one", "one-to-one", "one-to-many", "many-to-many"}:
        msg = "relationship must be 'many-to-one', 'one-to-one', 'one-to-many', or 'many-to-many'"
        raise JoinRelationshipError(msg)
    return right


def _validate_left_join_relationship(left: pd.DataFrame, by_columns: Sequence[str], *, relationship: str) -> None:
    if relationship in {"one-to-one", "one-to-many"} and left.duplicated(list(by_columns), keep=False).any():
        msg = f"join relationship {relationship!r} requires unique left-hand keys"
        raise JoinRelationshipError(msg)


def _merge_join_frames(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    by_columns: Sequence[str],
    how: str,
    suffixes: tuple[str, str],
    na_matches: str,
) -> pd.DataFrame:
    if na_matches != "never":
        return left.merge(
            right,
            how=cast(Any, how),
            on=list(by_columns),
            sort=False,
            suffixes=suffixes,
            indicator=True,
        )

    left_work = left.copy()
    right_work = right.copy()
    temp_keys = [f"__annplyr_join_key_{i}__" for i, _ in enumerate(by_columns)]
    for temp_key, column in zip(temp_keys, by_columns, strict=True):
        left_work[temp_key] = left_work[column].astype(object)
        right_work[temp_key] = right_work[column].astype(object)
        left_work.loc[left_work[column].isna(), temp_key] = f"__annplyr_left_na_{column}__"
        right_work.loc[right_work[column].isna(), temp_key] = f"__annplyr_right_na_{column}__"
    right_work = right_work.drop(columns=list(by_columns))
    joined = left_work.merge(
        right_work,
        how=cast(Any, how),
        on=temp_keys,
        sort=False,
        suffixes=suffixes,
        indicator=True,
    )
    return joined.drop(columns=temp_keys)


def _join_filter_mask(
    adata: AnnData,
    other: pd.DataFrame | Mapping[str, Any],
    *,
    by: str | Sequence[str] | None,
    axis: str,
    keep_matches: bool,
    na_matches: str,
) -> pd.Series:
    axis = _axis(axis)
    _validate_join_na_matches(na_matches)
    left_table = _axis_table(adata, axis)
    right = _coerce_join_frame(other)
    by_columns = _join_by_columns(left_table, right, by)
    left = left_table.loc[:, by_columns].reset_index(drop=True)
    right_keys = right.loc[:, by_columns]
    if len(by_columns) == 1:
        right_keys = right_keys.drop_duplicates()
        column = by_columns[0]
        matches = left[column].isin(right_keys[column].dropna())
        if na_matches == "na" and right_keys[column].isna().any():
            matches |= left[column].isna()
    else:
        right_keys = right_keys.groupby(by_columns, sort=False, observed=True, dropna=False).head(1)
        merged = left.merge(right_keys, how="left", on=by_columns, sort=False, indicator=True)
        matches = merged["_merge"].eq("both")
        if na_matches == "never":
            matches &= ~left.isna().any(axis=1).to_numpy()
    values = matches.to_numpy() if keep_matches else (~matches).to_numpy()
    return pd.Series(values, index=pd.RangeIndex(len(left)))


def pull_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    raw: Any = None,
    obsm: Mapping[str, Any] | None = None,
    varm: Mapping[str, Any] | None = None,
    obsp: Mapping[str, Any] | None = None,
    varp: Mapping[str, Any] | None = None,
    uns: Mapping[str, Any] | None = None,
    layer: str | None = None,
    max_matrix_values: int | None = None,
) -> pd.Series:
    provided = [value is not None for value in [obs, var, x, raw, obsm, varm, obsp, varp, uns]]
    if sum(provided) != 1:
        msg = "pull requires exactly one source"
        raise UnknownSourceError(msg)
    if obs is not None:
        RequestPlanner(max_matrix_values).validate()
        return _first_series(evaluate_select(obs_frame(adata), obs))
    if var is not None:
        RequestPlanner(max_matrix_values).validate()
        return _first_series(evaluate_select(var_frame(adata), var))
    matrix_spec: tuple[str, str | None, Any, str] | None = None
    if x is not None:
        matrix_spec = ("x", None, x, "obs")
    elif raw is not None:
        matrix_spec = ("raw", None, raw, "obs")
    elif obsm is not None:
        key, selector = _single_pull_mapping_item(obsm, source="obsm")
        matrix_spec = ("obsm", key, selector, "obs")
    elif varm is not None:
        key, selector = _single_pull_mapping_item(varm, source="varm")
        matrix_spec = ("varm", key, selector, "var")
    elif obsp is not None:
        key, selector = _single_pull_mapping_item(obsp, source="obsp")
        matrix_spec = ("obsp", key, selector, "obs")
    elif varp is not None:
        key, selector = _single_pull_mapping_item(varp, source="varp")
        matrix_spec = ("varp", key, selector, "var")
    if matrix_spec is not None:
        source, request_key, selector, axis = matrix_spec
        planner = RequestPlanner(max_matrix_values)
        token, _ = _add_matrix_request(
            planner,
            adata,
            source,
            selector,
            mode="selection",
            key=request_key,
            layer=layer if source == "x" else None,
            axis=axis,
            context=f"pull {source}",
        )
        frame = _decorate_projected_frame(planner.execute()[token], adata, axis=axis, request=selector)
        result = _first_series(evaluate_select(frame, selector))
        result.index = adata.obs_names if axis == "obs" else adata.var_names
        return result
    RequestPlanner(max_matrix_values).validate()
    key, selector = _single_pull_mapping_item(cast(Mapping[str, Any], uns), source="uns")
    return _first_series(evaluate_select(source_frame(adata, "uns", key=key), selector))


def _single_pull_mapping_item(mapping: Mapping[str, Any], *, source: str) -> tuple[str, Any]:
    if len(mapping) != 1:
        msg = f"pull {source} requires exactly one source key"
        raise UnknownSourceError(msg)
    return next(iter(mapping.items()))


def _first_series(frame: pd.DataFrame) -> pd.Series:
    if len(frame.columns) != 1:
        msg = "pull selectors must resolve to exactly one column"
        raise SelectionError(msg)
    return frame.iloc[:, 0]


def to_df_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    x: Any = None,
    raw: Any = None,
    obsm: Mapping[str, Any] | None = None,
    obsp: Mapping[str, Any] | None = None,
    layer: str | None = None,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    planner = RequestPlanner(max_matrix_values)
    descriptors: list[tuple[int, Any, str, str]] = []
    if obs is not None:
        obs_piece = evaluate_select(obs_frame(adata), obs).reset_index(drop=True)
    else:
        obs_piece = None
    if x is not None:
        token, _ = _add_matrix_request(
            planner, adata, "x", x, mode="selection", layer=layer, axis="obs", context="to_df x"
        )
        descriptors.append((token, x, "", "obs"))
    if raw is not None:
        token, _ = _add_matrix_request(planner, adata, "raw", raw, mode="selection", axis="obs", context="to_df raw")
        descriptors.append((token, raw, "raw_", "obs"))
    for key, selector in (obsm or {}).items():
        token, _ = _add_matrix_request(
            planner,
            adata,
            "obsm",
            selector,
            mode="selection",
            key=key,
            axis="obs",
            context=f"to_df obsm {key!r}",
        )
        descriptors.append((token, selector, f"{key}_", "obs"))
    for key, selector in (obsp or {}).items():
        token, _ = _add_matrix_request(
            planner,
            adata,
            "obsp",
            selector,
            mode="selection",
            key=key,
            axis="obs",
            context=f"to_df obsp {key!r}",
        )
        descriptors.append((token, selector, f"{key}_", "obs"))
    projected = planner.execute()
    pieces: list[pd.DataFrame] = [] if obs_piece is None else [obs_piece]
    for token, selector, prefix, axis in descriptors:
        frame = _decorate_projected_frame(projected[token], adata, axis=axis, request=selector)
        selected = evaluate_select(frame, selector).reset_index(drop=True)
        pieces.append(selected.add_prefix(prefix) if prefix else selected)
    if not pieces:
        return pd.DataFrame(index=adata.obs_names)
    out = pd.concat(pieces, axis=1)
    _ensure_unique([str(column) for column in out.columns], source="to_df")
    out.index = adata.obs_names
    return out


def to_tidy_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    x: Any = None,
    raw: Any = None,
    layer: str | None = None,
    obs_name: str = "obs_name",
    feature: str = "feature",
    value: str = "value",
    allow_all_features: bool = False,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    if x is not None and raw is not None:
        msg = "to_tidy accepts x or raw, not both"
        raise IncompatibleAxisError(msg)
    if x is None and raw is None and not allow_all_features:
        msg = "to_tidy requires explicit x feature selection; pass allow_all_features=True to export all features"
        raise SelectionError(msg)
    source = "raw" if raw is not None else "x"
    selector = raw if raw is not None else x
    planner = RequestPlanner(max_matrix_values)
    token, _ = _add_matrix_request(
        planner,
        adata,
        source,
        selector,
        mode="selection",
        layer=layer if source == "x" else None,
        axis="obs",
        context="to_tidy",
    )
    matrix = _decorate_projected_frame(planner.execute()[token], adata, axis="obs", request=selector)
    wide = evaluate_select(matrix, selector)
    meta = (
        evaluate_select(obs_frame(adata), obs).reset_index(drop=True)
        if obs is not None
        else pd.DataFrame(index=wide.index)
    )
    _validate_long_export_names(
        meta,
        wide,
        reserved={obs_name, feature, value},
        context="to_tidy",
    )
    combined = pd.concat([meta.reset_index(drop=True), wide.reset_index(drop=True)], axis=1)
    combined.insert(0, obs_name, adata.obs_names.to_numpy())
    tidy = _pivot_obs_major(
        combined,
        id_vars=[obs_name, *meta.columns],
        value_vars=list(wide.columns),
        var_name=feature,
        value_name=value,
    )
    return tidy[[obs_name, feature, value, *([col for col in tidy.columns if col not in {obs_name, feature, value}])]]


def _pivot_obs_major(
    data: pd.DataFrame,
    *,
    id_vars: Sequence[str],
    value_vars: Sequence[str],
    var_name: str,
    value_name: str,
) -> pd.DataFrame:
    """Pivot selected columns with rows varying before columns."""
    row_count = len(data)
    value_count = len(value_vars)
    positions = np.repeat(np.arange(row_count, dtype=np.intp), value_count)
    result = data.loc[:, list(id_vars)].iloc[positions].reset_index(drop=True)
    result[var_name] = np.tile(np.asarray(value_vars, dtype=object), row_count)
    values = data.loc[:, list(value_vars)]
    if value_count and all(isinstance(dtype, pd.SparseDtype) for dtype in values.dtypes):
        flattened = cast(Any, values).sparse.to_coo().reshape((row_count * value_count, 1), order="C")
        result[value_name] = cast(Any, pd.arrays.SparseArray).from_spmatrix(flattened)
    else:
        result[value_name] = values.to_numpy().reshape(-1)
    return result


def _check_reserved_names(columns: pd.Index | Sequence[str], reserved: set[str], *, context: str) -> None:
    collisions = [str(column) for column in columns if str(column) in reserved]
    if collisions:
        msg = f"{context} column(s) collide with reserved output name(s): {', '.join(collisions)}"
        raise NameRepairError(msg)


def _validate_long_export_names(
    meta: pd.DataFrame,
    values: pd.DataFrame,
    *,
    reserved: set[str],
    context: str,
) -> None:
    _check_reserved_names(values.columns, reserved, context=f"{context} feature")
    _check_reserved_names(meta.columns, reserved, context=f"{context} obs metadata")
    duplicated = sorted({str(column) for column in meta.columns} & {str(column) for column in values.columns})
    if duplicated:
        msg = f"Duplicate {context} column name(s) across sources: {', '.join(duplicated)}"
        raise NameRepairError(msg)


def pivot_longer_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    x: Any = None,
    raw: Any = None,
    layer: str | None = None,
    obs_name: str = "obs_name",
    names_to: str = "name",
    values_to: str = "value",
    allow_all_features: bool = False,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    if x is not None and raw is not None:
        msg = "pivot_longer accepts x or raw, not both"
        raise IncompatibleAxisError(msg)
    if x is None and raw is None and not allow_all_features:
        msg = "pivot_longer requires explicit x feature selection; pass allow_all_features=True to export all features"
        raise SelectionError(msg)
    source = "raw" if raw is not None else "x"
    selector = raw if raw is not None else x
    planner = RequestPlanner(max_matrix_values)
    token, _ = _add_matrix_request(
        planner,
        adata,
        source,
        selector,
        mode="selection",
        layer=layer if source == "x" else None,
        axis="obs",
        context="pivot_longer",
    )
    matrix = _decorate_projected_frame(planner.execute()[token], adata, axis="obs", request=selector)
    values = evaluate_select(matrix, selector)
    meta = evaluate_select(obs_frame(adata), obs) if obs is not None else pd.DataFrame(index=adata.obs_names)
    reserved = {obs_name, names_to, values_to}
    _validate_long_export_names(meta, values, reserved=reserved, context="pivot_longer")
    wide = pd.concat([meta.reset_index(drop=True), values.reset_index(drop=True)], axis=1)
    wide.insert(0, obs_name, adata.obs_names.to_numpy())
    return _pivot_obs_major(
        wide,
        id_vars=[obs_name, *meta.columns],
        value_vars=list(values.columns),
        var_name=names_to,
        value_name=values_to,
    )


def as_frame_adata(
    adata: AnnData,
    source: str,
    *,
    key: str | None = None,
    select: Any = None,
    layer: str | None = None,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    if source not in {"x", "raw", "obsm", "varm", "obsp", "varp"}:
        RequestPlanner(max_matrix_values).validate()
        return evaluate_select(source_frame(adata, source, key=key, layer=layer), select)
    axis = "var" if source in {"varm", "varp"} else "obs"
    context = f"as_frame {source}" if key is None else f"as_frame {source} {key!r}"
    planner = RequestPlanner(max_matrix_values)
    token, _ = _add_matrix_request(
        planner,
        adata,
        source,
        select,
        mode="selection",
        key=key,
        layer=layer,
        axis=axis,
        context=context,
    )
    frame = _decorate_projected_frame(planner.execute()[token], adata, axis=axis, request=select)
    selected = evaluate_select(frame, select)
    selected.index = adata.obs_names if axis == "obs" else adata.var_names
    return selected


def _check_matrix_materialization(frame: pd.DataFrame, max_matrix_values: int | None, *, context: str) -> None:
    if max_matrix_values is None:
        return
    if max_matrix_values < 0:
        msg = "max_matrix_values must be non-negative or None"
        raise AnnplyrError(msg)
    values = frame.shape[0] * frame.shape[1]
    if values > max_matrix_values:
        msg = f"{context} would materialize {values} matrix values, which exceeds max_matrix_values={max_matrix_values}"
        raise AnnplyrError(msg)


def pivot_wider(
    data: pd.DataFrame,
    *,
    id_cols: str | Sequence[str],
    names_from: str,
    values_from: str,
) -> pd.DataFrame:
    ids = [id_cols] if isinstance(id_cols, str) else list(id_cols)
    required = [*ids, names_from, values_from]
    missing = [column for column in required if column not in data.columns]
    if missing:
        msg = f"Unknown pivot_wider column(s): {', '.join(missing)}"
        raise UnknownColumnError(msg)
    if data.duplicated([*ids, names_from]).any():
        msg = "pivot_wider keys do not uniquely identify values"
        raise DuplicateNameError(msg)
    wide = data.pivot(index=ids, columns=names_from, values=values_from).reset_index()
    wide.columns.name = None
    return wide


def nest_by_adata(
    adata: AnnData,
    *,
    by: Any,
    obs: Any = None,
    var: Any = None,
    axis: str = "obs",
    name: str = "data",
) -> pd.DataFrame:
    axis = _axis(axis)
    if axis == "obs":
        frame = obs_frame(adata)
        data_selector = obs
    else:
        frame = var_frame(adata)
        data_selector = var
    keys = evaluate_select(frame, by)
    values = evaluate_select(frame, data_selector) if data_selector is not None else drop_axis_virtuals(frame)
    if name in keys.columns:
        msg = f"nest_by output column {name!r} collides with a grouping column"
        raise DuplicateNameError(msg)
    work = pd.concat([keys, values], axis=1)
    rows: list[dict[str, Any]] = []
    key_columns = list(keys.columns)
    value_columns = list(values.columns)
    for _, key_row in keys.drop_duplicates().iterrows():
        mask = pd.Series(True, index=keys.index)
        for column, value in key_row.items():
            mask &= keys[column].isna() if pd.isna(value) else keys[column].eq(value)
        row = {str(column): key_row[column] for column in key_columns}
        row[name] = work.loc[mask, value_columns].reset_index(drop=True)
        rows.append(row)
    return pd.DataFrame(rows, columns=[*key_columns, name])


def drop_axis_virtuals(frame: pd.DataFrame) -> pd.DataFrame:
    return evaluate_select(frame, None)


def unnest(data: pd.DataFrame, column: str) -> pd.DataFrame:
    if column not in data.columns:
        msg = f"Unknown unnest column: {column!r}"
        raise UnknownColumnError(msg)
    rows: list[pd.DataFrame] = []
    outer_columns = [name for name in data.columns if name != column]
    inner_columns: list[str] | None = None
    for _, row in data.iterrows():
        nested = row[column]
        if not isinstance(nested, pd.DataFrame):
            msg = f"unnest column {column!r} must contain pandas DataFrame objects"
            raise SelectionError(msg)
        if inner_columns is None:
            inner_columns = [str(c) for c in nested.columns]
        if nested.empty:
            continue
        prefix = pd.DataFrame({name: [row[name]] * len(nested) for name in outer_columns})
        rows.append(pd.concat([prefix.reset_index(drop=True), nested.reset_index(drop=True)], axis=1))
    if not rows:
        extra = inner_columns if inner_columns is not None else []
        return pd.DataFrame(columns=[*outer_columns, *extra])
    return pd.concat(rows, ignore_index=True)


def nest(
    data: pd.DataFrame,
    *,
    by: str | Sequence[str],
    columns: Sequence[str] | None = None,
    name: str = "data",
) -> pd.DataFrame:
    by_columns = [by] if isinstance(by, str) else list(by)
    _check_dataframe_columns(data, by_columns, context="nest")
    value_columns = (
        [column for column in data.columns if column not in by_columns] if columns is None else list(columns)
    )
    _check_dataframe_columns(data, value_columns, context="nest")
    if name in by_columns:
        msg = f"nest output column {name!r} collides with a grouping column"
        raise DuplicateNameError(msg)
    rows: list[dict[str, Any]] = []
    for key, group in data.groupby(by_columns, sort=False, dropna=False):
        key_values = _group_key_values(key, len(by_columns))
        row = dict(zip(by_columns, key_values, strict=True))
        row[name] = group.loc[:, value_columns].reset_index(drop=True)
        rows.append(row)
    return pd.DataFrame(rows, columns=[*by_columns, name])


def chop(
    data: pd.DataFrame,
    columns: str | Sequence[str],
    *,
    by: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    selected = [columns] if isinstance(columns, str) else list(columns)
    _check_dataframe_columns(data, selected, context="chop")
    by_columns = (
        [column for column in data.columns if column not in selected]
        if by is None
        else ([by] if isinstance(by, str) else list(by))
    )
    _check_dataframe_columns(data, by_columns, context="chop")
    rows: list[dict[str, Any]] = []
    for key, group in data.groupby(by_columns, sort=False, dropna=False):
        key_values = _group_key_values(key, len(by_columns))
        row = dict(zip(by_columns, key_values, strict=True))
        for column in selected:
            row[column] = group[column].tolist()
        rows.append(row)
    return pd.DataFrame(rows, columns=[*by_columns, *selected])


def unchop(data: pd.DataFrame, columns: str | Sequence[str], *, keep_empty: bool = False) -> pd.DataFrame:
    selected = [columns] if isinstance(columns, str) else list(columns)
    _check_dataframe_columns(data, selected, context="unchop")
    rows: list[dict[str, Any]] = []
    for _, row in data.iterrows():
        size = max((_list_like_len(row[column]) for column in selected), default=1)
        if size == 0 and not keep_empty:
            continue
        size = max(size, 1)
        for i in range(size):
            out = cast(dict[str, Any], row.drop(labels=selected).to_dict())
            for column in selected:
                values = _as_list_like(row[column])
                out[column] = values[i] if i < len(values) else pd.NA
            rows.append(out)
    return pd.DataFrame(rows, columns=list(data.columns))


def unnest_longer(
    data: pd.DataFrame,
    column: str,
    *,
    values_to: str | None = None,
    indices_to: str | None = None,
    keep_empty: bool = False,
) -> pd.DataFrame:
    _check_dataframe_columns(data, [column], context="unnest_longer")
    value_name = values_to or column
    rows: list[dict[str, Any]] = []
    for _, row in data.iterrows():
        values = _as_list_like(row[column])
        if not values and keep_empty:
            values = [pd.NA]
        for index, value in enumerate(values):
            out = cast(dict[str, Any], row.drop(labels=[column]).to_dict())
            out[value_name] = value
            if indices_to is not None:
                out[indices_to] = index
            rows.append(out)
    base_columns = [name for name in data.columns if name != column]
    return pd.DataFrame(rows, columns=[*base_columns, value_name, *([] if indices_to is None else [indices_to])])


def unnest_wider(data: pd.DataFrame, column: str, *, names_sep: str | None = None) -> pd.DataFrame:
    _check_dataframe_columns(data, [column], context="unnest_wider")
    pieces: list[pd.DataFrame] = []
    for value in data[column]:
        if isinstance(value, pd.DataFrame):
            piece = value.reset_index(drop=True).iloc[:1]
        elif isinstance(value, pd.Series):
            piece = value.to_frame().T.reset_index(drop=True)
        elif isinstance(value, Mapping):
            piece = pd.DataFrame([value])
        else:
            piece = pd.DataFrame([dict(enumerate(_as_list_like(value)))])
        pieces.append(piece.reset_index(drop=True))
    wider = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(index=data.index)
    wider.columns = [f"{column}{names_sep}{name}" if names_sep is not None else str(name) for name in wider.columns]
    out = data.drop(columns=[column]).reset_index(drop=True)
    _check_output_columns_available(out.columns, wider.columns, context="unnest_wider")
    return pd.concat([out, wider], axis=1)


def pack(data: pd.DataFrame, column: str, columns: str | Sequence[str]) -> pd.DataFrame:
    selected = [columns] if isinstance(columns, str) else list(columns)
    _check_dataframe_columns(data, selected, context="pack")
    out = data.drop(columns=selected).copy()
    _check_output_columns_available(out.columns, [column], context="pack")
    packed = data.loc[:, selected].apply(lambda row: row.to_dict(), axis=1)
    insert_at = _first_column_position(data, selected, fallback=len(out.columns))
    out.insert(min(insert_at, len(out.columns)), column, packed)
    return out


def unpack(
    data: pd.DataFrame,
    column: str,
    *,
    names_sep: str | None = None,
    remove: bool = True,
) -> pd.DataFrame:
    _check_dataframe_columns(data, [column], context="unpack")
    wider = unnest_wider(data.loc[:, [column]], column, names_sep=names_sep)
    out = data.drop(columns=[column]).copy() if remove else data.copy()
    _check_output_columns_available(out.columns, wider.columns, context="unpack")
    insert_at = _first_column_position(data, [column], fallback=len(out.columns))
    for offset, wider_column in enumerate(wider.columns):
        out.insert(min(insert_at + offset, len(out.columns)), str(wider_column), wider[wider_column].to_numpy())
    return out


def hoist(data: pd.DataFrame, column: str, **paths: str | int) -> pd.DataFrame:
    _check_dataframe_columns(data, [column], context="hoist")
    out = data.copy()
    for name, path in paths.items():
        out[name] = out[column].map(lambda value, path=path: _pluck(value, path))
    return out


def drop_na(data: pd.DataFrame, columns: str | Sequence[str] | None = None) -> pd.DataFrame:
    subset = None if columns is None else ([columns] if isinstance(columns, str) else list(columns))
    _check_dataframe_columns(data, subset or list(data.columns), context="drop_na")
    return data.dropna(subset=subset).reset_index(drop=True)


def fill(
    data: pd.DataFrame,
    columns: str | Sequence[str],
    *,
    direction: str = "down",
) -> pd.DataFrame:
    selected = [columns] if isinstance(columns, str) else list(columns)
    _check_dataframe_columns(data, selected, context="fill")
    out = data.copy()
    if direction == "down":
        out[selected] = out[selected].ffill()
    elif direction == "up":
        out[selected] = out[selected].bfill()
    elif direction == "downup":
        out[selected] = out[selected].ffill().bfill()
    elif direction == "updown":
        out[selected] = out[selected].bfill().ffill()
    else:
        msg = "direction must be 'down', 'up', 'downup', or 'updown'"
        raise SelectionError(msg)
    return out


def separate(
    data: pd.DataFrame,
    column: str,
    *,
    into: Sequence[str],
    sep: str = "_",
    remove: bool = True,
) -> pd.DataFrame:
    """Split a single column into multiple columns using a regular expression separator.

    Parameters
    ----------
    data:
        Input DataFrame.
    column:
        Name of the column to split.
    into:
        Names for the output columns produced by the split.
    sep:
        Regular expression used to split the column value. Defaults to ``"_"``.
        Both :func:`separate` and :func:`separate_rows` treat *sep* as a
        regular expression, matching tidyr's documented behaviour.
    remove:
        If ``True`` (the default), remove the source column from the result.
    """
    _check_dataframe_columns(data, [column], context="separate")
    out = data.copy()
    max_splits = len(into) - 1
    split_series = out[column].map(
        lambda value: re.split(sep, str(value), maxsplit=max_splits) if not pd.isna(value) else [pd.NA] * len(into)
    )
    split = pd.DataFrame(split_series.tolist(), index=out.index)
    for index, name in enumerate(into):
        out[str(name)] = split[index] if index in split.columns else pd.NA
    if remove:
        out = out.drop(columns=[column])
    return out


def separate_rows(data: pd.DataFrame, columns: str | Sequence[str], *, sep: str = ",") -> pd.DataFrame:
    """Separate delimited values in one or more columns into multiple rows.

    Parameters
    ----------
    data:
        Input DataFrame.
    columns:
        Column name or sequence of column names whose values should be expanded.
    sep:
        Regular expression used to split each value. Defaults to ``","``
        (a literal comma). Both :func:`separate` and :func:`separate_rows` treat
        *sep* as a regular expression, matching tidyr's documented behaviour.
    """
    selected = [columns] if isinstance(columns, str) else list(columns)
    _check_dataframe_columns(data, selected, context="separate_rows")
    out = data.copy()
    for column in selected:
        out[column] = out[column].map(lambda value: [] if pd.isna(value) else re.split(sep, str(value)))
    return out.explode(selected, ignore_index=True)


def extract(
    data: pd.DataFrame,
    column: str,
    *,
    into: Sequence[str],
    regex: str,
    remove: bool = True,
) -> pd.DataFrame:
    _check_dataframe_columns(data, [column], context="extract")
    extracted = data[column].astype("string").str.extract(regex, expand=True)
    out = data.copy()
    for index, name in enumerate(into):
        out[str(name)] = extracted[index] if index in extracted.columns else pd.NA
    if remove:
        out = out.drop(columns=[column])
    return out


def unite(
    data: pd.DataFrame,
    column: str,
    columns: Sequence[str],
    *,
    sep: str = "_",
    remove: bool = True,
    na_rm: bool = False,
) -> pd.DataFrame:
    selected = list(columns)
    _check_dataframe_columns(data, selected, context="unite")
    out = data.copy()

    def _join(row: pd.Series) -> str:
        values = [row[name] for name in selected]
        if na_rm:
            values = [value for value in values if not pd.isna(value)]
        return sep.join("" if pd.isna(value) else str(value) for value in values)

    out[column] = out.apply(_join, axis=1)
    if remove:
        out = out.drop(columns=selected)
    insert_at = _first_column_position(data, selected, fallback=len(out.columns))
    series = out.pop(column)
    out.insert(min(insert_at, len(out.columns)), column, series)
    return out


def _first_column_position(data: pd.DataFrame, columns: Sequence[str], *, fallback: int) -> int:
    positions: list[int] = []
    for name in columns:
        loc = data.columns.get_loc(name)
        if not isinstance(loc, int | np.integer):
            msg = f"Column {name!r} is duplicated; position-sensitive operation requires unique columns"
            raise DuplicateNameError(msg)
        positions.append(int(loc))
    return min(positions, default=fallback)


def _as_list_like(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, pd.Series):
        return value.tolist()
    if pd.isna(value):
        return []
    return [value]


def _list_like_len(value: Any) -> int:
    return len(_as_list_like(value))


def _pluck(value: Any, path: str | int) -> Any:
    try:
        if isinstance(value, Mapping):
            return value.get(path, pd.NA)
        if isinstance(value, pd.Series):
            return value.get(path, pd.NA)
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            return value[int(path)]
    except (IndexError, KeyError, TypeError, ValueError):
        return pd.NA
    return pd.NA


def _group_key_values(key: Any, n_columns: int) -> tuple[Any, ...]:
    if n_columns == 1:
        return (key[0],) if isinstance(key, tuple) else (key,)
    return tuple(key)


def _check_dataframe_columns(data: pd.DataFrame, columns: Sequence[str], *, context: str) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        msg = f"Unknown {context} column(s): {', '.join(missing)}"
        raise UnknownColumnError(msg)


def _check_output_columns_available(
    existing: pd.Index | Sequence[str], new: pd.Index | Sequence[str], *, context: str
) -> None:
    collisions = sorted({str(column) for column in existing} & {str(column) for column in new})
    if collisions:
        msg = f"{context} output column(s) would duplicate existing column(s): {', '.join(collisions)}"
        raise DuplicateNameError(msg)
