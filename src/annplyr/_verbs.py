from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any, cast

import narwhals as nw
import numpy as np
import pandas as pd
from anndata import AnnData

from annplyr._errors import (
    DuplicateNameError,
    IncompatibleAxisError,
    SelectionError,
    UnknownColumnError,
    UnknownSourceError,
)
from annplyr._expr import Desc
from annplyr._frames import (
    evaluate_assignments,
    evaluate_filter,
    evaluate_select,
    intersect_ordered,
    obs_frame,
    obsm_frame,
    var_frame,
    varm_frame,
    with_row_number,
    x_frame,
)


def _subset(adata: AnnData, obs_idx: Any, var_idx: Any, *, copy: bool = False) -> AnnData:
    out = adata[obs_idx, var_idx]
    return out.copy() if copy else out


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
    obs_names: Any = None,
    var_names: Any = None,
    obsm: Mapping[str, Any] | None = None,
    varm: Mapping[str, Any] | None = None,
    layer: str | None = None,
    copy: bool = False,
) -> AnnData:
    obs_indices: list[pd.Index] = []
    var_indices: list[pd.Index] = []

    if obs is not None:
        obs_indices.append(evaluate_filter(obs_frame(adata), obs))
    if obs_names is not None:
        obs_indices.append(
            evaluate_filter(pd.DataFrame({"obs_names": adata.obs_names}, index=adata.obs_names), obs_names)
        )
    if x is not None:
        obs_indices.append(evaluate_filter(x_frame(adata, layer=layer), x))
    for key, predicates in (obsm or {}).items():
        obs_indices.append(evaluate_filter(obsm_frame(adata, key), predicates))

    if var is not None:
        var_indices.append(evaluate_filter(var_frame(adata), var))
    if var_names is not None:
        var_indices.append(
            evaluate_filter(pd.DataFrame({"var_names": adata.var_names}, index=adata.var_names), var_names)
        )
    for key, predicates in (varm or {}).items():
        var_indices.append(evaluate_filter(varm_frame(adata, key), predicates))

    obs_idx = intersect_ordered(adata.obs_names, *obs_indices)
    var_idx = intersect_ordered(adata.var_names, *var_indices)
    return _subset(adata, obs_idx, var_idx, copy=copy)


def select_adata(adata: AnnData, obs: Any = None, var: Any = None, x: Any = None, copy: bool = False) -> AnnData:
    obs_columns = (
        _real_columns(evaluate_select(obs_frame(adata), obs).columns, _obs_table(adata).columns)
        if obs is not None
        else _obs_table(adata).columns
    )
    var_columns = (
        _real_columns(evaluate_select(var_frame(adata), var).columns, _var_table(adata).columns)
        if var is not None
        else _var_table(adata).columns
    )
    var_names = (
        _real_columns(evaluate_select(x_frame(adata), x).columns, adata.var_names) if x is not None else adata.var_names
    )

    out = adata[:, var_names].copy() if copy else adata[:, var_names]
    out.obs = _obs_table(out).loc[:, list(obs_columns)].copy()
    out.var = _var_table(out).loc[:, list(var_columns)].copy()
    return out


def _real_columns(selected: pd.Index, available: pd.Index) -> list[str]:
    return [column for column in selected if column in available]


def rename_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, str] | None = None,
    var: Mapping[str, str] | None = None,
    x: Mapping[str, str] | None = None,
    copy: bool = True,
) -> AnnData:
    out = adata.copy() if copy else adata
    if obs:
        out.obs = _obs_table(out).rename(columns=_rename_mapping(_obs_table(out).columns, obs, source="obs")).copy()
    if var:
        out.var = _var_table(out).rename(columns=_rename_mapping(_var_table(out).columns, var, source="var")).copy()
    if x:
        out.var_names = _renamed_names(out.var_names, x, source="x")
    return out


def rename_with_adata(
    adata: AnnData,
    func: Callable[[str], str],
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    copy: bool = True,
) -> AnnData:
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
    return rename_adata(adata, obs=obs_mapping, var=var_mapping, x=x_mapping, copy=copy)


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
    copy: bool = True,
) -> AnnData:
    out = adata.copy() if copy else adata
    if obs is not None:
        out.obs = (
            _obs_table(out)
            .loc[
                :,
                _relocated_order(
                    _obs_table(out).columns,
                    _selected_columns(obs_frame(out), obs, _obs_table(out).columns),
                    before,
                    after,
                ),
            ]
            .copy()
        )
    if var is not None:
        out.var = (
            _var_table(out)
            .loc[
                :,
                _relocated_order(
                    _var_table(out).columns,
                    _selected_columns(var_frame(out), var, _var_table(out).columns),
                    before,
                    after,
                ),
            ]
            .copy()
        )
    if x is not None:
        out = out[
            :, _relocated_order(out.var_names, _selected_columns(x_frame(out), x, out.var_names), before, after)
        ].copy()
    return out


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
        msg = f"Unknown relocate anchor: {before!r}"
        raise UnknownColumnError(msg)
    elif after is not None:
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
) -> AnnData:
    axis = _axis(axis)
    if axis == "obs":
        frame, selector, available = _distinct_source(adata, obs=obs, x=x, axis=axis)
        selected = evaluate_select(frame, selector)
        obs_idx = selected.drop_duplicates(keep="first").index
        out = _subset(adata, obs_idx, slice(None), copy=copy)
        if not keep_all and obs is not None:
            out.obs = _obs_table(out).loc[:, _real_columns(selected.columns, available)].copy()
        return out

    frame, selector, available = _distinct_source(adata, var=var, axis=axis)
    selected = evaluate_select(frame, selector)
    var_idx = selected.drop_duplicates(keep="first").index
    out = _subset(adata, slice(None), var_idx, copy=copy)
    if not keep_all and var is not None:
        out.var = _var_table(out).loc[:, _real_columns(selected.columns, available)].copy()
    return out


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

    work = with_row_number(frame)
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
    return expr


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
    obsm: Mapping[str, Any] | None = None,
    varm: Mapping[str, Any] | None = None,
    layer: str | None = None,
    copy: bool = False,
) -> AnnData:
    obs_idx = adata.obs_names
    var_idx = adata.var_names

    for frame, by in _obs_sort_frames(adata, obs=obs, x=x, obsm=obsm, layer=layer):
        obs_idx = _sort_values_for_frame(frame.loc[obs_idx], by)
    for frame, by in _var_sort_frames(adata, var=var, varm=varm):
        var_idx = _sort_values_for_frame(frame.loc[var_idx], by)

    return _subset(adata, obs_idx, var_idx, copy=copy)


def _obs_sort_frames(
    adata: AnnData,
    *,
    obs: Any = None,
    x: Any = None,
    obsm: Mapping[str, Any] | None = None,
    layer: str | None = None,
):
    if obs is not None:
        yield obs_frame(adata), obs
    if x is not None:
        yield x_frame(adata, layer=layer), x
    for key, by in (obsm or {}).items():
        yield obsm_frame(adata, key), by


def _var_sort_frames(adata: AnnData, *, var: Any = None, varm: Mapping[str, Any] | None = None):
    if var is not None:
        yield var_frame(adata), var
    for key, by in (varm or {}).items():
        yield varm_frame(adata, key), by


def slice_adata(adata: AnnData, *indices: Any, axis: str = "obs", copy: bool = False) -> AnnData:
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


def slice_head_adata(adata: AnnData, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
    return slice_adata(adata, slice(0, n), axis=axis, copy=copy)


def slice_tail_adata(adata: AnnData, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
    return slice_adata(adata, slice(-n, None), axis=axis, copy=copy)


def slice_min_adata(adata: AnnData, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
    axis = _axis(axis)
    arranged = arrange_adata(adata, obs=by) if axis == "obs" else arrange_adata(adata, var=by)
    return slice_head_adata(arranged, n=n, axis=axis, copy=copy)


def slice_max_adata(adata: AnnData, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
    axis = _axis(axis)
    arranged = arrange_adata(adata, obs=Desc(by)) if axis == "obs" else arrange_adata(adata, var=Desc(by))
    return slice_head_adata(arranged, n=n, axis=axis, copy=copy)


def slice_sample_adata(
    adata: AnnData,
    n: int | None = None,
    *,
    prop: float | None = None,
    replace: bool = False,
    random_state: int | None = None,
    axis: str = "obs",
    copy: bool = False,
) -> AnnData:
    axis = _axis(axis)
    size = adata.n_obs if axis == "obs" else adata.n_vars
    if n is None:
        n = int(round(size * prop)) if prop is not None else min(size, 1)
    rng = np.random.default_rng(random_state)
    selected = rng.choice(size, size=n, replace=replace)
    return slice_adata(adata, selected.tolist(), axis=axis, copy=copy)


def mutate_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    var: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
    layer: str | None = None,
    inplace: bool = False,
) -> AnnData:
    out = adata if inplace else adata.copy()
    for frame, assignments in _obs_assignment_frames(out, obs=obs, x=x, obsm=obsm, layer=layer):
        values = evaluate_assignments(frame, assignments)
        for column in values.columns:
            _obs_table(out)[column] = values[column].to_numpy()
    for frame, assignments in _var_assignment_frames(out, var=var, varm=varm):
        values = evaluate_assignments(frame, assignments)
        for column in values.columns:
            _var_table(out)[column] = values[column].to_numpy()
    return out


def transmute_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    var: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
    layer: str | None = None,
) -> AnnData:
    out = mutate_adata(adata, obs=obs, var=var, x=x, obsm=obsm, varm=varm, layer=layer, inplace=False)
    obs_columns = _assignment_names(obs, x, obsm)
    var_columns = _assignment_names(var, None, varm)
    if obs_columns:
        out.obs = _obs_table(out).loc[:, obs_columns].copy()
    if var_columns:
        out.var = _var_table(out).loc[:, var_columns].copy()
    return out


def _assignment_names(
    direct: Mapping[str, Any] | None,
    matrix: Mapping[str, Any] | None,
    keyed: Mapping[str, Mapping[str, Any]] | None,
) -> list[str]:
    names: list[str] = []
    names.extend((direct or {}).keys())
    names.extend((matrix or {}).keys())
    for assignments in (keyed or {}).values():
        names.extend(assignments.keys())
    return names


def _obs_assignment_frames(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    layer: str | None = None,
):
    if obs:
        yield obs_frame(adata), obs
    if x:
        yield x_frame(adata, layer=layer), x
    for key, assignments in (obsm or {}).items():
        yield obsm_frame(adata, key), assignments


def _var_assignment_frames(
    adata: AnnData,
    *,
    var: Mapping[str, Any] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
):
    if var:
        yield var_frame(adata), var
    for key, assignments in (varm or {}).items():
        yield varm_frame(adata, key), assignments


def summarize_adata(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    var: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
    by: Any = None,
    layer: str | None = None,
) -> pd.DataFrame:
    obs_axis_requested = any(source is not None for source in (obs, x, obsm))
    var_axis_requested = any(source is not None for source in (var, varm))
    if obs_axis_requested and var_axis_requested:
        msg = "summarize accepts obs-axis or var-axis sources, not both at once"
        raise IncompatibleAxisError(msg)

    if var_axis_requested:
        sources = _var_summary_sources(adata, var=var, varm=varm)
        return summarize_sources(var_frame(adata), sources, by=by)

    sources = _obs_summary_sources(adata, obs=obs, x=x, obsm=obsm, layer=layer)
    return summarize_sources(obs_frame(adata), sources, by=by)


def _obs_summary_sources(
    adata: AnnData,
    *,
    obs: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    obsm: Mapping[str, Mapping[str, Any]] | None = None,
    layer: str | None = None,
) -> list[tuple[pd.DataFrame, Mapping[str, Any]]]:
    sources: list[tuple[pd.DataFrame, Mapping[str, Any]]] = []
    if obs:
        sources.append((obs_frame(adata), obs))
    if x:
        sources.append((x_frame(adata, layer=layer), x))
    for key, assignments in (obsm or {}).items():
        sources.append((obsm_frame(adata, key), assignments))
    return sources


def _var_summary_sources(
    adata: AnnData,
    *,
    var: Mapping[str, Any] | None = None,
    varm: Mapping[str, Mapping[str, Any]] | None = None,
) -> list[tuple[pd.DataFrame, Mapping[str, Any]]]:
    sources: list[tuple[pd.DataFrame, Mapping[str, Any]]] = []
    if var:
        sources.append((var_frame(adata), var))
    for key, assignments in (varm or {}).items():
        sources.append((varm_frame(adata, key), assignments))
    return sources


def summarize_sources(
    by_source: pd.DataFrame,
    sources: list[tuple[pd.DataFrame, Mapping[str, Any]]],
    *,
    by: Any = None,
) -> pd.DataFrame:
    _, by_columns = prepare_by_frame(by_source, by)
    by_values = evaluate_select(by_source, by) if by_columns else pd.DataFrame(index=by_source.index)
    pieces: list[pd.DataFrame] = []

    for frame, assignments in sources:
        work = frame.copy()
        for column in by_columns:
            work[column] = by_values[column].to_numpy()
        pieces.append(summarize_frame(work, assignments=assignments, by=by_columns))

    if not pieces:
        return count_frame(by_source, by=by)
    return _merge_summary_pieces(pieces, by_columns)


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
    assignments = assignments or {}
    work, by_columns = prepare_by_frame(frame, by)
    exprs = [
        expr.alias(name) if hasattr(expr, "alias") else nw.col(expr).alias(name) for name, expr in assignments.items()
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
        work[column] = by_frame[column].to_numpy()
    return work, by_columns


def count_adata(adata: AnnData, by: Any = None, *, axis: str = "obs", name: str = "n") -> pd.DataFrame:
    frame = obs_frame(adata) if _axis(axis) == "obs" else var_frame(adata)
    return count_frame(frame, by=by, name=name)


def tally_adata(adata: AnnData, by: Any = None, *, axis: str = "obs", name: str = "n") -> pd.DataFrame:
    return count_adata(adata, by=by, axis=axis, name=name)


def add_count_adata(
    adata: AnnData,
    by: Any = None,
    *,
    axis: str = "obs",
    name: str = "n",
    inplace: bool = False,
) -> AnnData:
    axis = _axis(axis)
    out = adata if inplace else adata.copy()
    frame = obs_frame(out) if axis == "obs" else var_frame(out)
    values = _count_values(frame, by=by)
    if axis == "obs":
        _obs_table(out)[name] = values.to_numpy()
    else:
        _var_table(out)[name] = values.to_numpy()
    return out


def _count_values(frame: pd.DataFrame, by: Any = None) -> pd.Series:
    by_frame = evaluate_select(frame, by) if by is not None else pd.DataFrame(index=frame.index)
    by_columns = list(by_frame.columns)
    if not by_columns:
        return pd.Series(len(frame), index=frame.index)
    work = by_frame.copy()
    work["__annplyr_count_row__"] = np.arange(len(work))
    return work.groupby(by_columns, sort=False, dropna=False)["__annplyr_count_row__"].transform("size")


def count_frame(frame: pd.DataFrame, by: Any = None, *, name: str = "n") -> pd.DataFrame:
    work, by_columns = prepare_by_frame(frame, by)
    if not by_columns:
        return pd.DataFrame({name: [len(frame)]})
    return work.groupby(by_columns, sort=False, dropna=False).size().reset_index(name=name)


def pull_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    var: Any = None,
    x: Any = None,
    obsm: Mapping[str, Any] | None = None,
    varm: Mapping[str, Any] | None = None,
    layer: str | None = None,
) -> pd.Series:
    provided = [value is not None for value in [obs, var, x, obsm, varm]]
    if sum(provided) != 1:
        msg = "pull requires exactly one source"
        raise UnknownSourceError(msg)
    if obs is not None:
        return _first_series(evaluate_select(obs_frame(adata), obs))
    if var is not None:
        return _first_series(evaluate_select(var_frame(adata), var))
    if x is not None:
        return _first_series(evaluate_select(x_frame(adata, layer=layer), x))
    if obsm is not None:
        key, selector = next(iter(obsm.items()))
        return _first_series(evaluate_select(obsm_frame(adata, key), selector))
    key, selector = next(iter((varm or {}).items()))
    return _first_series(evaluate_select(varm_frame(adata, key), selector))


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
    obsm: Mapping[str, Any] | None = None,
    layer: str | None = None,
) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    if obs is not None:
        pieces.append(evaluate_select(obs_frame(adata), obs))
    if x is not None:
        pieces.append(evaluate_select(x_frame(adata, layer=layer), x))
    for key, selector in (obsm or {}).items():
        selected = evaluate_select(obsm_frame(adata, key), selector)
        selected = selected.add_prefix(f"{key}_")
        pieces.append(selected)
    if not pieces:
        return pd.DataFrame(index=adata.obs_names)
    out = pd.concat(pieces, axis=1)
    out.index = adata.obs_names
    return out


def to_tidy_adata(
    adata: AnnData,
    *,
    obs: Any = None,
    x: Any = None,
    layer: str | None = None,
    obs_name: str = "obs_name",
    feature: str = "feature",
    value: str = "value",
) -> pd.DataFrame:
    wide = evaluate_select(x_frame(adata, layer=layer), x if x is not None else list(adata.var_names))
    wide = wide.copy()
    wide[obs_name] = adata.obs_names.to_numpy()
    tidy = wide.melt(id_vars=obs_name, var_name=feature, value_name=value)
    if obs is not None:
        meta = evaluate_select(obs_frame(adata), obs).copy()
        meta[obs_name] = adata.obs_names.to_numpy()
        tidy = tidy.merge(meta, on=obs_name, how="left")
    return tidy[[obs_name, feature, value, *([col for col in tidy.columns if col not in {obs_name, feature, value}])]]
