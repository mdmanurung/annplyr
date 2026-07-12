from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import narwhals as nw
import numpy as np
import pandas as pd
from anndata import AnnData

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


def _axis(axis: str) -> str:
    if axis in {"obs", "observation", "observations", "0", 0}:  # type: ignore[comparison-overlap]
        return "obs"
    if axis in {"var", "variable", "variables", "1", 1}:  # type: ignore[comparison-overlap]
        return "var"
    msg = "axis must be 'obs' or 'var'"
    raise ValueError(msg)


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
        _real_columns(evaluate_select(obs_frame(adata), obs).columns, adata.obs.columns)
        if obs is not None
        else adata.obs.columns
    )
    var_columns = (
        _real_columns(evaluate_select(var_frame(adata), var).columns, adata.var.columns)
        if var is not None
        else adata.var.columns
    )
    var_names = (
        _real_columns(evaluate_select(x_frame(adata), x).columns, adata.var_names) if x is not None else adata.var_names
    )

    out = adata[:, var_names].copy() if copy else adata[:, var_names]
    out.obs = out.obs.loc[:, list(obs_columns)].copy()
    out.var = out.var.loc[:, list(var_columns)].copy()
    return out


def _real_columns(selected: pd.Index, available: pd.Index) -> list[str]:
    return [column for column in selected if column in available]


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
    arranged = arrange_adata(adata, **{_axis(axis): by})
    return slice_head_adata(arranged, n=n, axis=axis, copy=copy)


def slice_max_adata(adata: AnnData, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
    arranged = arrange_adata(adata, **{_axis(axis): Desc(by)})
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
            out.obs[column] = values[column].to_numpy()
    for frame, assignments in _var_assignment_frames(out, var=var, varm=varm):
        values = evaluate_assignments(frame, assignments)
        for column in values.columns:
            out.var[column] = values[column].to_numpy()
    return out


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
        raise ValueError(msg)

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
        raise ValueError(msg)
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
        raise ValueError(msg)
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
