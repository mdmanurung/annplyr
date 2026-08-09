from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any, TypeVar, cast

import numpy as np
import pandas as pd
from anndata import AnnData

from annplyr._errors import AnnplyrError, IncompatibleAxisError, SelectionError
from annplyr._expr import col
from annplyr._frames import evaluate_assignments, evaluate_select, obs_frame, var_frame
from annplyr._groups import GroupPlan, GroupSpec
from annplyr._verbs import (
    _axis_positions,
    _desc_order_by,
    _slice_selector,
    _sort_values_for_frame,
    _subset_positions,
    _validate_slice_n,
    add_count_adata,
    anti_join_adata,
    arrange_adata,
    distinct_adata,
    filter_adata,
    full_join_adata,
    inner_join_adata,
    left_join_adata,
    mutate_adata,
    relocate_adata,
    rename_adata,
    rename_with_adata,
    right_join_adata,
    select_adata,
    semi_join_adata,
    summarize_adata,
    transmute_adata,
)

T = TypeVar("T")


class GroupedAnnData:
    """Persistent positional grouping over one AnnData metadata axis."""

    def __init__(
        self,
        adata: AnnData,
        *,
        obs: Any = None,
        var: Any = None,
        spec: GroupSpec | None = None,
    ):
        self._adata = adata
        self._spec = spec if spec is not None else GroupSpec.resolve(adata, obs=obs, var=var)
        self._spec.validate(adata)

    @property
    def _axis(self) -> str:
        return self._spec.axis

    @property
    def _by(self) -> list[str]:
        return list(self._spec.columns)

    def _plan(self) -> GroupPlan:
        return GroupPlan.build(self._adata, self._spec)

    def _wrap(self, adata: AnnData, *, spec: GroupSpec | None = None) -> GroupedAnnData:
        return GroupedAnnData(adata, spec=spec or self._spec)

    def _same_shape_result(
        self,
        adata: AnnData,
        *,
        inplace: bool,
        spec: GroupSpec | None = None,
    ) -> GroupedAnnData:
        updated = spec or self._spec
        if inplace:
            if adata is not self._adata:
                raise AnnplyrError("grouped in-place verb returned a replacement AnnData object")
            self._spec = updated
            return self
        return self._wrap(adata, spec=updated)

    def __iter__(self) -> Iterator[tuple[dict[str, Any], AnnData]]:
        plan = self._plan()
        for row, positions in zip(range(len(plan.keys)), plan.positions, strict=True):
            key = {column: plan.keys.iloc[row][column] for column in self._spec.columns}
            if self._axis == "obs":
                group = _subset_positions(self._adata, positions, None, copy=False)
            else:
                group = _subset_positions(self._adata, None, positions, copy=False)
            yield key, group

    def group_vars(self) -> list[str]:
        self._spec.validate(self._adata)
        return list(self._spec.columns)

    def group_keys(self) -> pd.DataFrame:
        return GroupPlan.keys_for(self._adata, self._spec)

    def group_data(self) -> pd.DataFrame:
        return self._plan().group_data()

    def ungroup(self) -> AnnData:
        return self._adata

    def filter(
        self,
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
    ) -> GroupedAnnData:
        plan = self._plan()
        result = filter_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obs_names=obs_names,
            var_names=var_names,
            obsm=obsm,
            varm=varm,
            layer=layer,
            copy=copy,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def select(self, obs: Any = None, var: Any = None, x: Any = None, copy: bool = True) -> GroupedAnnData:
        self._spec.validate(self._adata)
        if self._axis == "obs" and obs is not None:
            obs = self._retain_group_keys(obs, axis="obs")
        elif self._axis == "var" and var is not None:
            var = self._retain_group_keys(var, axis="var")
        result = select_adata(self._adata, obs=obs, var=var, x=x, copy=copy)
        return self._wrap(result)

    def _retain_group_keys(self, selector: Any, *, axis: str) -> list[str]:
        frame = obs_frame(self._adata) if axis == "obs" else var_frame(self._adata)
        selected = [str(column) for column in evaluate_select(frame, selector).columns]
        omitted = [column for column in self._spec.columns if column not in selected]
        return [*omitted, *selected]

    def rename(
        self,
        obs: Mapping[str, str] | None = None,
        var: Mapping[str, str] | None = None,
        x: Mapping[str, str] | None = None,
        *,
        inplace: bool = False,
    ) -> GroupedAnnData:
        axis_mapping = obs if self._axis == "obs" else var
        update = {old: new for new, old in (axis_mapping or {}).items()}
        spec = self._spec.renamed(update)
        result = rename_adata(self._adata, obs=obs, var=var, x=x, inplace=inplace)
        return self._same_shape_result(result, inplace=inplace, spec=spec)

    def rename_with(
        self,
        func: Callable[[str], str],
        *,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        inplace: bool = False,
    ) -> GroupedAnnData:
        before = list(cast(pd.DataFrame, self._adata.obs if self._axis == "obs" else self._adata.var).columns)
        result = rename_with_adata(
            self._adata,
            func,
            obs=obs,
            var=var,
            x=x,
            inplace=inplace,
        )
        after = list(cast(pd.DataFrame, result.obs if self._axis == "obs" else result.var).columns)
        spec = self._spec.renamed(dict(zip(before, after, strict=True)))
        return self._same_shape_result(result, inplace=inplace, spec=spec)

    def relocate(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        *,
        before: str | None = None,
        after: str | None = None,
        inplace: bool = False,
    ) -> GroupedAnnData:
        result = relocate_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            before=before,
            after=after,
            inplace=inplace,
        )
        return self._same_shape_result(result, inplace=inplace)

    def transmute(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        raw: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        plan = self._plan()
        if self._axis == "obs":
            retained = {column: col(column) for column in self._spec.columns}
            retained.update(obs or {})
            obs = retained
        else:
            retained = {column: col(column) for column in self._spec.columns}
            retained.update(var or {})
            var = retained
        result = transmute_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            layer=layer,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def mutate(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        raw: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        layer: str | None = None,
        inplace: bool = False,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        plan = self._plan()
        result = mutate_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            layer=layer,
            inplace=inplace,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._same_shape_result(result, inplace=inplace)

    def summarize(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        raw: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        plan = self._plan()
        result = summarize_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            by=list(self._spec.columns),
            layer=layer,
            max_matrix_values=max_matrix_values,
        ).reset_index(drop=True)
        if len(result) == len(plan.keys):
            for column in self._spec.columns:
                result[column] = pd.Series(plan.keys[column].array, index=result.index)
        return result

    summarise = summarize

    def count(self, *, wt: Any = None, sort: bool = False, name: str = "n") -> pd.DataFrame:
        weights = None
        if wt is not None:
            frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
            weights = evaluate_assignments(frame, {"__annplyr_wt__": wt})["__annplyr_wt__"]
        result = GroupPlan.count_for(self._adata, self._spec, weights=weights, name=name)
        if sort:
            result = result.sort_values(name, ascending=False, kind="mergesort").reset_index(drop=True)
        return result

    def tally(self, *, wt: Any = None, sort: bool = False, name: str = "n") -> pd.DataFrame:
        return self.count(wt=wt, sort=sort, name=name)

    def add_count(
        self,
        *,
        wt: Any = None,
        sort: bool = False,
        name: str = "n",
        inplace: bool = False,
    ) -> GroupedAnnData:
        self._plan()
        result = add_count_adata(
            self._adata,
            by=list(self._spec.columns),
            wt=wt,
            sort=sort,
            axis=self._axis,
            name=name,
            inplace=inplace,
        )
        return self._same_shape_result(result, inplace=inplace)

    def add_tally(
        self,
        *,
        wt: Any = None,
        sort: bool = False,
        name: str = "n",
        inplace: bool = False,
    ) -> GroupedAnnData:
        return self.add_count(wt=wt, sort=sort, name=name, inplace=inplace)

    def arrange(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        raw: Any = None,
        obsm: Mapping[str, Any] | None = None,
        varm: Mapping[str, Any] | None = None,
        layer: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        plan = self._plan()
        result = arrange_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            layer=layer,
            copy=copy,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def distinct(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        *,
        keep_all: bool = False,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        if self._axis == "obs" and var is not None:
            raise IncompatibleAxisError("obs-grouped distinct cannot use var selectors")
        if self._axis == "var" and (obs is not None or x is not None):
            raise IncompatibleAxisError("var-grouped distinct cannot use obs or x selectors")
        if not keep_all and self._axis == "obs" and obs is not None:
            obs = self._retain_group_keys(obs, axis="obs")
        if not keep_all and self._axis == "var" and var is not None:
            var = self._retain_group_keys(var, axis="var")
        plan = self._plan()
        result = distinct_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            axis=self._axis,
            keep_all=keep_all,
            copy=copy,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def slice(self, *indices: Any, copy: bool = True) -> GroupedAnnData:
        plan = self._plan()
        selector = _slice_selector(indices)
        pieces: list[np.ndarray] = []
        for positions in plan.positions:
            if isinstance(selector, slice):
                local = np.arange(len(positions), dtype=np.intp)[selector]
            else:
                requested = np.asarray(selector)
                if np.issubdtype(requested.dtype, np.bool_):
                    local = _axis_positions(requested, len(positions), axis=self._axis)
                else:
                    normalized = requested.astype(np.intp, copy=True)
                    normalized[normalized < 0] += len(positions)
                    local = normalized[(normalized >= 0) & (normalized < len(positions))]
            pieces.append(positions[local])
        return self._subset_group_positions(pieces, copy=copy)

    def slice_head(self, n: int = 5, *, copy: bool = True) -> GroupedAnnData:
        _validate_slice_n(n)
        return self.slice(slice(0, n), copy=copy)

    def slice_tail(self, n: int = 5, *, copy: bool = True) -> GroupedAnnData:
        _validate_slice_n(n)
        return self.slice(slice(0, 0) if n == 0 else slice(-n, None), copy=copy)

    def slice_min(self, by: Any, n: int = 5, *, copy: bool = True) -> GroupedAnnData:
        return self._slice_ordered(by, n=n, descending=False, copy=copy)

    def slice_max(self, by: Any, n: int = 5, *, copy: bool = True) -> GroupedAnnData:
        return self._slice_ordered(by, n=n, descending=True, copy=copy)

    def _slice_ordered(self, by: Any, *, n: int, descending: bool, copy: bool) -> GroupedAnnData:
        _validate_slice_n(n)
        plan = self._plan()
        frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
        order_by = _desc_order_by(by) if descending else by
        pieces: list[np.ndarray] = []
        for positions in plan.positions:
            local = frame.iloc[positions, :].reset_index(drop=True)
            ordered = _sort_values_for_frame(local, order_by).to_numpy(dtype=np.intp)[:n]
            pieces.append(positions[ordered])
        return self._subset_group_positions(pieces, copy=copy)

    def slice_sample(
        self,
        n: int | None = None,
        *,
        prop: float | None = None,
        replace: bool = False,
        random_state: int | None = None,
        copy: bool = True,
    ) -> GroupedAnnData:
        if n is not None and prop is not None:
            raise SelectionError("slice_sample accepts n and prop as mutually exclusive arguments")
        if n is not None and n < 0:
            raise SelectionError("slice_sample n must be non-negative")
        if prop is not None and prop < 0:
            raise SelectionError("slice_sample prop must be non-negative")
        plan = self._plan()
        rng = np.random.default_rng(random_state)
        pieces: list[np.ndarray] = []
        for positions in plan.positions:
            take = int(round(len(positions) * prop)) if n is None and prop is not None else (n if n is not None else 1)
            if not replace and take > len(positions):
                raise SelectionError("slice_sample n cannot be larger than a group size unless replace=True")
            pieces.append(positions[rng.choice(len(positions), size=take, replace=replace)])
        return self._subset_group_positions(pieces, copy=copy)

    def _subset_group_positions(self, pieces: list[np.ndarray], *, copy: bool) -> GroupedAnnData:
        positions = np.concatenate(pieces) if pieces else np.empty(0, dtype=np.intp)
        result = (
            _subset_positions(self._adata, positions, None, copy=copy)
            if self._axis == "obs"
            else _subset_positions(self._adata, None, positions, copy=copy)
        )
        return self._wrap(result)

    def left_join(
        self,
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
    ) -> GroupedAnnData:
        return self._join(
            "left",
            other,
            by=by,
            axis=axis,
            relationship=relationship,
            multiple=multiple,
            unmatched=unmatched,
            na_matches=na_matches,
            suffixes=suffixes,
            copy=copy,
        )

    def inner_join(
        self,
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
    ) -> GroupedAnnData:
        return self._join(
            "inner",
            other,
            by=by,
            axis=axis,
            relationship=relationship,
            multiple=multiple,
            unmatched=unmatched,
            na_matches=na_matches,
            suffixes=suffixes,
            copy=copy,
        )

    def right_join(
        self,
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
    ) -> GroupedAnnData:
        return self._join(
            "right",
            other,
            by=by,
            axis=axis,
            relationship=relationship,
            multiple=multiple,
            unmatched=unmatched,
            na_matches=na_matches,
            suffixes=suffixes,
            copy=copy,
        )

    def full_join(
        self,
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
    ) -> GroupedAnnData:
        return self._join(
            "full",
            other,
            by=by,
            axis=axis,
            relationship=relationship,
            multiple=multiple,
            unmatched=unmatched,
            na_matches=na_matches,
            suffixes=suffixes,
            copy=copy,
        )

    def semi_join(
        self,
        other: pd.DataFrame | Mapping[str, Any],
        *,
        by: str | Sequence[str] | None = None,
        axis: str = "obs",
        na_matches: str = "na",
        copy: bool = True,
    ) -> GroupedAnnData:
        return self._join("semi", other, by=by, axis=axis, na_matches=na_matches, copy=copy)

    def anti_join(
        self,
        other: pd.DataFrame | Mapping[str, Any],
        *,
        by: str | Sequence[str] | None = None,
        axis: str = "obs",
        na_matches: str = "na",
        copy: bool = True,
    ) -> GroupedAnnData:
        return self._join("anti", other, by=by, axis=axis, na_matches=na_matches, copy=copy)

    def _join(
        self,
        kind: str,
        other: pd.DataFrame | Mapping[str, Any],
        *,
        by: str | Sequence[str] | None = None,
        axis: str = "obs",
        relationship: str = "many-to-one",
        multiple: str = "error",
        unmatched: str | None = None,
        na_matches: str = "na",
        suffixes: tuple[str, str] = ("", "_right"),
        copy: bool = True,
    ) -> GroupedAnnData:
        self._plan()
        functions = {
            "left": left_join_adata,
            "inner": inner_join_adata,
            "right": right_join_adata,
            "full": full_join_adata,
            "semi": semi_join_adata,
            "anti": anti_join_adata,
        }
        kwargs: dict[str, Any] = {"by": by, "axis": axis, "na_matches": na_matches, "copy": copy}
        if kind not in {"semi", "anti"}:
            kwargs.update(
                relationship=relationship,
                multiple=multiple,
                unmatched=unmatched if unmatched is not None else ("error" if kind in {"right", "full"} else "drop"),
                suffixes=suffixes,
            )
        result = functions[kind](self._adata, other, **kwargs)
        spec = (
            self._join_spec(other, by=by, axis=axis, suffixes=suffixes, result=result)
            if kind not in {"semi", "anti"}
            else self._spec
        )
        spec.validate(result)
        return self._wrap(result, spec=spec)

    def _join_spec(
        self,
        other: pd.DataFrame | Mapping[str, Any],
        *,
        by: str | Sequence[str] | None,
        axis: str,
        suffixes: tuple[str, str],
        result: AnnData,
    ) -> GroupSpec:
        if axis != self._axis:
            return self._spec
        right = other if isinstance(other, pd.DataFrame) else pd.DataFrame(other)
        left = cast(pd.DataFrame, self._adata.obs if axis == "obs" else self._adata.var)
        by_columns = (
            [column for column in left.columns if column in right.columns]
            if by is None
            else ([by] if isinstance(by, str) else list(by))
        )
        update = {
            column: f"{column}{suffixes[0]}"
            for column in self._spec.columns
            if column in right.columns and column not in by_columns
        }
        spec = self._spec.renamed(update)
        spec.validate(result)
        return spec

    def pipe(self, func: Callable[..., T] | tuple[Callable[..., T], str], *args: Any, **kwargs: Any) -> T:
        if isinstance(func, tuple):
            call, data_keyword = func
            if data_keyword in kwargs:
                raise AnnplyrError(f"{data_keyword!r} is both the pipe target and a keyword argument")
            kwargs[data_keyword] = self
            return call(*args, **kwargs)
        return func(self, *args, **kwargs)
