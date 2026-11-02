from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any, Literal, TypeVar, cast

import numpy as np
import pandas as pd
from anndata import AnnData

from annplyr._aggregate import aggregate_cells_adata
from annplyr._errors import AnnplyrError, IncompatibleAxisError, SelectionError
from annplyr._expr import AnnplyrExpr, col, everything
from annplyr._frames import evaluate_assignments, evaluate_select, obs_frame, var_frame
from annplyr._groups import GroupPlan, GroupSpec
from annplyr._resolve import META_AND_X, bare_key, container_kwargs, relocate_anchor, rename_kwargs
from annplyr._tidy import TidyAnnData
from annplyr._typing import (
    AnnDataWithAnnplyr,
    Axis,
    Expression,
    JoinBy,
    JoinInput,
    JoinMultiple,
    JoinRelationship,
    JoinUnmatched,
    NaMatches,
    Selector,
)
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
    merge_named,
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
        obs: Selector | None = None,
        var: Selector | None = None,
        spec: GroupSpec | None = None,
    ):
        self._adata = cast(AnnDataWithAnnplyr, adata)
        self._spec = spec if spec is not None else GroupSpec.resolve(adata, obs=obs, var=var)
        self._spec.validate(adata)

    @property
    def _axis(self) -> Axis:
        return self._spec.axis

    @property
    def _by(self) -> list[str]:
        return list(self._spec.columns)

    def _plan(self) -> GroupPlan:
        return GroupPlan.build(self._adata, self._spec)

    def _wrap(self, adata: AnnData, *, spec: GroupSpec | None = None) -> GroupedAnnData:
        return GroupedAnnData(adata, spec=spec or self._spec)

    def __iter__(self) -> Iterator[tuple[dict[str, Any], AnnDataWithAnnplyr]]:
        plan = self._plan()
        for row, positions in zip(range(len(plan.keys)), plan.positions, strict=True):
            key = {column: plan.keys.iloc[row][column] for column in self._spec.columns}
            if self._axis == "obs":
                group = self._adata[positions]
            else:
                group = self._adata[:, positions]
            yield key, cast(AnnDataWithAnnplyr, TidyAnnData.wrap(group))

    def group_vars(self) -> list[str]:
        self._spec.validate(self._adata)
        return list(self._spec.columns)

    def group_keys(self) -> pd.DataFrame:
        return GroupPlan.keys_for(self._adata, self._spec)

    def group_data(self) -> pd.DataFrame:
        return self._plan().group_data()

    def ungroup(self) -> AnnDataWithAnnplyr:
        return cast(AnnDataWithAnnplyr, TidyAnnData.wrap(self._adata))

    def filter(
        self,
        *predicates: Expression | Sequence[Expression],
        source: str | None = None,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        """Subset within groups, one predicate per argument, combined with AND.

        Same container syntax as :meth:`annplyr.AnnplyrAccessor.filter`; the
        predicates are applied inside each group, and groups are concatenated
        in first-seen order.

        Parameters
        ----------
        *predicates
            See the annotated public signature.
        max_matrix_values
            See the annotated public signature.

        Returns
        -------
        GroupedAnnData

        Ownership
        ---------
        Returns an independent object; a backed input returns the view.

        Raises
        ------
        SelectionError
            If the request violates the corresponding public contract.
        UnknownColumnError
            If the request violates the corresponding public contract.
        UnknownSourceError
            If the request violates the corresponding public contract.
        IncompatibleAxisError
            If the request violates the corresponding public contract.
        AnnplyrError
            If the request violates the corresponding public contract.
        """
        plan = self._plan()
        result = filter_adata(
            self._adata,
            *predicates,
            source=source,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def select(self, *names: Selector) -> GroupedAnnData:
        self._spec.validate(self._adata)
        kwargs = container_kwargs(self._adata, names, verb="select", allowed=META_AND_X, negation=True)
        if self._axis in kwargs:
            kwargs[self._axis] = self._retain_group_keys(kwargs[self._axis], axis=self._axis)
        return self._wrap(select_adata(self._adata, **kwargs))

    def _retain_group_keys(self, selector: Selector, *, axis: Axis) -> list[str]:
        frame = obs_frame(self._adata) if axis == "obs" else var_frame(self._adata)
        selected = [str(column) for column in evaluate_select(frame, selector).columns]
        omitted = [column for column in self._spec.columns if column not in selected]
        return [*omitted, *selected]

    def rename(self, mapping: Mapping[str, Any] | None = None, /, **new_to_old: Any) -> GroupedAnnData:
        kwargs = rename_kwargs(self._adata, mapping, new_to_old)
        update = {old: new for new, old in kwargs.get(self._axis, {}).items()}
        result = rename_adata(self._adata, **kwargs)
        return self._wrap(result, spec=self._spec.renamed(update))

    def rename_with(
        self,
        func: Callable[[str], str],
        *names: Selector,
    ) -> GroupedAnnData:
        before = list(cast(pd.DataFrame, self._adata.obs if self._axis == "obs" else self._adata.var).columns)
        kwargs = container_kwargs(self._adata, names, verb="rename_with", allowed=META_AND_X) or {"obs": everything()}
        result = rename_with_adata(self._adata, func, **kwargs)
        after = list(cast(pd.DataFrame, result.obs if self._axis == "obs" else result.var).columns)
        spec = self._spec.renamed(dict(zip(before, after, strict=True)))
        return self._wrap(result, spec=spec)

    def relocate(
        self,
        *names: Selector,
        before: str | AnnplyrExpr | None = None,
        after: str | AnnplyrExpr | None = None,
    ) -> GroupedAnnData:
        kwargs = container_kwargs(self._adata, names, verb="relocate", allowed=META_AND_X, negation=True)
        result = relocate_adata(
            self._adata,
            **kwargs,
            before=relocate_anchor(self._adata, before, kwargs),
            after=relocate_anchor(self._adata, after, kwargs),
        )
        return self._wrap(result)

    def transmute(
        self,
        assignments: Mapping[str, Any] | None = None,
        /,
        *,
        into: str = "obs",
        max_matrix_values: int | None = None,
        **named: Any,
    ) -> GroupedAnnData:
        assignments = merge_named(assignments, named, verb="transmute")
        plan = self._plan()
        # group keys are always retained so future grouped calls keep working
        destination = into if into in {"obs", "var"} else self._axis
        if destination == self._axis:
            retained: dict[str, Any] = {column: col(column) for column in self._spec.columns}
            retained.update(assignments or {})
            assignments = retained
        result = transmute_adata(
            self._adata,
            assignments,
            into=destination,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def mutate(
        self,
        assignments: Mapping[str, Any] | None = None,
        /,
        *,
        into: str = "obs",
        max_matrix_values: int | None = None,
        **named: Any,
    ) -> GroupedAnnData:
        plan = self._plan()
        result = mutate_adata(
            self._adata,
            merge_named(assignments, named, verb="mutate"),
            into=into,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def summarize(
        self,
        stats: Mapping[str, Any] | None = None,
        /,
        *,
        source: str | None = None,
        max_matrix_values: int | None = None,
        **named: Any,
    ) -> pd.DataFrame:
        if "by" in named:
            msg = "group_by already sets the keys, so the grouped summarize takes no by="
            raise TypeError(msg)
        plan = self._plan()
        result = summarize_adata(
            self._adata,
            merge_named(stats, named, verb="summarize"),
            by=list(self._spec.columns),
            source=source,
            max_matrix_values=max_matrix_values,
            _group_plan=plan,
        ).reset_index(drop=True)
        if len(result) == len(plan.keys):
            for column in self._spec.columns:
                result[column] = pd.Series(plan.keys[column].array, index=result.index)
        return result

    summarise = summarize

    def count(self, *, wt: Expression | None = None, sort: bool = False, name: str = "n") -> pd.DataFrame:
        weights = None
        if wt is not None:
            frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
            weights = evaluate_assignments(frame, {"__annplyr_wt__": wt})["__annplyr_wt__"]
        result = GroupPlan.count_for(self._adata, self._spec, weights=weights, name=name)
        if sort:
            result = result.sort_values(name, ascending=False, kind="mergesort").reset_index(drop=True)
        return result

    def add_count(
        self,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        name: str = "n",
    ) -> GroupedAnnData:
        self._plan()
        result = add_count_adata(
            self._adata,
            by=list(self._spec.columns),
            wt=wt,
            sort=sort,
            axis=self._axis,
            name=name,
        )
        return self._wrap(result)

    def arrange(
        self,
        *by: Expression | Sequence[Expression],
        source: str | None = None,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        plan = self._plan()
        return self._wrap(
            arrange_adata(
                self._adata,
                *by,
                source=source,
                max_matrix_values=max_matrix_values,
                _group_positions=plan.positions,
                _group_axis=self._axis,
            )
        )

    def distinct(
        self,
        *names: Selector,
        keep_all: bool = False,
        max_matrix_values: int | None = None,
    ) -> GroupedAnnData:
        kwargs = container_kwargs(self._adata, names, verb="distinct", allowed=META_AND_X, negation=True)
        wrong = {"obs": ("var", "x"), "var": ("obs", "x")}[self._axis]
        if any(kind in kwargs for kind in wrong):
            raise IncompatibleAxisError(f"{self._axis}-grouped distinct cannot use {' or '.join(wrong)} names")
        if not keep_all and self._axis in kwargs:
            kwargs[self._axis] = self._retain_group_keys(kwargs[self._axis], axis=self._axis)
        plan = self._plan()
        result = distinct_adata(
            self._adata,
            **kwargs,
            axis=self._axis,
            keep_all=keep_all,
            max_matrix_values=max_matrix_values,
            _group_positions=plan.positions,
            _group_axis=self._axis,
        )
        return self._wrap(result)

    def slice(self, *indices: Any) -> GroupedAnnData:
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
        return self._subset_group_positions(pieces)

    def slice_head(self, n: int = 5) -> GroupedAnnData:
        _validate_slice_n(n)
        return self.slice(slice(0, n))

    def slice_tail(self, n: int = 5) -> GroupedAnnData:
        _validate_slice_n(n)
        return self.slice(slice(0, 0) if n == 0 else slice(-n, None))

    def slice_min(self, by: Expression, n: int = 5) -> GroupedAnnData:
        return self._slice_ordered(by, n=n, descending=False)

    def slice_max(self, by: Expression, n: int = 5) -> GroupedAnnData:
        return self._slice_ordered(by, n=n, descending=True)

    def _slice_ordered(self, by: Expression, *, n: int, descending: bool) -> GroupedAnnData:
        _validate_slice_n(n)
        plan = self._plan()
        frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
        by = bare_key(by, axis=self._axis)
        order_by = _desc_order_by(by) if descending else by
        pieces: list[np.ndarray] = []
        for positions in plan.positions:
            local = frame.iloc[positions, :].reset_index(drop=True)
            ordered = _sort_values_for_frame(local, order_by).to_numpy(dtype=np.intp)[:n]
            pieces.append(positions[ordered])
        return self._subset_group_positions(pieces)

    def slice_sample(
        self,
        n: int | None = None,
        *,
        prop: float | None = None,
        replace: bool = False,
        random_state: int | None = None,
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
        return self._subset_group_positions(pieces)

    def _subset_group_positions(self, pieces: list[np.ndarray]) -> GroupedAnnData:
        positions = np.concatenate(pieces) if pieces else np.empty(0, dtype=np.intp)
        result = (
            _subset_positions(self._adata, positions, None)
            if self._axis == "obs"
            else _subset_positions(self._adata, None, positions)
        )
        return self._wrap(result)

    def left_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "drop",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
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
        )

    def inner_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "drop",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
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
        )

    def right_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "error",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
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
        )

    def full_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "error",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
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
        )

    def semi_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        na_matches: NaMatches = "na",
    ) -> GroupedAnnData:
        return self._join("semi", other, by=by, axis=axis, na_matches=na_matches)

    def anti_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        na_matches: NaMatches = "na",
    ) -> GroupedAnnData:
        return self._join("anti", other, by=by, axis=axis, na_matches=na_matches)

    def _join(
        self,
        kind: str,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched | None = None,
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
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
        kwargs: dict[str, Any] = {"by": by, "axis": axis, "na_matches": na_matches}
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
        other: JoinInput,
        *,
        by: JoinBy,
        axis: Axis,
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

    def aggregate_cells(self, *, source: Any, fun: Literal["sum", "mean"] = "sum") -> TidyAnnData:
        """Pseudobulk over the grouping keys; see :meth:`AnnplyrAccessor.aggregate_cells`."""
        if self._axis != "obs":
            raise IncompatibleAxisError("aggregate_cells groups cells, so group_by obs metadata")
        return TidyAnnData.wrap(aggregate_cells_adata(self._adata, *self._spec.columns, source=source, fun=fun))

    def pipe(self, func: Callable[..., T] | tuple[Callable[..., T], str], *args: Any, **kwargs: Any) -> T:
        if isinstance(func, tuple):
            call, data_keyword = func
            if data_keyword in kwargs:
                raise AnnplyrError(f"{data_keyword!r} is both the pipe target and a keyword argument")
            kwargs[data_keyword] = self
            return call(*args, **kwargs)
        return func(self, *args, **kwargs)
