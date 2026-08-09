from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any, TypeVar, cast, overload

import pandas as pd
from anndata import AnnData

from annplyr._errors import AnnplyrError
from annplyr._extensions import register_anndata_accessor
from annplyr._grouped import GroupedAnnData
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
    Source,
    SourceSelectors,
)
from annplyr._verbs import (
    add_count_adata,
    add_tally_adata,
    anti_join_adata,
    arrange_adata,
    as_frame_adata,
    count_adata,
    distinct_adata,
    filter_adata,
    full_join_adata,
    inner_join_adata,
    left_join_adata,
    mutate_adata,
    nest_by_adata,
    pivot_longer_adata,
    pull_adata,
    relocate_adata,
    rename_adata,
    rename_with_adata,
    right_join_adata,
    select_adata,
    semi_join_adata,
    slice_adata,
    slice_head_adata,
    slice_max_adata,
    slice_min_adata,
    slice_sample_adata,
    slice_tail_adata,
    summarize_adata,
    tally_adata,
    to_df_adata,
    to_tidy_adata,
    transmute_adata,
)

T = TypeVar("T")


def _with_annplyr_accessor(adata: AnnData) -> AnnDataWithAnnplyr:
    """Narrow an AnnData result after the namespace has been registered."""
    return cast(AnnDataWithAnnplyr, adata)


@register_anndata_accessor("ap")
class AnnplyrAccessor:
    """Dataframe-style AnnData wrangling accessor."""

    def __init__(self, adata: AnnData):
        self._obj = cast(AnnDataWithAnnplyr, adata)

    def filter(
        self,
        obs: Expression | Sequence[Expression] | None = None,
        var: Expression | Sequence[Expression] | None = None,
        x: Expression | Sequence[Expression] | None = None,
        raw: Expression | Sequence[Expression] | None = None,
        obs_names: Expression | Sequence[Expression] | None = None,
        var_names: Expression | Sequence[Expression] | None = None,
        obsm: SourceSelectors | None = None,
        varm: SourceSelectors | None = None,
        layer: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            filter_adata(
                self._obj,
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
            )
        )

    def select(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(select_adata(self._obj, obs=obs, var=var, x=x, copy=copy))

    def rename(
        self,
        obs: Mapping[str, str] | None = None,
        var: Mapping[str, str] | None = None,
        x: Mapping[str, str] | None = None,
        *,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(rename_adata(self._obj, obs=obs, var=var, x=x, inplace=inplace))

    def rename_with(
        self,
        func: Callable[[str], str],
        *,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(rename_with_adata(self._obj, func, obs=obs, var=var, x=x, inplace=inplace))

    def relocate(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        *,
        before: str | None = None,
        after: str | None = None,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            relocate_adata(self._obj, obs=obs, var=var, x=x, before=before, after=after, inplace=inplace)
        )

    def distinct(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        *,
        axis: Axis = "obs",
        keep_all: bool = False,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            distinct_adata(
                self._obj,
                obs=obs,
                var=var,
                x=x,
                axis=axis,
                keep_all=keep_all,
                copy=copy,
                max_matrix_values=max_matrix_values,
            )
        )

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
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            left_join_adata(
                self._obj,
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
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            inner_join_adata(
                self._obj,
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
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            right_join_adata(
                self._obj,
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
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            full_join_adata(
                self._obj,
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
        )

    def semi_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        na_matches: NaMatches = "na",
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            semi_join_adata(self._obj, other, by=by, axis=axis, na_matches=na_matches, copy=copy)
        )

    def anti_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        na_matches: NaMatches = "na",
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            anti_join_adata(self._obj, other, by=by, axis=axis, na_matches=na_matches, copy=copy)
        )

    def arrange(
        self,
        obs: Expression | Sequence[Expression] | None = None,
        var: Expression | Sequence[Expression] | None = None,
        x: Expression | Sequence[Expression] | None = None,
        raw: Expression | Sequence[Expression] | None = None,
        obsm: SourceSelectors | None = None,
        varm: SourceSelectors | None = None,
        layer: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            arrange_adata(
                self._obj,
                obs=obs,
                var=var,
                x=x,
                raw=raw,
                obsm=obsm,
                varm=varm,
                layer=layer,
                copy=copy,
                max_matrix_values=max_matrix_values,
            )
        )

    def slice(self, *indices: Any, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_adata(self._obj, *indices, axis=axis, copy=copy))

    def slice_head(self, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_head_adata(self._obj, n=n, axis=axis, copy=copy))

    def slice_tail(self, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_tail_adata(self._obj, n=n, axis=axis, copy=copy))

    def slice_min(self, by: Expression, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_min_adata(self._obj, by=by, n=n, axis=axis, copy=copy))

    def slice_max(self, by: Expression, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_max_adata(self._obj, by=by, n=n, axis=axis, copy=copy))

    def slice_sample(
        self,
        n: int | None = None,
        *,
        prop: float | None = None,
        replace: bool = False,
        random_state: int | None = None,
        axis: Axis = "obs",
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            slice_sample_adata(
                self._obj,
                n=n,
                prop=prop,
                replace=replace,
                random_state=random_state,
                axis=axis,
                copy=copy,
            )
        )

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
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            mutate_adata(
                self._obj,
                obs=obs,
                var=var,
                x=x,
                raw=raw,
                obsm=obsm,
                varm=varm,
                layer=layer,
                inplace=inplace,
                max_matrix_values=max_matrix_values,
            )
        )

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
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            transmute_adata(
                self._obj,
                obs=obs,
                var=var,
                x=x,
                raw=raw,
                obsm=obsm,
                varm=varm,
                layer=layer,
                max_matrix_values=max_matrix_values,
            )
        )

    @overload
    def group_by(self, obs: None = None, var: None = None) -> AnnDataWithAnnplyr: ...

    @overload
    def group_by(self, obs: Selector, var: None = None) -> GroupedAnnData: ...

    @overload
    def group_by(self, obs: None = None, var: Selector = ...) -> GroupedAnnData: ...

    def group_by(self, obs: Selector | None = None, var: Selector | None = None) -> AnnDataWithAnnplyr | GroupedAnnData:
        if obs is None and var is None:
            return self._obj
        return GroupedAnnData(self._obj, obs=obs, var=var)

    def summarize(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        raw: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        *,
        by: Selector | None = None,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return summarize_adata(
            self._obj,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            by=by,
            layer=layer,
            max_matrix_values=max_matrix_values,
        )

    summarise = summarize

    def count(
        self,
        by: Selector | None = None,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
    ) -> pd.DataFrame:
        return count_adata(self._obj, by=by, wt=wt, sort=sort, axis=axis, name=name)

    def tally(
        self,
        by: Selector | None = None,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
    ) -> pd.DataFrame:
        return tally_adata(self._obj, by=by, wt=wt, sort=sort, axis=axis, name=name)

    def add_count(
        self,
        by: Selector | None = None,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            add_count_adata(self._obj, by=by, wt=wt, sort=sort, axis=axis, name=name, inplace=inplace)
        )

    def add_tally(
        self,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            add_tally_adata(self._obj, wt=wt, sort=sort, axis=axis, name=name, inplace=inplace)
        )

    def pull(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        obsm: SourceSelectors | None = None,
        varm: SourceSelectors | None = None,
        obsp: SourceSelectors | None = None,
        varp: SourceSelectors | None = None,
        uns: SourceSelectors | None = None,
        *,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.Series:
        return pull_adata(
            self._obj,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            obsp=obsp,
            varp=varp,
            uns=uns,
            layer=layer,
            max_matrix_values=max_matrix_values,
        )

    def to_df(
        self,
        obs: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        obsm: SourceSelectors | None = None,
        obsp: SourceSelectors | None = None,
        *,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return to_df_adata(
            self._obj,
            obs=obs,
            x=x,
            raw=raw,
            obsm=obsm,
            obsp=obsp,
            layer=layer,
            max_matrix_values=max_matrix_values,
        )

    def to_tidy(
        self,
        obs: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        *,
        layer: str | None = None,
        obs_name: str = "obs_name",
        feature: str = "feature",
        value: str = "value",
        allow_all_features: bool = False,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return to_tidy_adata(
            self._obj,
            obs=obs,
            x=x,
            raw=raw,
            layer=layer,
            obs_name=obs_name,
            feature=feature,
            value=value,
            allow_all_features=allow_all_features,
            max_matrix_values=max_matrix_values,
        )

    def pivot_longer(
        self,
        obs: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        *,
        layer: str | None = None,
        obs_name: str = "obs_name",
        names_to: str = "name",
        values_to: str = "value",
        allow_all_features: bool = False,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return pivot_longer_adata(
            self._obj,
            obs=obs,
            x=x,
            raw=raw,
            layer=layer,
            obs_name=obs_name,
            names_to=names_to,
            values_to=values_to,
            allow_all_features=allow_all_features,
            max_matrix_values=max_matrix_values,
        )

    def as_frame(
        self,
        source: Source,
        *,
        key: str | None = None,
        select: Selector | None = None,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return as_frame_adata(
            self._obj, source=source, key=key, select=select, layer=layer, max_matrix_values=max_matrix_values
        )

    def nest_by(
        self,
        *,
        by: Selector,
        obs: Selector | None = None,
        var: Selector | None = None,
        axis: Axis = "obs",
        name: str = "data",
    ) -> pd.DataFrame:
        return nest_by_adata(self._obj, by=by, obs=obs, var=var, axis=axis, name=name)

    def pipe(self, func: Callable[..., T] | tuple[Callable[..., T], str], *args: Any, **kwargs: Any) -> T:
        if isinstance(func, tuple):
            call, data_keyword = func
            if data_keyword in kwargs:
                msg = f"{data_keyword!r} is both the pipe target and a keyword argument"
                raise AnnplyrError(msg)
            kwargs[data_keyword] = self._obj
            return call(*args, **kwargs)
        return func(self._obj, *args, **kwargs)
