from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any, TypeVar

from anndata import AnnData

from annplyr._extensions import register_anndata_accessor
from annplyr._grouped import GroupedAnnData
from annplyr._verbs import (
    arrange_adata,
    count_adata,
    filter_adata,
    mutate_adata,
    pull_adata,
    select_adata,
    slice_adata,
    slice_head_adata,
    slice_max_adata,
    slice_min_adata,
    slice_sample_adata,
    slice_tail_adata,
    summarize_adata,
    to_df_adata,
    to_tidy_adata,
)

T = TypeVar("T")


@register_anndata_accessor("ap")
class AnnplyrAccessor:
    """Dataframe-style AnnData wrangling accessor."""

    def __init__(self, adata: AnnData):
        self._obj = adata

    def filter(
        self,
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
        return filter_adata(self._obj, obs, var, x, obs_names, var_names, obsm, varm, layer, copy)

    def select(self, obs: Any = None, var: Any = None, x: Any = None, copy: bool = False) -> AnnData:
        return select_adata(self._obj, obs=obs, var=var, x=x, copy=copy)

    def arrange(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        obsm: Mapping[str, Any] | None = None,
        varm: Mapping[str, Any] | None = None,
        layer: str | None = None,
        copy: bool = False,
    ) -> AnnData:
        return arrange_adata(self._obj, obs=obs, var=var, x=x, obsm=obsm, varm=varm, layer=layer, copy=copy)

    def slice(self, *indices: Any, axis: str = "obs", copy: bool = False) -> AnnData:
        return slice_adata(self._obj, *indices, axis=axis, copy=copy)

    def slice_head(self, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
        return slice_head_adata(self._obj, n=n, axis=axis, copy=copy)

    def slice_tail(self, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
        return slice_tail_adata(self._obj, n=n, axis=axis, copy=copy)

    def slice_min(self, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
        return slice_min_adata(self._obj, by=by, n=n, axis=axis, copy=copy)

    def slice_max(self, by: Any, n: int = 5, *, axis: str = "obs", copy: bool = False) -> AnnData:
        return slice_max_adata(self._obj, by=by, n=n, axis=axis, copy=copy)

    def slice_sample(
        self,
        n: int | None = None,
        *,
        prop: float | None = None,
        replace: bool = False,
        random_state: int | None = None,
        axis: str = "obs",
        copy: bool = False,
    ) -> AnnData:
        return slice_sample_adata(
            self._obj,
            n=n,
            prop=prop,
            replace=replace,
            random_state=random_state,
            axis=axis,
            copy=copy,
        )

    def mutate(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        layer: str | None = None,
        inplace: bool = False,
    ) -> AnnData:
        return mutate_adata(self._obj, obs=obs, var=var, x=x, obsm=obsm, varm=varm, layer=layer, inplace=inplace)

    def group_by(self, obs: Any = None, var: Any = None) -> AnnData | GroupedAnnData:
        if obs is None and var is None:
            return self._obj
        return GroupedAnnData(self._obj, obs=obs, var=var)

    def summarize(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        *,
        by: Any = None,
        layer: str | None = None,
    ):
        return summarize_adata(self._obj, obs=obs, var=var, x=x, obsm=obsm, varm=varm, by=by, layer=layer)

    summarise = summarize

    def count(self, by: Any = None, *, axis: str = "obs", name: str = "n"):
        return count_adata(self._obj, by=by, axis=axis, name=name)

    def pull(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        obsm: Mapping[str, Any] | None = None,
        varm: Mapping[str, Any] | None = None,
        *,
        layer: str | None = None,
    ):
        return pull_adata(self._obj, obs=obs, var=var, x=x, obsm=obsm, varm=varm, layer=layer)

    def to_df(
        self,
        obs: Any = None,
        x: Any = None,
        obsm: Mapping[str, Any] | None = None,
        *,
        layer: str | None = None,
    ):
        return to_df_adata(self._obj, obs=obs, x=x, obsm=obsm, layer=layer)

    def to_tidy(
        self,
        obs: Any = None,
        x: Any = None,
        *,
        layer: str | None = None,
        obs_name: str = "obs_name",
        feature: str = "feature",
        value: str = "value",
    ):
        return to_tidy_adata(
            self._obj,
            obs=obs,
            x=x,
            layer=layer,
            obs_name=obs_name,
            feature=feature,
            value=value,
        )

    def pipe(self, func: Callable[..., T] | tuple[Callable[..., T], str], *args: Any, **kwargs: Any) -> T:
        if isinstance(func, tuple):
            call, data_keyword = func
            if data_keyword in kwargs:
                msg = f"{data_keyword!r} is both the pipe target and a keyword argument"
                raise ValueError(msg)
            kwargs[data_keyword] = self._obj
            return call(*args, **kwargs)
        return func(self._obj, *args, **kwargs)
