from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any

import pandas as pd
from anndata import AnnData

from annplyr._frames import evaluate_select, obs_frame, var_frame
from annplyr._verbs import count_frame, filter_adata, mutate_adata, prepare_by_frame, summarize_frame


class GroupedAnnData:
    """A grouped AnnData wrapper."""

    def __init__(self, adata: AnnData, *, obs: Any = None, var: Any = None):
        if obs is not None and var is not None:
            msg = "group_by currently accepts one axis at a time"
            raise ValueError(msg)
        self._adata = adata
        self._axis = "obs" if obs is not None else "var"
        self._by = obs if obs is not None else var

    def __iter__(self) -> Iterator[tuple[dict[str, Any], AnnData]]:
        frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
        by_frame = evaluate_select(frame, self._by)
        keys = by_frame.drop_duplicates()
        for _, key_row in keys.iterrows():
            mask = pd.Series(True, index=by_frame.index)
            for column, value in key_row.items():
                mask &= by_frame[column].eq(value)
            key = key_row.to_dict()
            if self._axis == "obs":
                yield key, self._adata[mask.to_numpy(), :]
            else:
                yield key, self._adata[:, mask.to_numpy()]

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
        return filter_adata(self._adata, obs, var, x, obs_names, var_names, obsm, varm, layer, copy)

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
        return mutate_adata(self._adata, obs=obs, var=var, x=x, obsm=obsm, varm=varm, layer=layer, inplace=inplace)

    def summarize(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
    ) -> pd.DataFrame:
        frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
        assignments = obs if self._axis == "obs" else var
        work, by_columns = prepare_by_frame(frame, self._by)
        return summarize_frame(work, assignments=assignments, by=by_columns)

    summarise = summarize

    def count(self, *, name: str = "n") -> pd.DataFrame:
        frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
        return count_frame(frame, by=self._by, name=name)
