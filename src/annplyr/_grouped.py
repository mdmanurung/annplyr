from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any, cast

import pandas as pd
from anndata import AnnData

from annplyr._errors import IncompatibleAxisError
from annplyr._frames import evaluate_assignments, evaluate_select, obs_frame, var_frame
from annplyr._verbs import count_frame, filter_adata, mutate_adata, prepare_by_frame, summarize_frame


class GroupedAnnData:
    """A grouped AnnData wrapper."""

    def __init__(self, adata: AnnData, *, obs: Any = None, var: Any = None):
        if obs is not None and var is not None:
            msg = "group_by currently accepts one axis at a time"
            raise IncompatibleAxisError(msg)
        self._adata = adata
        self._axis = "obs" if obs is not None else "var"
        self._by = obs if obs is not None else var

    def __iter__(self) -> Iterator[tuple[dict[str, Any], AnnData]]:
        yield from self._iter_groups(self._adata)

    def _iter_groups(self, adata: AnnData) -> Iterator[tuple[dict[str, Any], AnnData]]:
        by_frame = self._by_frame(adata)
        for _, key_row in by_frame.drop_duplicates().iterrows():
            mask = _key_mask(by_frame, key_row)
            key = _key_dict(key_row)
            if self._axis == "obs":
                yield key, adata[mask.to_numpy(), :]
            else:
                yield key, adata[:, mask.to_numpy()]

    def _axis_frame(self, adata: AnnData) -> pd.DataFrame:
        return obs_frame(adata) if self._axis == "obs" else var_frame(adata)

    def _by_frame(self, adata: AnnData) -> pd.DataFrame:
        return evaluate_select(self._axis_frame(adata), self._by)

    def group_vars(self) -> list[str]:
        """Return grouping variable names."""
        return list(self._by_frame(self._adata).columns)

    def group_keys(self) -> pd.DataFrame:
        """Return one row per group key."""
        return self._by_frame(self._adata).drop_duplicates().reset_index(drop=True)

    def group_data(self) -> pd.DataFrame:
        """Return group keys plus axis labels in each group."""
        by_frame = self._by_frame(self._adata)
        rows: list[dict[str, Any]] = []
        for _, key_row in by_frame.drop_duplicates().iterrows():
            mask = _key_mask(by_frame, key_row)
            row = _key_dict(key_row)
            row[".rows"] = by_frame.index[mask.to_numpy()].tolist()
            rows.append(row)
        return pd.DataFrame(rows, columns=[*by_frame.columns, ".rows"])

    def ungroup(self) -> AnnData:
        """Return the underlying AnnData object."""
        return self._adata

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
        if self._axis == "obs" and obs:
            out = self._adata if inplace else self._adata.copy()
            for _, group in self._iter_groups(out):
                values = evaluate_assignments(obs_frame(group), obs)
                for column in values.columns:
                    cast(pd.DataFrame, out.obs).loc[group.obs_names, column] = values[column].to_numpy()
            return mutate_adata(out, var=var, x=x, obsm=obsm, varm=varm, layer=layer, inplace=True)
        if self._axis == "var" and var:
            out = self._adata if inplace else self._adata.copy()
            for _, group in self._iter_groups(out):
                values = evaluate_assignments(var_frame(group), var)
                for column in values.columns:
                    cast(pd.DataFrame, out.var).loc[group.var_names, column] = values[column].to_numpy()
            return mutate_adata(out, obs=obs, x=x, obsm=obsm, varm=varm, layer=layer, inplace=True)
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


def _key_mask(by_frame: pd.DataFrame, key_row: pd.Series) -> pd.Series:
    mask = pd.Series(True, index=by_frame.index)
    for column, value in key_row.items():
        if pd.isna(value):
            mask &= by_frame[column].isna()
        else:
            mask &= by_frame[column].eq(value)
    return mask


def _key_dict(key_row: pd.Series) -> dict[str, Any]:
    return {str(column): value for column, value in key_row.items()}
