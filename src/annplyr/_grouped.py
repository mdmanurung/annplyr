from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any, cast

import numpy as np
import pandas as pd
from anndata import AnnData

from annplyr._errors import AnnplyrError, IncompatibleAxisError, SelectionError
from annplyr._frames import evaluate_assignments, evaluate_select, obs_frame, var_frame
from annplyr._verbs import (
    _obs_assignment_frames,
    _validate_slice_n,
    _var_assignment_frames,
    add_count_adata,
    arrange_adata,
    count_frame,
    distinct_adata,
    filter_adata,
    mutate_adata,
    slice_max_adata,
    slice_min_adata,
    summarize_adata,
)


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
        raw: Any = None,
        obs_names: Any = None,
        var_names: Any = None,
        obsm: Mapping[str, Any] | None = None,
        varm: Mapping[str, Any] | None = None,
        layer: str | None = None,
        copy: bool = False,
    ) -> AnnData:
        adata = self._adata
        if self._axis == "obs" and any(value is not None for value in (obs, x, raw, obs_names, obsm)):
            labels: list[str] = []
            for _, group in self._iter_groups(adata):
                filtered = filter_adata(
                    group,
                    obs=obs,
                    x=x,
                    raw=raw,
                    obs_names=obs_names,
                    obsm=obsm,
                    layer=layer,
                    copy=False,
                )
                labels.extend(filtered.obs_names.tolist())
            adata = adata[labels, :]
            obs = x = raw = obs_names = obsm = None
        elif self._axis == "var" and any(value is not None for value in (var, var_names, varm)):
            labels = []
            for _, group in self._iter_groups(adata):
                filtered = filter_adata(group, var=var, var_names=var_names, varm=varm, copy=False)
                labels.extend(filtered.var_names.tolist())
            adata = adata[:, labels]
            var = var_names = varm = None
        return filter_adata(
            adata,
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
    ) -> AnnData:
        if self._adata.isbacked:
            msg = "grouped mutate cannot modify an AnnData object in backed mode; call .to_memory() first"
            raise AnnplyrError(msg)
        if self._axis == "obs" and any(value is not None for value in (obs, x, raw, obsm)):
            out = self._adata if inplace else self._adata.copy()
            by_frame = self._by_frame(out)
            obs_col_pieces: dict[str, list[tuple[np.ndarray, pd.Series]]] = {}
            for _, key_row in by_frame.drop_duplicates().iterrows():
                mask = _key_mask(by_frame, key_row)
                positions = np.where(mask.to_numpy())[0]
                group = out[mask.to_numpy(), :]
                for frame, assignments in _obs_assignment_frames(group, obs=obs, x=x, raw=raw, obsm=obsm, layer=layer):
                    group_values = evaluate_assignments(frame, assignments)
                    for column in group_values.columns:
                        if column not in obs_col_pieces:
                            obs_col_pieces[column] = []
                        obs_col_pieces[column].append((positions, group_values[column].reset_index(drop=True)))
            obs_df = cast(pd.DataFrame, out.obs)
            for column, pieces in obs_col_pieces.items():
                all_pos = np.concatenate([p for p, _ in pieces])
                order = np.argsort(all_pos, kind="stable")
                merged = pd.concat([s for _, s in pieces], ignore_index=True).iloc[order]
                obs_df[column] = merged.to_numpy()
            return mutate_adata(out, var=var, varm=varm, layer=layer, inplace=True)
        if self._axis == "var" and any(value is not None for value in (var, varm)):
            out = self._adata if inplace else self._adata.copy()
            by_frame = self._by_frame(out)
            var_col_pieces: dict[str, list[tuple[np.ndarray, pd.Series]]] = {}
            for _, key_row in by_frame.drop_duplicates().iterrows():
                mask = _key_mask(by_frame, key_row)
                positions = np.where(mask.to_numpy())[0]
                group = out[:, mask.to_numpy()]
                for frame, assignments in _var_assignment_frames(group, var=var, varm=varm):
                    group_values = evaluate_assignments(frame, assignments)
                    for column in group_values.columns:
                        if column not in var_col_pieces:
                            var_col_pieces[column] = []
                        var_col_pieces[column].append((positions, group_values[column].reset_index(drop=True)))
            var_df = cast(pd.DataFrame, out.var)
            for column, pieces in var_col_pieces.items():
                all_pos = np.concatenate([p for p, _ in pieces])
                order = np.argsort(all_pos, kind="stable")
                merged = pd.concat([s for _, s in pieces], ignore_index=True).iloc[order]
                var_df[column] = merged.to_numpy()
            return mutate_adata(out, obs=obs, x=x, raw=raw, obsm=obsm, layer=layer, inplace=True)
        return mutate_adata(
            self._adata, obs=obs, var=var, x=x, raw=raw, obsm=obsm, varm=varm, layer=layer, inplace=inplace
        )

    def summarize(
        self,
        obs: Mapping[str, Any] | None = None,
        var: Mapping[str, Any] | None = None,
        x: Mapping[str, Any] | None = None,
        raw: Mapping[str, Any] | None = None,
        obsm: Mapping[str, Mapping[str, Any]] | None = None,
        varm: Mapping[str, Mapping[str, Any]] | None = None,
        layer: str | None = None,
    ) -> pd.DataFrame:
        return summarize_adata(
            self._adata,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            by=self._by,
            layer=layer,
        )

    summarise = summarize

    def count(self, *, wt: Any = None, sort: bool = False, name: str = "n") -> pd.DataFrame:
        frame = obs_frame(self._adata) if self._axis == "obs" else var_frame(self._adata)
        return count_frame(frame, by=self._by, wt=wt, sort=sort, name=name)

    def tally(self, *, wt: Any = None, sort: bool = False, name: str = "n") -> pd.DataFrame:
        return self.count(wt=wt, sort=sort, name=name)

    def add_count(self, *, wt: Any = None, sort: bool = False, name: str = "n", inplace: bool = False) -> AnnData:
        return add_count_adata(self._adata, by=self._by, wt=wt, sort=sort, axis=self._axis, name=name, inplace=inplace)

    def add_tally(self, *, wt: Any = None, sort: bool = False, name: str = "n", inplace: bool = False) -> AnnData:
        return add_count_adata(self._adata, by=self._by, wt=wt, sort=sort, axis=self._axis, name=name, inplace=inplace)

    def arrange(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        raw: Any = None,
        obsm: Mapping[str, Any] | None = None,
        varm: Mapping[str, Any] | None = None,
        layer: str | None = None,
        copy: bool = False,
    ) -> AnnData:
        adata = self._adata
        if self._axis == "obs" and any(value is not None for value in (obs, x, raw, obsm)):
            labels: list[str] = []
            for _, group in self._iter_groups(adata):
                arranged = arrange_adata(group, obs=obs, x=x, raw=raw, obsm=obsm, layer=layer, copy=False)
                labels.extend(arranged.obs_names.tolist())
            adata = adata[labels, :]
            obs = x = raw = obsm = None
        elif self._axis == "var" and any(value is not None for value in (var, varm)):
            labels = []
            for _, group in self._iter_groups(adata):
                arranged = arrange_adata(group, var=var, varm=varm, copy=False)
                labels.extend(arranged.var_names.tolist())
            adata = adata[:, labels]
            var = varm = None
        return arrange_adata(adata, obs=obs, var=var, x=x, raw=raw, obsm=obsm, varm=varm, layer=layer, copy=copy)

    def distinct(
        self,
        obs: Any = None,
        var: Any = None,
        x: Any = None,
        *,
        keep_all: bool = False,
        copy: bool = True,
    ) -> AnnData:
        if self._axis == "obs":
            if var is not None:
                msg = "obs-grouped distinct cannot use var selectors"
                raise IncompatibleAxisError(msg)
            labels: list[str] = []
            for _, group in self._iter_groups(self._adata):
                distinct = distinct_adata(group, obs=obs, x=x, axis="obs", keep_all=keep_all, copy=False)
                labels.extend(distinct.obs_names.tolist())
            return self._subset_labels(labels, copy=copy)
        if obs is not None or x is not None:
            msg = "var-grouped distinct cannot use obs or x selectors"
            raise IncompatibleAxisError(msg)
        labels = []
        for _, group in self._iter_groups(self._adata):
            distinct = distinct_adata(group, var=var, axis="var", keep_all=keep_all, copy=False)
            labels.extend(distinct.var_names.tolist())
        return self._subset_labels(labels, copy=copy)

    def slice_head(self, n: int = 5, *, copy: bool = False) -> AnnData:
        return self._slice_group_positions(n=n, tail=False, copy=copy)

    def slice_tail(self, n: int = 5, *, copy: bool = False) -> AnnData:
        return self._slice_group_positions(n=n, tail=True, copy=copy)

    def slice_min(self, by: Any, n: int = 5, *, copy: bool = False) -> AnnData:
        labels: list[str] = []
        for _, group in self._iter_groups(self._adata):
            sliced = slice_min_adata(group, by=by, n=n, axis=self._axis)
            labels.extend(_axis_names(sliced, self._axis))
        return self._subset_labels(labels, copy=copy)

    def slice_max(self, by: Any, n: int = 5, *, copy: bool = False) -> AnnData:
        labels: list[str] = []
        for _, group in self._iter_groups(self._adata):
            sliced = slice_max_adata(group, by=by, n=n, axis=self._axis)
            labels.extend(_axis_names(sliced, self._axis))
        return self._subset_labels(labels, copy=copy)

    def slice_sample(
        self,
        n: int | None = None,
        *,
        prop: float | None = None,
        replace: bool = False,
        random_state: int | None = None,
        copy: bool = False,
    ) -> AnnData:
        if n is not None and prop is not None:
            msg = "slice_sample accepts n and prop as mutually exclusive arguments"
            raise SelectionError(msg)
        if n is not None and n < 0:
            msg = "slice_sample n must be non-negative"
            raise SelectionError(msg)
        if prop is not None and prop < 0:
            msg = "slice_sample prop must be non-negative"
            raise SelectionError(msg)
        rng = np.random.default_rng(random_state)
        labels: list[str] = []
        for _, group in self._iter_groups(self._adata):
            axis_names = _axis_names(group, self._axis)
            size = len(axis_names)
            take = int(round(size * prop)) if n is None and prop is not None else (n if n is not None else min(size, 1))
            if not replace and take > size:
                msg = "slice_sample n cannot be larger than a group size unless replace=True"
                raise SelectionError(msg)
            selected = rng.choice(axis_names, size=take, replace=replace)
            labels.extend([str(label) for label in selected])
        return self._subset_labels(labels, copy=copy)

    def _slice_group_positions(self, *, n: int, tail: bool, copy: bool) -> AnnData:
        _validate_slice_n(n)
        labels: list[str] = []
        for _, group in self._iter_groups(self._adata):
            axis_names = group.obs_names if self._axis == "obs" else group.var_names
            selected = axis_names[0:0] if n == 0 else (axis_names[-n:] if tail else axis_names[:n])
            labels.extend(selected.tolist())
        if self._axis == "obs":
            out = self._adata[labels, :]
        else:
            out = self._adata[:, labels]
        return out.copy() if copy else out

    def _subset_labels(self, labels: list[str], *, copy: bool) -> AnnData:
        if self._axis == "obs":
            out = self._adata[labels, :]
        else:
            out = self._adata[:, labels]
        return out.copy() if copy else out


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


def _axis_names(adata: AnnData, axis: str) -> list[str]:
    return adata.obs_names.tolist() if axis == "obs" else adata.var_names.tolist()
