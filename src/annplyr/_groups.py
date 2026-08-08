from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, cast

import narwhals as nw
import numpy as np
import pandas as pd
from anndata import AnnData

from annplyr._errors import IncompatibleAxisError, SelectionError, UnknownColumnError
from annplyr._expr import AnnplyrExpr
from annplyr._frames import VIRTUAL_COLUMNS, evaluate_select, obs_frame, var_frame

Axis = Literal["obs", "var"]
_MISSING_OBJECT_KEY = object()


@dataclass(frozen=True)
class GroupSpec:
    """Concrete grouping axis and existing metadata column names."""

    axis: Axis
    columns: tuple[str, ...]

    @classmethod
    def resolve(cls, adata: AnnData, *, obs: Any = None, var: Any = None) -> GroupSpec:
        if obs is not None and var is not None:
            raise IncompatibleAxisError("group_by accepts one axis at a time")
        if obs is None and var is None:
            raise SelectionError("group_by requires obs or var metadata columns")
        axis: Axis = "obs" if obs is not None else "var"
        selector = obs if obs is not None else var
        if isinstance(selector, (AnnplyrExpr, nw.Expr)):
            raise SelectionError("group_by accepts existing metadata names or tidy selectors, not computed expressions")
        frame = obs_frame(adata) if axis == "obs" else var_frame(adata)
        selected = evaluate_select(frame, selector)
        columns = tuple(str(column) for column in selected.columns)
        available = cast(pd.DataFrame, adata.obs if axis == "obs" else adata.var).columns
        if not columns:
            raise SelectionError("group_by selection must contain at least one metadata column")
        if any(column in VIRTUAL_COLUMNS or column not in available for column in columns):
            raise SelectionError("group_by accepts only existing obs or var metadata columns")
        if len(set(columns)) != len(columns):
            raise SelectionError("group_by columns must be unique")
        return cls(axis, columns)

    def validate(self, adata: AnnData) -> None:
        table = cast(pd.DataFrame, adata.obs if self.axis == "obs" else adata.var)
        missing = [column for column in self.columns if column not in table.columns]
        if missing:
            raise UnknownColumnError(f"Stored grouping column(s) are missing: {', '.join(missing)}")

    def renamed(self, mapping: dict[str, str]) -> GroupSpec:
        columns = tuple(mapping.get(column, column) for column in self.columns)
        if len(set(columns)) != len(columns):
            raise SelectionError("grouped rename would make grouping columns ambiguous")
        return GroupSpec(self.axis, columns)


@dataclass(frozen=True)
class GroupPlan:
    """One O(n + g) positional plan for a grouped public call."""

    spec: GroupSpec
    group_ids: np.ndarray
    first_positions: np.ndarray
    positions: tuple[np.ndarray, ...]
    keys: pd.DataFrame

    @staticmethod
    def _key_frame(adata: AnnData, spec: GroupSpec) -> pd.DataFrame:
        spec.validate(adata)
        table = cast(pd.DataFrame, adata.obs if spec.axis == "obs" else adata.var)
        return table.loc[:, list(spec.columns)].reset_index(drop=True)

    @staticmethod
    def _first_key_positions(key_frame: pd.DataFrame) -> np.ndarray:
        """Find first-seen keys while treating all object null sentinels alike."""
        comparison = key_frame
        for column in key_frame.columns:
            series = key_frame[column]
            if not pd.api.types.is_object_dtype(series.dtype):
                continue
            missing = series.isna().to_numpy()
            if not missing.any():
                continue
            if comparison is key_frame:
                comparison = key_frame.copy()
            normalized = series.to_numpy(dtype=object, copy=True)
            normalized[missing] = _MISSING_OBJECT_KEY
            comparison[column] = normalized
        return np.flatnonzero(~comparison.duplicated(keep="first").to_numpy())

    @classmethod
    def keys_for(cls, adata: AnnData, spec: GroupSpec) -> pd.DataFrame:
        """Return observed keys in first-seen order without building row groups."""
        key_frame = cls._key_frame(adata, spec)
        positions = cls._first_key_positions(key_frame)
        return key_frame.iloc[positions, :].reset_index(drop=True).copy()

    @classmethod
    def count_for(
        cls,
        adata: AnnData,
        spec: GroupSpec,
        *,
        weights: pd.Series | None,
        name: str,
    ) -> pd.DataFrame:
        """Aggregate observed groups without constructing positional arrays."""
        key_frame = cls._key_frame(adata, spec)
        positions = cls._first_key_positions(key_frame)
        keys = key_frame.iloc[positions, :].reset_index(drop=True).copy()
        if key_frame.empty:
            keys[name] = pd.Series(dtype="float64")
            return keys
        work = key_frame.copy()
        weight_column = "__annplyr_group_weight__"
        while weight_column in work.columns:
            weight_column += "_"
        if weights is not None:
            work[weight_column] = weights.array
        grouped = work.groupby(list(spec.columns), sort=False, observed=True, dropna=False)
        values = grouped.size() if weights is None else grouped[weight_column].sum()
        if len(values) != len(keys):
            raise RuntimeError("group count keys do not match observed groups")
        keys[name] = pd.Series(values.array, index=keys.index)
        return keys

    @classmethod
    def build(cls, adata: AnnData, spec: GroupSpec) -> GroupPlan:
        key_frame = cls._key_frame(adata, spec)
        if key_frame.empty:
            return cls(
                spec,
                np.empty(0, dtype=np.intp),
                np.empty(0, dtype=np.intp),
                (),
                key_frame.iloc[:0, :].copy(),
            )

        grouped = key_frame.groupby(
            list(spec.columns),
            sort=False,
            observed=True,
            dropna=False,
        )
        raw_ids = grouped.ngroup().to_numpy(dtype=np.intp)
        raw_group_count = int(raw_ids.max()) + 1
        row_positions = np.arange(len(raw_ids), dtype=np.intp)
        first_by_raw_group = np.full(raw_group_count, len(raw_ids), dtype=np.intp)
        np.minimum.at(first_by_raw_group, raw_ids, row_positions)
        raw_group_order = np.argsort(first_by_raw_group, kind="stable")
        dense_by_raw_group = np.empty(raw_group_count, dtype=np.intp)
        dense_by_raw_group[raw_group_order] = np.arange(raw_group_count, dtype=np.intp)
        group_ids = dense_by_raw_group[raw_ids]
        counts = np.bincount(group_ids, minlength=raw_group_count)
        grouped_positions = np.argsort(group_ids, kind="stable")
        position_arrays = tuple(np.split(grouped_positions, np.cumsum(counts)[:-1]))
        first_positions = first_by_raw_group[raw_group_order]
        keys = key_frame.iloc[first_positions, :].reset_index(drop=True).copy()
        return cls(spec, group_ids, first_positions, position_arrays, keys)

    def group_data(self) -> pd.DataFrame:
        result = self.keys.copy()
        result[".rows"] = pd.Series([positions.tolist() for positions in self.positions], dtype=object)
        return result
