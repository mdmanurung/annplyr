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

    @classmethod
    def build(cls, adata: AnnData, spec: GroupSpec) -> GroupPlan:
        spec.validate(adata)
        table = cast(pd.DataFrame, adata.obs if spec.axis == "obs" else adata.var)
        key_frame = table.loc[:, list(spec.columns)].reset_index(drop=True)
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
        buckets: dict[int, list[int]] = {}
        for position, raw_group_id in enumerate(raw_ids):
            buckets.setdefault(int(raw_group_id), []).append(position)
        ordered = sorted(buckets.values(), key=lambda positions: positions[0])
        position_arrays = tuple(np.asarray(positions, dtype=np.intp) for positions in ordered)
        first_positions = np.asarray([positions[0] for positions in ordered], dtype=np.intp)
        group_ids = np.empty(len(raw_ids), dtype=np.intp)
        for dense_group_id, positions in enumerate(position_arrays):
            group_ids[positions] = dense_group_id
        keys = key_frame.iloc[first_positions, :].reset_index(drop=True).copy()
        return cls(spec, group_ids, first_positions, position_arrays, keys)

    def group_data(self) -> pd.DataFrame:
        result = self.keys.copy()
        result[".rows"] = pd.Series([positions.tolist() for positions in self.positions], dtype=object)
        return result
