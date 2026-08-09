from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import narwhals as nw
import numpy as np
import pandas as pd

from annplyr._errors import DuplicateNameError, SelectionError
from annplyr._expr import _ReductionSpec, reduction_spec
from annplyr._frames import evaluate_select, expand_assignments, with_row_number
from annplyr._groups import GroupPlan
from annplyr._sources import ChunkPlan

_MISSING_UNIQUE = object()


@dataclass(frozen=True)
class SummaryGroupPlan:
    """Evaluated grouping keys and positional ids shared by summary sources."""

    group_ids: np.ndarray
    keys: pd.DataFrame
    columns: tuple[str, ...]
    internal_columns: tuple[str, ...]

    @classmethod
    def build(
        cls,
        by_source: pd.DataFrame,
        by: Any,
        *,
        grouped: GroupPlan | None = None,
    ) -> SummaryGroupPlan:
        if grouped is not None:
            columns = tuple(str(column) for column in grouped.keys.columns)
            internal = tuple(f"__annplyr_by_{position}__" for position in range(len(columns)))
            keys = grouped.keys.copy()
            keys.columns = list(internal)
            return cls(grouped.group_ids, keys, columns, internal)
        if by is None:
            return cls(
                np.zeros(len(by_source), dtype=np.intp),
                pd.DataFrame(index=pd.RangeIndex(1)),
                (),
                (),
            )

        by_values = evaluate_select(by_source, by).reset_index(drop=True)
        columns = tuple(str(column) for column in by_values.columns)
        internal = tuple(f"__annplyr_by_{position}__" for position in range(len(columns)))
        if by_values.empty:
            keys = by_values.copy()
            keys.columns = list(internal)
            return cls(np.empty(0, dtype=np.intp), keys, columns, internal)

        grouped_values = by_values.groupby(list(by_values.columns), sort=False, observed=True, dropna=False)
        raw_ids = grouped_values.ngroup().to_numpy(dtype=np.intp)
        raw_group_count = int(raw_ids.max()) + 1
        row_positions = np.arange(len(raw_ids), dtype=np.intp)
        first_by_raw_group = np.full(raw_group_count, len(raw_ids), dtype=np.intp)
        np.minimum.at(first_by_raw_group, raw_ids, row_positions)
        raw_group_order = np.argsort(first_by_raw_group, kind="stable")
        dense_by_raw_group = np.empty(raw_group_count, dtype=np.intp)
        dense_by_raw_group[raw_group_order] = np.arange(raw_group_count, dtype=np.intp)
        group_ids = dense_by_raw_group[raw_ids]
        keys = by_values.iloc[first_by_raw_group[raw_group_order], :].reset_index(drop=True).copy()
        keys.columns = list(internal)
        return cls(group_ids, keys, columns, internal)

    @property
    def count(self) -> int:
        return len(self.keys)

    def add_to_frame(self, frame: pd.DataFrame, row_positions: np.ndarray) -> pd.DataFrame:
        work = frame.copy(deep=False)
        for position, internal in enumerate(self.internal_columns):
            values = self.keys.iloc[self.group_ids[row_positions], position]
            work[internal] = pd.Series(values.array, index=work.index)
        return work

    def restore_keys(self, result: pd.DataFrame) -> pd.DataFrame:
        renamed = result.rename(columns=dict(zip(self.internal_columns, self.columns, strict=True)))
        if len(renamed) != self.count:
            return renamed
        for source, target in zip(self.internal_columns, self.columns, strict=True):
            if target in renamed.columns:
                renamed[target] = pd.Series(self.keys[source].array, index=renamed.index)
        return renamed


@dataclass(frozen=True)
class ReductionPlan:
    """Canonical scalar assignments eligible for chunked execution."""

    assignments: tuple[tuple[str, _ReductionSpec], ...]

    @classmethod
    def resolve(cls, schema: pd.DataFrame, assignments: Any) -> ReductionPlan | None:
        expanded = expand_assignments(schema, assignments)
        planned: list[tuple[str, _ReductionSpec]] = []
        for name, expression in expanded.items():
            spec = reduction_spec(expression)
            if spec is None:
                return None
            planned.append((str(name), spec))
        return cls(tuple(planned))


class _PairwiseSumState:
    """Stream NumPy's pairwise summation tree without retaining source values."""

    def __init__(self, size: int, dtype: np.dtype[Any]) -> None:
        self.size = size
        self.dtype = dtype
        self._leaf_sizes = iter(_pairwise_leaf_sizes(size))
        self._target = next(self._leaf_sizes, 0)
        self._parts: list[np.ndarray] = []
        self._filled = 0
        self._leaves: list[Any] = []

    def update(self, values: np.ndarray) -> None:
        array = np.asarray(values, dtype=self.dtype)
        if self.dtype.kind in "fc" and np.isnan(array).any():
            array = array.copy()
            array[np.isnan(array)] = 0
        offset = 0
        while offset < len(array):
            take = min(self._target - self._filled, len(array) - offset)
            self._parts.append(array[offset : offset + take])
            self._filled += take
            offset += take
            if self._filled == self._target:
                leaf = self._parts[0] if len(self._parts) == 1 else np.concatenate(self._parts)
                self._leaves.append(np.sum(leaf, dtype=self.dtype))
                self._parts = []
                self._filled = 0
                self._target = next(self._leaf_sizes, 0)

    def finish(self) -> Any:
        if self._target or self._parts:
            raise RuntimeError("pairwise reduction received the wrong number of rows")
        if self.size == 0:
            return self.dtype.type(0)
        leaves = iter(self._leaves)
        value = _combine_pairwise_leaves(self.size, leaves, self.dtype)
        try:
            next(leaves)
        except StopIteration:
            return value
        raise RuntimeError("pairwise reduction retained unexpected leaf values")


class _Accumulator:
    def __init__(
        self,
        spec: _ReductionSpec,
        groups: int,
        *,
        full_pass: bool = False,
        group_sizes: np.ndarray | None = None,
    ) -> None:
        self.spec = spec
        self.groups = groups
        self.full_pass = full_pass
        self.dtype: Any = None
        self.counts = np.zeros(groups, dtype=np.int64)
        self.totals: list[Any] = [None] * groups
        self.values: list[Any] = [None] * groups
        self.seen = np.zeros(groups, dtype=bool)
        self.means = np.zeros(groups, dtype=np.longdouble)
        self.m2 = np.zeros(groups, dtype=np.longdouble)
        self.collections: list[Any] = [None] * groups
        self.group_sizes = group_sizes
        self.pairwise: list[_PairwiseSumState] | None = None

    def observe_dtype(self, series: pd.Series) -> None:
        if self.dtype is None:
            self.dtype = series.dtype
            dtype = series.dtype.subtype if isinstance(series.dtype, pd.SparseDtype) else series.dtype
            try:
                numpy_dtype = np.dtype(cast(Any, dtype))
            except TypeError:
                return
            if (
                not self.full_pass
                and self.spec.operation in {"mean", "sum"}
                and numpy_dtype.kind in "fc"
                and self.group_sizes is not None
            ):
                self.pairwise = [_PairwiseSumState(int(size), numpy_dtype) for size in self.group_sizes]

    def update(self, series: pd.Series | None, group_ids: np.ndarray) -> None:
        operation = self.spec.operation
        if operation == "len":
            self.counts += np.bincount(group_ids, minlength=self.groups)
            return
        if series is None:
            raise RuntimeError(f"reduction {operation!r} requires an input series")
        self.observe_dtype(series)
        for group in np.unique(group_ids):
            local = series.iloc[np.flatnonzero(group_ids == group)]
            self._update_group(int(group), local)

    def _update_group(self, group: int, local: pd.Series) -> None:
        operation = self.spec.operation
        if operation == "sum":
            if self.pairwise is not None:
                self.pairwise[group].update(local.to_numpy())
            else:
                value = _local_sum(local, preserve_eager=self.full_pass)
                self.totals[group] = value if self.totals[group] is None else self.totals[group] + value
        elif operation == "mean":
            count = int(local.count())
            if self.pairwise is not None:
                self.pairwise[group].update(local.to_numpy())
                self.counts[group] += count
            elif count:
                value = _local_sum(local, preserve_eager=self.full_pass)
                self.totals[group] = value if self.totals[group] is None else self.totals[group] + value
                self.counts[group] += count
        elif operation == "std":
            if self.full_pass and not isinstance(local.dtype, pd.SparseDtype):
                self.values[group] = local.std(skipna=True, ddof=self.spec.ddof)
                self.seen[group] = True
                return
            dense = np.asarray(local.dropna().to_numpy(), dtype=np.longdouble)
            count = len(dense)
            if count:
                mean = dense.mean(dtype=np.longdouble)
                m2 = np.square(dense - mean).sum(dtype=np.longdouble)
                previous = int(self.counts[group])
                total = previous + count
                delta = mean - self.means[group]
                self.m2[group] += m2 + delta * delta * previous * count / total
                self.means[group] += delta * count / total
                self.counts[group] = total
        elif operation in {"min", "max"}:
            value = getattr(local, operation)(skipna=True)
            if not _is_missing(value):
                if not self.seen[group]:
                    self.values[group] = value
                    self.seen[group] = True
                else:
                    combine = np.minimum if operation == "min" else np.maximum
                    self.values[group] = combine(self.values[group], value)
        elif operation == "first":
            if not self.seen[group] and len(local):
                self.values[group] = local.iloc[0]
                self.seen[group] = True
        elif operation == "last":
            if len(local):
                self.values[group] = local.iloc[-1]
                self.seen[group] = True
        elif operation == "median":
            dense = np.asarray(local.dropna().to_numpy())
            if len(dense):
                if self.collections[group] is None:
                    self.collections[group] = []
                self.collections[group].append(dense.copy())
        elif operation == "n_unique":
            if self.collections[group] is None:
                self.collections[group] = set()
            unique = cast(set[Any], self.collections[group])
            for value in pd.unique(local):
                unique.add(_unique_key(value))
        elif operation == "all":
            value = bool(local.all(skipna=True))
            self.values[group] = value if not self.seen[group] else bool(self.values[group]) and value
            self.seen[group] = True
        elif operation == "any":
            value = bool(local.any(skipna=True))
            self.values[group] = value if not self.seen[group] else bool(self.values[group]) or value
            self.seen[group] = True
        else:  # pragma: no cover - guarded by expression metadata
            raise RuntimeError(f"unknown chunked reduction: {operation!r}")

    def finish(self) -> pd.Series:
        operation = self.spec.operation
        values: list[Any] = []
        for group in range(self.groups):
            if operation == "len":
                value = self.counts[group]
            elif operation == "sum":
                value = self.pairwise[group].finish() if self.pairwise is not None else self.totals[group]
                value = 0 if value is None else value
            elif operation == "mean":
                total = self.pairwise[group].finish() if self.pairwise is not None else self.totals[group]
                value = np.nan if not self.counts[group] else total / self.counts[group]
            elif operation == "std":
                if self.seen[group]:
                    value = self.values[group]
                else:
                    denominator = int(self.counts[group]) - self.spec.ddof
                    value = np.nan if denominator <= 0 else np.sqrt(self.m2[group] / denominator)
            elif operation in {"min", "max", "first", "last"}:
                value = self.values[group] if self.seen[group] else np.nan
            elif operation == "median":
                chunks = self.collections[group]
                value = np.nan if not chunks else np.median(np.concatenate(chunks))
            elif operation == "n_unique":
                value = 0 if self.collections[group] is None else len(self.collections[group])
            elif operation == "all":
                value = self.values[group] if self.seen[group] else True
            elif operation == "any":
                value = self.values[group] if self.seen[group] else False
            else:  # pragma: no cover - guarded by expression metadata
                raise RuntimeError(f"unknown chunked reduction: {operation!r}")
            values.append(value)
        return _typed_result(values, self.dtype, operation)


def summarize_chunked(
    chunks: ChunkPlan,
    reductions: ReductionPlan,
    groups: SummaryGroupPlan,
    *,
    virtual_name: str | None = None,
    virtual_values: np.ndarray | None = None,
) -> pd.DataFrame:
    """Evaluate canonical scalar reductions with bounded source reads."""
    if _can_chunk_columns(chunks, reductions, virtual_name=virtual_name):
        return _summarize_column_chunks(chunks, reductions, groups)
    group_sizes = np.bincount(groups.group_ids, minlength=groups.count)
    accumulators = [_Accumulator(spec, groups.count, group_sizes=group_sizes) for _, spec in reductions.assignments]
    schema = chunks.request.adapter.schema.iloc[:0, :].copy()
    if virtual_name is not None and virtual_values is not None:
        schema[virtual_name] = pd.Series([], dtype=virtual_values.dtype)
    schema_inputs = _evaluate_inputs(schema, reductions)
    for accumulator, series in zip(accumulators, schema_inputs, strict=True):
        if series is not None:
            accumulator.observe_dtype(series)

    for row_positions, frame in chunks.read_chunks():
        if virtual_name is not None and virtual_values is not None:
            frame = frame.copy(deep=False)
            frame[virtual_name] = virtual_values[row_positions]
        inputs = _evaluate_inputs(frame, reductions)
        group_ids = groups.group_ids[row_positions]
        for accumulator, series in zip(accumulators, inputs, strict=True):
            accumulator.update(series, group_ids)

    result = groups.keys.copy()
    for (name, _), accumulator in zip(reductions.assignments, accumulators, strict=True):
        result[name] = accumulator.finish().array
    return result


def _can_chunk_columns(
    chunks: ChunkPlan,
    reductions: ReductionPlan,
    *,
    virtual_name: str | None,
) -> bool:
    if len(chunks.request.row_positions) > chunks.target_values:
        return False
    available = [chunks.request.adapter.names[position] for position in chunks.request.column_positions]
    for _, spec in reductions.assignments:
        if spec.input_expr is None:
            continue
        if spec.source_column is None or spec.source_column == virtual_name:
            return False
        matches = available.count(spec.source_column)
        if matches > 1:
            raise DuplicateNameError(f"Reduction source column {spec.source_column!r} is duplicated")
        if matches == 0:
            return False
    return True


def _summarize_column_chunks(
    chunks: ChunkPlan,
    reductions: ReductionPlan,
    groups: SummaryGroupPlan,
) -> pd.DataFrame:
    if not any(isinstance(dtype, pd.SparseDtype) for dtype in chunks.request.adapter.schema.dtypes):
        return _summarize_dense_column_chunks(chunks, reductions, groups)
    accumulators: list[_Accumulator | None] = [
        _Accumulator(spec, groups.count, full_pass=True) for _, spec in reductions.assignments
    ]
    finished: list[pd.Series | None] = [None] * len(accumulators)
    names = chunks.request.adapter.names
    name_counts: dict[str, int] = {}
    for name in names:
        name_counts[name] = name_counts.get(name, 0) + 1
    source_positions = {name: position for position, name in enumerate(names) if name_counts[name] == 1}
    reduction_positions: dict[int, list[int]] = {}
    length_positions: list[int] = []
    for reduction_position, (_, spec) in enumerate(reductions.assignments):
        if spec.input_expr is None:
            length_positions.append(reduction_position)
        elif spec.source_column is not None:
            reduction_positions.setdefault(source_positions[spec.source_column], []).append(reduction_position)

    rows = chunks.request.row_positions
    group_ids = groups.group_ids[rows]
    for position in length_positions:
        accumulator = accumulators[position]
        if accumulator is None:  # pragma: no cover - unique assignment positions
            raise RuntimeError("length reduction was already finalized")
        accumulator.update(None, group_ids)
        finished[position] = accumulator.finish()
        accumulators[position] = None

    for columns, frame in chunks.read_column_chunks():
        local_by_source = {int(source): offset for offset, source in enumerate(columns)}
        batch_positions = [
            position
            for source, positions in reduction_positions.items()
            if source in local_by_source
            for position in positions
        ]
        for source, positions in reduction_positions.items():
            if source not in local_by_source:
                continue
            series = frame.iloc[:, local_by_source[source]].reset_index(drop=True)
            for position in positions:
                accumulator = accumulators[position]
                if accumulator is None:  # pragma: no cover - each source is read once
                    raise RuntimeError("source reduction was already finalized")
                accumulator.update(series, group_ids)
        for position in batch_positions:
            accumulator = accumulators[position]
            if accumulator is None:  # pragma: no cover - each source is read once
                raise RuntimeError("source reduction was already finalized")
            finished[position] = accumulator.finish()
            accumulators[position] = None

    result = groups.keys.copy()
    for (name, _), values in zip(reductions.assignments, finished, strict=True):
        if values is None:  # pragma: no cover - guarded by reduction/source planning
            raise RuntimeError(f"reduction {name!r} was not evaluated")
        result[name] = values.array
    return result


def _summarize_dense_column_chunks(
    chunks: ChunkPlan,
    reductions: ReductionPlan,
    groups: SummaryGroupPlan,
) -> pd.DataFrame:
    names = chunks.request.adapter.names
    source_positions = {name: position for position, name in enumerate(names)}
    reductions_by_source: dict[int, list[tuple[str, _ReductionSpec]]] = {}
    length_reductions: list[tuple[str, _ReductionSpec]] = []
    for name, spec in reductions.assignments:
        if spec.input_expr is None:
            length_reductions.append((name, spec))
        elif spec.source_column is not None:
            reductions_by_source.setdefault(source_positions[spec.source_column], []).append((name, spec))

    pieces: list[pd.DataFrame] = []
    include_lengths = True
    rows = chunks.request.row_positions
    for columns, frame in chunks.read_column_chunks():
        batch = list(length_reductions) if include_lengths else []
        include_lengths = False
        for source in columns:
            batch.extend(reductions_by_source.get(int(source), []))
        if not batch:
            continue
        work = groups.add_to_frame(frame, rows)
        expressions = [_reduction_expression(name, spec) for name, spec in batch]
        if groups.internal_columns:
            piece = nw.from_native(work).group_by(*groups.internal_columns).agg(*expressions).to_native()
        else:
            piece = nw.from_native(with_row_number(work)).select(*expressions).to_native()
        pieces.append(piece)

    if not pieces:
        raise RuntimeError("dense column reduction plan produced no batches")
    result = pieces[0].reset_index(drop=True)
    for piece in pieces[1:]:
        piece = piece.reset_index(drop=True)
        if groups.internal_columns:
            result = result.merge(piece, on=list(groups.internal_columns), how="outer", sort=False)
        else:
            result = pd.concat([result, piece], axis=1)
    if len(result) == groups.count:
        for internal in groups.internal_columns:
            result[internal] = pd.Series(groups.keys[internal].array, index=result.index)
    ordered = [*groups.internal_columns, *(name for name, _ in reductions.assignments)]
    return result.loc[:, ordered]


def _reduction_expression(name: str, spec: _ReductionSpec) -> nw.Expr:
    if spec.input_expr is None:
        return nw.len().alias(name)
    method = getattr(spec.input_expr, spec.operation)
    expression = method(ddof=spec.ddof) if spec.operation == "std" else method()
    return expression.alias(name)


def _evaluate_inputs(
    frame: pd.DataFrame,
    reductions: ReductionPlan,
) -> list[pd.Series | None]:
    inputs: list[pd.Series | None] = [None] * len(reductions.assignments)
    derived: list[tuple[int, str, nw.Expr]] = []
    direct: dict[str, pd.Series] = {}
    for position, (_, spec) in enumerate(reductions.assignments):
        if spec.input_expr is None:
            continue
        if spec.source_column is not None:
            selected = direct.get(spec.source_column)
            if selected is not None:
                inputs[position] = selected
                continue
            selected = frame.loc[:, spec.source_column]
            if isinstance(selected, pd.DataFrame):
                raise DuplicateNameError(f"Reduction source column {spec.source_column!r} is duplicated")
            selected = selected.reset_index(drop=True)
            direct[spec.source_column] = selected
            inputs[position] = selected
            continue
        name = f"__annplyr_reduce_{position}__"
        derived.append((position, name, spec.input_expr.alias(name)))
    if derived:
        try:
            evaluated = nw.from_native(with_row_number(frame)).select(*(expr for _, _, expr in derived)).to_native()
        except Exception as exc:
            raise SelectionError(f"Chunked reduction input evaluation failed: {exc}") from exc
        for position, name, _ in derived:
            inputs[position] = evaluated[name].reset_index(drop=True)
    return inputs


def _typed_result(values: list[Any], dtype: Any, operation: str) -> pd.Series:
    if operation in {"len", "n_unique"}:
        return pd.Series(values, dtype="int64")
    if operation in {"all", "any"}:
        return pd.Series(values, dtype="bool")
    dense_dtype = dtype.subtype if isinstance(dtype, pd.SparseDtype) else dtype
    try:
        numpy_dtype = np.dtype(dense_dtype)
    except TypeError:
        return pd.Series(values, dtype=dense_dtype)
    if operation in {"mean", "median", "std"} and numpy_dtype.kind not in "fc":
        numpy_dtype = np.dtype("float64")
    elif operation == "sum":
        if numpy_dtype.kind in "bi" and numpy_dtype.itemsize < 8:
            numpy_dtype = np.dtype("int64")
        elif numpy_dtype.kind == "u" and numpy_dtype.itemsize < 8:
            numpy_dtype = np.dtype("uint64")
    if numpy_dtype.kind in "biu" and any(_is_missing(value) for value in values):
        numpy_dtype = np.dtype("float64")
    return pd.Series(np.asarray(values, dtype=numpy_dtype))


def _local_sum(series: pd.Series, *, preserve_eager: bool = False) -> Any:
    if preserve_eager:
        return series.sum(skipna=True)
    dtype = series.dtype.subtype if isinstance(series.dtype, pd.SparseDtype) else series.dtype
    try:
        numpy_dtype = np.dtype(cast(Any, dtype))
    except TypeError:
        return series.sum(skipna=True)
    if numpy_dtype.kind == "f":
        accumulator_dtype = np.float64 if numpy_dtype.itemsize < 8 else np.longdouble
        return np.nansum(series.to_numpy(), dtype=accumulator_dtype)
    if numpy_dtype.kind == "c":
        return np.nansum(series.to_numpy(), dtype=np.clongdouble)
    return series.sum(skipna=True)


def _pairwise_leaf_sizes(size: int) -> list[int]:
    if size <= 128:
        return [size] if size else []
    left = size // 2
    left -= left % 8
    return [*_pairwise_leaf_sizes(left), *_pairwise_leaf_sizes(size - left)]


def _combine_pairwise_leaves(size: int, leaves: Any, dtype: np.dtype[Any]) -> Any:
    if size <= 128:
        return next(leaves)
    left = size // 2
    left -= left % 8
    left_value = _combine_pairwise_leaves(left, leaves, dtype)
    right_value = _combine_pairwise_leaves(size - left, leaves, dtype)
    return dtype.type(left_value + right_value)


def _is_missing(value: Any) -> bool:
    missing = pd.isna(value)
    return bool(missing) if isinstance(missing, (bool, np.bool_)) else False


def _unique_key(value: Any) -> Any:
    if _is_missing(value):
        return _MISSING_UNIQUE
    return value.item() if isinstance(value, np.generic) else value
