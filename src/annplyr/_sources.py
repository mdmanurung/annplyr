from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Protocol, cast

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from annplyr._errors import AnnplyrError, SelectionError, UnknownColumnError, UnknownSourceError
from annplyr._expr import AnnplyrExpr, Desc, expression_dependencies

DEFAULT_REDUCTION_CHUNK_VALUES = 25_165_824


def _positions(indexer: Any, size: int, *, dimension: str) -> np.ndarray:
    base = np.arange(size, dtype=np.intp)
    if indexer is None:
        return base
    if isinstance(indexer, slice):
        return base[indexer]
    values = np.asarray(indexer)
    if values.ndim != 1:
        raise SelectionError(f"{dimension} projection must be one-dimensional")
    if np.issubdtype(values.dtype, np.bool_):
        if len(values) != size:
            raise SelectionError(f"{dimension} boolean projection has the wrong length")
        return base[values]
    if not np.issubdtype(values.dtype, np.integer):
        raise SelectionError(f"{dimension} projection must use integer positions")
    values = values.astype(np.intp, copy=True)
    values[values < 0] += size
    if ((values < 0) | (values >= size)).any():
        raise SelectionError(f"{dimension} projection is out of bounds")
    return values


def _schema_frame(names: Sequence[str], dtypes: Sequence[Any]) -> pd.DataFrame:
    series = [pd.Series([], dtype=dtype) for dtype in dtypes]
    frame = pd.concat(series, axis=1) if series else pd.DataFrame()
    frame.columns = list(names)
    return frame


class SourceAdapter(Protocol):
    """Schema-first boundary for an aligned tabular or matrix source."""

    @property
    def shape(self) -> tuple[int, int]: ...

    @property
    def names(self) -> tuple[str, ...]: ...

    @property
    def dtypes(self) -> tuple[Any, ...]: ...

    @property
    def schema(self) -> pd.DataFrame: ...

    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame: ...


@dataclass(frozen=True)
class PandasAdapter:
    value: pd.DataFrame

    @property
    def shape(self) -> tuple[int, int]:
        return self.value.shape

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(str(name) for name in self.value.columns)

    @property
    def dtypes(self) -> tuple[Any, ...]:
        return tuple(self.value.dtypes)

    @property
    def schema(self) -> pd.DataFrame:
        return self.value.iloc[:0, :].reset_index(drop=True).copy()

    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        return self.value.iloc[row_positions, column_positions].reset_index(drop=True).copy()


@dataclass(frozen=True)
class DenseAdapter:
    value: Any
    _names: tuple[str, ...]

    @property
    def shape(self) -> tuple[int, int]:
        shape = self.value.shape
        return int(shape[0]), int(shape[1] if len(shape) > 1 else 1)

    @property
    def names(self) -> tuple[str, ...]:
        return self._names

    @property
    def dtypes(self) -> tuple[Any, ...]:
        return (np.dtype(self.value.dtype),) * self.shape[1]

    @property
    def schema(self) -> pd.DataFrame:
        return _schema_frame(self.names, self.dtypes)

    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        values = np.asarray(self.value)
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        projected = values[np.ix_(row_positions, column_positions)]
        return pd.DataFrame(projected, columns=[self.names[position] for position in column_positions])


@dataclass(frozen=True)
class SparseAdapter:
    value: Any
    _names: tuple[str, ...]

    @property
    def shape(self) -> tuple[int, int]:
        return int(self.value.shape[0]), int(self.value.shape[1])

    @property
    def names(self) -> tuple[str, ...]:
        return self._names

    @property
    def dtypes(self) -> tuple[Any, ...]:
        return (np.dtype(self.value.dtype),) * self.shape[1]

    @property
    def schema(self) -> pd.DataFrame:
        return _schema_frame(self.names, [pd.SparseDtype(dtype, fill_value=0) for dtype in self.dtypes])

    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        projected = self.value[row_positions, :][:, column_positions]
        return _sparse_frame(projected, [self.names[position] for position in column_positions])


def _sparse_frame(value: Any, names: Sequence[str]) -> pd.DataFrame:
    matrix = sparse.csc_matrix(value)
    from_spmatrix = cast(Any, pd.arrays.SparseArray).from_spmatrix
    columns = [pd.Series(from_spmatrix(matrix[:, position]), name=name) for position, name in enumerate(names)]
    return pd.concat(columns, axis=1) if columns else pd.DataFrame(columns=list(names))


def _restore_sorted_projection(
    values: Any,
    row_positions: np.ndarray,
    column_positions: np.ndarray,
) -> np.ndarray:
    """Read sorted-backed indices and restore duplicate/requested order."""
    row_unique, row_inverse = np.unique(row_positions, return_inverse=True)
    column_unique, column_inverse = np.unique(column_positions, return_inverse=True)
    n_rows, n_columns = int(values.shape[0]), int(values.shape[1])

    if np.array_equal(row_unique, np.arange(n_rows)):
        projected = np.asarray(values[:, column_unique])
    elif np.array_equal(column_unique, np.arange(n_columns)):
        projected = np.asarray(values[row_unique, :])
    else:
        rows = [np.asarray(values[int(row), column_unique]).reshape(1, -1) for row in row_unique]
        projected = np.concatenate(rows, axis=0) if rows else np.empty((0, len(column_unique)), dtype=values.dtype)
    return projected[row_inverse, :][:, column_inverse]


@dataclass(frozen=True)
class BackedDenseAdapter(DenseAdapter):
    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        projected = _restore_sorted_projection(self.value, row_positions, column_positions)
        return pd.DataFrame(projected, columns=[self.names[position] for position in column_positions])


@dataclass(frozen=True)
class BackedSparseAdapter(SparseAdapter):
    def read(self, row_positions: np.ndarray, column_positions: np.ndarray) -> pd.DataFrame:
        storage = str(getattr(self.value, "format", "csr")).lower()
        if "csc" in storage:
            projected = self.value[:, column_positions][row_positions, :]
        else:
            projected = self.value[row_positions, :][:, column_positions]
        if hasattr(projected, "to_memory"):
            projected = projected.to_memory()
        return _sparse_frame(projected, [self.names[position] for position in column_positions])


def adapter_for(value: Any, *, names: Sequence[str] | None = None) -> SourceAdapter:
    """Return the adapter matching an in-memory, view, or backed source."""
    if isinstance(value, pd.DataFrame):
        return PandasAdapter(value)
    shape = getattr(value, "shape", None)
    if shape is None:
        raise UnknownSourceError("source has no two-dimensional shape")
    width = int(shape[1] if len(shape) > 1 else 1)
    resolved_names = tuple(str(name) for name in (names if names is not None else range(width)))
    if len(resolved_names) != width:
        raise SelectionError("source column names do not match its width")
    if sparse.issparse(value):
        return SparseAdapter(value, resolved_names)
    module = type(value).__module__
    type_name = type(value).__name__.lower()
    if "sparse_dataset" in module or type_name in {"csrdataset", "cscdataset"}:
        return BackedSparseAdapter(value, resolved_names)
    if isinstance(value, np.ndarray):
        return DenseAdapter(value, resolved_names)
    return BackedDenseAdapter(value, resolved_names)


def source_adapter(adata: AnnData, source: str, *, key: str | None = None, layer: str | None = None) -> SourceAdapter:
    if source == "x":
        try:
            value = adata.layers[layer] if layer is not None else adata.X
        except KeyError as exc:
            raise UnknownSourceError(f"Unknown layer: {layer!r}") from exc
        return adapter_for(value, names=tuple(str(name) for name in adata.var_names))
    if source == "raw":
        if layer is not None:
            raise UnknownSourceError("raw source does not support layer")
        if adata.raw is None:
            raise UnknownSourceError("AnnData object has no raw matrix")
        return adapter_for(adata.raw.X, names=tuple(str(name) for name in adata.raw.var_names))
    if source in {"obsm", "varm", "obsp", "varp"}:
        if key is None:
            raise UnknownSourceError(f"{source} source requires a key")
        mapping = getattr(adata, source)
        try:
            value = mapping[key]
        except KeyError as exc:
            raise UnknownSourceError(f"Unknown {source} key: {key!r}") from exc
        if isinstance(value, pd.DataFrame):
            return adapter_for(value)
        if source == "obsp":
            names = tuple(str(name) for name in adata.obs_names)
        elif source == "varp":
            names = tuple(str(name) for name in adata.var_names)
        else:
            names = None
        return adapter_for(value, names=names)
    if source == "obs":
        return PandasAdapter(cast(pd.DataFrame, adata.obs))
    if source == "var":
        return PandasAdapter(cast(pd.DataFrame, adata.var))
    if source == "uns":
        if key is None:
            raise UnknownSourceError("uns source requires a key")
        try:
            value = adata.uns[key]
        except KeyError as exc:
            raise UnknownSourceError(f"Unknown uns key: {key!r}") from exc
        if isinstance(value, pd.Series):
            value = value.to_frame()
        elif isinstance(value, Mapping):
            value = pd.DataFrame(value)
        if not isinstance(value, pd.DataFrame):
            raise UnknownSourceError(f"uns key {key!r} cannot be represented as a pandas DataFrame")
        return PandasAdapter(value)
    raise UnknownSourceError(f"Unknown AnnData source: {source!r}")


def resolve_columns(adapter: SourceAdapter, request: Any, *, mode: str = "expression") -> np.ndarray:
    """Resolve exact source positions without reading source values."""
    names = list(adapter.names)
    if request is None:
        return np.arange(len(names), dtype=np.intp)
    if mode == "selection":
        if isinstance(request, str):
            return _name_positions(names, [request])
        if isinstance(request, Sequence) and all(isinstance(item, str) for item in request):
            return _name_positions(names, request)
    schema = adapter.schema
    dependencies = _request_dependencies(request, schema)
    if dependencies is None:
        return np.arange(len(names), dtype=np.intp)
    missing = [name for name in dependencies if name not in names and not name.startswith("__annplyr_")]
    if missing:
        raise UnknownColumnError(f"Unknown source column(s): {', '.join(sorted(missing))}")
    if mode == "selection":
        selected = _selection_names(request, schema)
        return _name_positions(names, selected)
    dependency_set = set(dependencies)
    return np.asarray([position for position, name in enumerate(names) if name in dependency_set], dtype=np.intp)


def _selection_names(request: Any, schema: pd.DataFrame) -> list[str]:
    from annplyr._frames import evaluate_select

    return [str(name) for name in evaluate_select(schema, request).columns]


def _name_positions(available: Sequence[str], selected: Sequence[str]) -> np.ndarray:
    positions_by_name: dict[str, list[int]] = {}
    for position, name in enumerate(available):
        positions_by_name.setdefault(name, []).append(position)
    offsets: dict[str, int] = {}
    positions: list[int] = []
    for name in selected:
        offset = offsets.get(name, 0)
        matches = positions_by_name.get(name, [])
        if offset >= len(matches):
            raise UnknownColumnError(f"Unknown or over-selected source column: {name!r}")
        positions.append(matches[offset])
        offsets[name] = offset + 1
    return np.asarray(positions, dtype=np.intp)


def _request_dependencies(request: Any, schema: pd.DataFrame) -> frozenset[str] | None:
    if isinstance(request, Desc):
        return _request_dependencies(request.expr, schema)
    if isinstance(request, str):
        return frozenset({request})
    if isinstance(request, AnnplyrExpr):
        return request.dependencies
    dependencies = expression_dependencies(request)
    if dependencies is None:
        return None
    if isinstance(request, Mapping):
        from annplyr._frames import expand_assignments

        expanded = expand_assignments(schema, request)
        return _merge_request_dependencies(expanded.values(), schema)
    if isinstance(request, Sequence) and not isinstance(request, (str, bytes)):
        return _merge_request_dependencies(request, schema)
    if hasattr(request, "to_expr"):
        return _request_dependencies(request.to_expr(schema), schema)
    if hasattr(request, "resolve"):
        return frozenset(_selection_names(request, schema))
    return dependencies


def request_dependencies(request: Any, schema: pd.DataFrame) -> frozenset[str] | None:
    """Resolve conservative request dependencies from schema alone."""
    return _request_dependencies(request, schema)


def _merge_request_dependencies(values: Any, schema: pd.DataFrame) -> frozenset[str] | None:
    merged: set[str] = set()
    for value in values:
        dependencies = _request_dependencies(value, schema)
        if dependencies is None:
            return None
        merged.update(dependencies)
    return frozenset(merged)


@dataclass(frozen=True)
class PlannedRead:
    adapter: SourceAdapter
    row_positions: np.ndarray
    column_positions: np.ndarray
    context: str
    charge: bool = True

    @property
    def projected_cells(self) -> int:
        return len(self.row_positions) * len(self.column_positions) if self.charge else 0

    def read(self, row_positions: np.ndarray | None = None) -> pd.DataFrame:
        """Read one validated projection while retaining positional order."""
        rows = self.row_positions if row_positions is None else row_positions
        if len(self.column_positions) == 0:
            return pd.DataFrame(index=pd.RangeIndex(len(rows)))
        frame = self.adapter.read(rows, self.column_positions)
        frame.index = pd.RangeIndex(len(frame))
        return frame


@dataclass(frozen=True)
class ChunkPlan:
    """Deterministic row or column chunks for one resolved projection."""

    request: PlannedRead
    target_values: int = DEFAULT_REDUCTION_CHUNK_VALUES

    def __post_init__(self) -> None:
        if self.target_values <= 0:
            raise AnnplyrError("chunk target must be positive")

    @property
    def rows_per_chunk(self) -> int:
        width = max(1, len(self.request.column_positions))
        return max(1, self.target_values // width)

    @property
    def chunk_count(self) -> int:
        rows = len(self.request.row_positions)
        return (rows + self.rows_per_chunk - 1) // self.rows_per_chunk

    @property
    def columns_per_chunk(self) -> int:
        rows = max(1, len(self.request.row_positions))
        return max(1, self.target_values // rows)

    @property
    def column_chunk_count(self) -> int:
        if len(self.request.row_positions) > self.target_values:
            return 0
        columns = len(self.request.column_positions)
        if columns == 0:
            return 1
        return (columns + self.columns_per_chunk - 1) // self.columns_per_chunk

    def row_chunks(self) -> Iterator[np.ndarray]:
        """Yield contiguous request-order row positions without copying them."""
        step = self.rows_per_chunk
        for start in range(0, len(self.request.row_positions), step):
            yield self.request.row_positions[start : start + step]

    def read_chunks(self) -> Iterator[tuple[np.ndarray, pd.DataFrame]]:
        """Read each planned chunk after the enclosing planner is validated."""
        for rows in self.row_chunks():
            yield rows, self.request.read(rows)

    def column_chunks(self) -> Iterator[np.ndarray]:
        """Yield full-row column batches when one column fits the target."""
        if len(self.request.row_positions) > self.target_values:
            return
        columns = self.request.column_positions
        if len(columns) == 0:
            yield columns
            return
        step = self.columns_per_chunk
        for start in range(0, len(columns), step):
            yield columns[start : start + step]

    def read_column_chunks(self) -> Iterator[tuple[np.ndarray, pd.DataFrame]]:
        """Read bounded full-row column batches in source-projection order."""
        for columns in self.column_chunks():
            request = PlannedRead(
                self.request.adapter,
                self.request.row_positions,
                columns,
                self.request.context,
                self.request.charge,
            )
            yield columns, request.read()


@dataclass
class RequestPlanner:
    """Plan and budget every source before performing the first adapter read.

    ``max_matrix_values`` counts cumulative logical projected cells delivered
    to expression evaluation or pandas export. Opaque expressions resolve to
    all source columns. It does not estimate bytes or ownership-only AnnData
    copying.
    """

    max_matrix_values: int | None = None
    requests: list[PlannedRead] = field(default_factory=list)

    def add(
        self,
        adapter: SourceAdapter,
        *,
        request: Any = None,
        mode: str = "expression",
        row_positions: Any = None,
        context: str,
        charge: bool = True,
    ) -> int:
        rows = _positions(row_positions, adapter.shape[0], dimension="row")
        columns = resolve_columns(adapter, request, mode=mode)
        self.requests.append(PlannedRead(adapter, rows, columns, context, charge))
        return len(self.requests) - 1

    @property
    def projected_cells(self) -> int:
        return sum(request.projected_cells for request in self.requests)

    def validate(self) -> None:
        if self.max_matrix_values is not None and self.max_matrix_values < 0:
            raise AnnplyrError("max_matrix_values must be non-negative or None")
        if self.max_matrix_values is not None and self.projected_cells > self.max_matrix_values:
            contexts = ", ".join(request.context for request in self.requests if request.projected_cells)
            raise AnnplyrError(
                f"planned source(s) ({contexts}) would materialize {self.projected_cells} matrix values, "
                f"which exceeds max_matrix_values={self.max_matrix_values}"
            )

    def execute(self) -> list[pd.DataFrame]:
        self.validate()
        return [request.read() for request in self.requests]

    def chunk_plan(self, token: int, *, target_values: int = DEFAULT_REDUCTION_CHUNK_VALUES) -> ChunkPlan:
        """Return one chunk plan after validating every cumulative request."""
        self.validate()
        try:
            request = self.requests[token]
        except IndexError as exc:
            raise AnnplyrError(f"unknown planned request token: {token}") from exc
        return ChunkPlan(request, target_values)
