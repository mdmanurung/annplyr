from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, Protocol

import narwhals as nw
import numpy as np
import pandas as pd

from annplyr._errors import DuplicateNameError, SelectionError, UnknownColumnError

_MISSING = object()


@dataclass(frozen=True)
class Desc:
    expr: Any


Cardinality = Literal["row", "scalar", "unknown"]


def _merge_dependencies(*values: Any) -> frozenset[str] | None:
    dependencies: set[str] = set()
    for value in values:
        if isinstance(value, AnnplyrExpr):
            if value.dependencies is None:
                return None
            dependencies.update(value.dependencies)
        elif isinstance(value, nw.Expr):
            return None
        elif isinstance(value, Mapping):
            nested = _merge_dependencies(*value.values())
            if nested is None:
                return None
            dependencies.update(nested)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            nested = _merge_dependencies(*value)
            if nested is None:
                return None
            dependencies.update(nested)
    return frozenset(dependencies)


def _unwrap(value: Any) -> Any:
    if isinstance(value, AnnplyrExpr):
        return value.expr
    if isinstance(value, Mapping):
        return {key: _unwrap(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return tuple(_unwrap(item) for item in value)
    if isinstance(value, list):
        return [_unwrap(item) for item in value]
    return value


def _combined_cardinality(*values: Any) -> Cardinality:
    cardinalities = [value.cardinality for value in values if isinstance(value, AnnplyrExpr)]
    if any(cardinality == "unknown" for cardinality in cardinalities):
        return "unknown"
    if any(cardinality == "row" for cardinality in cardinalities):
        return "row"
    return "scalar"


@dataclass(frozen=True)
class AnnplyrExpr:
    """Narwhals expression plus conservative source-planning metadata."""

    expr: nw.Expr
    dependencies: frozenset[str] | None
    output_width: int | None = 1
    cardinality: Cardinality = "unknown"

    __array_priority__ = 1000

    def to_narwhals(self) -> nw.Expr:
        """Return the wrapped public Narwhals expression."""
        return self.expr

    def _binary(self, other: Any, method: str, *, reflected: bool = False) -> AnnplyrExpr:
        if reflected:
            reverse_method = f"__r{method[2:]}"
            result = getattr(self.expr, reverse_method)(_unwrap(other))
        else:
            result = getattr(self.expr, method)(_unwrap(other))
        return AnnplyrExpr(
            result,
            _merge_dependencies(self, other),
            1,
            _combined_cardinality(self, other),
        )

    def __add__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__add__")

    def __radd__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__add__", reflected=True)

    def __sub__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__sub__")

    def __rsub__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__sub__", reflected=True)

    def __mul__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__mul__")

    def __rmul__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__mul__", reflected=True)

    def __truediv__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__truediv__")

    def __rtruediv__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__truediv__", reflected=True)

    def __floordiv__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__floordiv__")

    def __rfloordiv__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__floordiv__", reflected=True)

    def __mod__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__mod__")

    def __rmod__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__mod__", reflected=True)

    def __pow__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__pow__")

    def __rpow__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__pow__", reflected=True)

    def __and__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__and__")

    def __rand__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__and__", reflected=True)

    def __or__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__or__")

    def __ror__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__or__", reflected=True)

    def __eq__(self, other: Any) -> AnnplyrExpr:  # type: ignore[override]
        return self._binary(other, "__eq__")

    def __ne__(self, other: Any) -> AnnplyrExpr:  # type: ignore[override]
        return self._binary(other, "__ne__")

    def __lt__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__lt__")

    def __le__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__le__")

    def __gt__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__gt__")

    def __ge__(self, other: Any) -> AnnplyrExpr:
        return self._binary(other, "__ge__")

    def __invert__(self) -> AnnplyrExpr:
        return AnnplyrExpr(~self.expr, self.dependencies, self.output_width, self.cardinality)

    def __neg__(self) -> AnnplyrExpr:
        return AnnplyrExpr(-self.expr, self.dependencies, self.output_width, self.cardinality)

    def __pos__(self) -> AnnplyrExpr:
        return self

    def __bool__(self) -> bool:
        raise TypeError("AnnplyrExpr cannot be used as a Python boolean; combine predicates with & or |")

    def __getattr__(self, name: str) -> Any:
        attribute = getattr(self.expr, name)
        if callable(attribute):
            return _expression_method(self, attribute, name)
        return _ExprNamespace(self, attribute)


@dataclass(frozen=True)
class _ExprNamespace:
    owner: AnnplyrExpr
    value: Any

    def __getattr__(self, name: str) -> Any:
        attribute = getattr(self.value, name)
        if callable(attribute):
            return _expression_method(self.owner, attribute, name)
        return _ExprNamespace(self.owner, attribute)


_SCALAR_METHODS = {"all", "any", "first", "last", "len", "max", "mean", "median", "min", "n_unique", "std", "sum"}


def _expression_method(owner: AnnplyrExpr, method: Callable[..., Any], name: str) -> Callable[..., Any]:
    def call(*args: Any, **kwargs: Any) -> Any:
        result = method(*(_unwrap(arg) for arg in args), **{key: _unwrap(value) for key, value in kwargs.items()})
        if not isinstance(result, nw.Expr):
            return result
        dependencies = _merge_dependencies(owner, args, kwargs)
        if name == "over":
            cardinality: Cardinality = "row"
        elif name in _SCALAR_METHODS or kwargs.get("returns_scalar") is True:
            cardinality = "scalar"
        else:
            cardinality = _combined_cardinality(owner, *args)
        return AnnplyrExpr(result, dependencies, owner.output_width, cardinality)

    return call


def to_narwhals(value: Any) -> Any:
    """Unwrap an annplyr expression; leave non-wrapper values unchanged."""
    return _unwrap(value)


def expression_dependencies(value: Any) -> frozenset[str] | None:
    """Return exact dependencies or ``None`` for an opaque expression."""
    if isinstance(value, AnnplyrExpr):
        return value.dependencies
    if isinstance(value, nw.Expr):
        return None
    if isinstance(value, Desc):
        return expression_dependencies(value.expr)
    if hasattr(value, "to_expr"):
        return None
    return frozenset()


class AnnplyrSelector(Protocol):
    """Selector resolved by annplyr against a pandas frame."""

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        """Resolve selected column names."""


def _public_columns(frame: pd.DataFrame) -> list[str]:
    virtual = set(frame.attrs.get("annplyr_virtual_columns", set())) | {
        "__annplyr_obs_names__",
        "__annplyr_var_names__",
        "__annplyr_row_number__",
    }
    return [str(column) for column in frame.columns if str(column) not in virtual]


def _all_columns(frame: pd.DataFrame) -> list[str]:
    return [str(column) for column in frame.columns]


def _resolve_columns(selector: Any, frame: pd.DataFrame) -> list[str]:
    columns = _all_columns(frame)
    public_columns = _public_columns(frame)
    if selector is None:
        return public_columns
    if isinstance(selector, str):
        selector = all_of(selector)
    elif isinstance(selector, Sequence) and not isinstance(selector, (str, bytes)):
        if all(isinstance(item, str) for item in selector):
            selector = all_of(selector)
    if hasattr(selector, "resolve"):
        return selector.resolve(frame, columns, public_columns)
    msg = "tidyselect helper requires a string, sequence of strings, or annplyr selector"
    raise UnknownColumnError(msg)


@dataclass(frozen=True)
class _NameSelector:
    names: tuple[str, ...]
    strict: bool = True

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        missing = [name for name in self.names if name not in columns]
        if missing and self.strict:
            msg = f"Unknown column(s): {', '.join(missing)}"
            raise UnknownColumnError(msg)
        return [name for name in self.names if name in columns]

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _EverythingSelector:
    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        return public_columns

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _PatternSelector:
    pattern: str

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        regex = re.compile(self.pattern)
        return [column for column in public_columns if regex.search(column)]

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _WhereSelector:
    predicate: Any

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        selected: list[str] = []
        for column in public_columns:
            probe = pd.Series([], dtype=frame[column].dtype, name=column)
            try:
                decision = self.predicate(probe)
            except Exception as exc:
                msg = (
                    "where() predicates receive a zero-length typed Series in v0.3; "
                    "value-dependent predicates are unsupported"
                )
                raise SelectionError(msg) from exc
            if not isinstance(decision, (bool, np.bool_)):
                msg = "where() predicates must return one bool from schema/dtype information"
                raise SelectionError(msg)
            if bool(decision):
                selected.append(column)
        return selected

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _LastColSelector:
    offset: int = 0

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        index = len(public_columns) - 1 - self.offset
        if index < 0 or index >= len(public_columns):
            msg = f"last_col offset {self.offset!r} is outside the available columns"
            raise UnknownColumnError(msg)
        return [public_columns[index]]

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _UnionSelector:
    selectors: tuple[AnnplyrSelector, ...]

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        resolved: list[str] = []
        for selector in self.selectors:
            for name in selector.resolve(frame, columns, public_columns):
                if name not in resolved:
                    resolved.append(name)
        return resolved

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((*self.selectors, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _IntersectionSelector:
    left: AnnplyrSelector
    right: AnnplyrSelector

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        right = set(self.right.resolve(frame, columns, public_columns))
        return [name for name in self.left.resolve(frame, columns, public_columns) if name in right]

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _ComplementSelector:
    selector: AnnplyrSelector

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        excluded = set(self.selector.resolve(frame, columns, public_columns))
        return [name for name in public_columns if name not in excluded]

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> AnnplyrSelector:
        return self.selector


@dataclass(frozen=True)
class Across:
    selector: Any
    fns: Any = None
    names: str | None = None

    def expand(self, frame: pd.DataFrame) -> dict[str, Any]:
        selected = _resolve_columns(self.selector, frame)
        functions = _normalize_across_functions(self.fns)
        output: dict[str, Any] = {}
        for column in selected:
            for function_name, function in functions:
                template = self.names or ("{col}" if len(functions) == 1 else "{col}_{fn}")
                name = template.format(col=column, fn=function_name)
                if name in output:
                    msg = f"across generated duplicate output name: {name!r}"
                    raise DuplicateNameError(msg)
                output[name] = function(column)
        return output


@dataclass(frozen=True)
class _PickSelector:
    selector: Any

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        return _resolve_columns(self.selector, frame)

    def __or__(self, other: AnnplyrSelector) -> _UnionSelector:
        return _UnionSelector((self, other))

    def __and__(self, other: AnnplyrSelector) -> _IntersectionSelector:
        return _IntersectionSelector(self, other)

    def __invert__(self) -> _ComplementSelector:
        return _ComplementSelector(self)


@dataclass(frozen=True)
class _IfAnyAll:
    selector: Any
    predicate: Callable[[str], Any]
    how: str

    def to_expr(self, frame: pd.DataFrame) -> Any:
        selected = _resolve_columns(self.selector, frame)
        if not selected:
            row = col("__annplyr_row_number__")
            return row < 0 if self.how == "any" else row > 0
        expr = self.predicate(selected[0])
        for column in selected[1:]:
            next_expr = self.predicate(column)
            expr = (expr | next_expr) if self.how == "any" else (expr & next_expr)
        return expr


def _normalize_across_functions(fns: Any) -> list[tuple[str, Callable[[str], Any]]]:
    if fns is None:
        return [("", col)]
    if isinstance(fns, Mapping):
        functions = [(str(name), fn) for name, fn in fns.items()]
        if all(callable(function) for _, function in functions):
            return functions
    if callable(fns):
        return [(_function_label(fns, 1), fns)]
    if isinstance(fns, Sequence) and not isinstance(fns, (str, bytes)):
        functions = [(_function_label(fn, i), fn) for i, fn in enumerate(fns, start=1)]
        if all(callable(function) for _, function in functions):
            return functions
    msg = "across fns must be a callable, sequence of callables, mapping, or None"
    raise UnknownColumnError(msg)


def _function_label(function: Any, index: int) -> str:
    name = getattr(function, "__name__", "")
    if not name or name == "<lambda>":
        return f"fn{index}"
    return name


def col(*names: str | Iterable[str]) -> AnnplyrExpr:
    flattened: list[str] = []
    for name in names:
        if isinstance(name, str):
            flattened.append(name)
        else:
            flattened.extend(str(item) for item in name)
    return AnnplyrExpr(nw.col(*flattened), frozenset(flattened), len(flattened), "row")


def lit(value: Any) -> AnnplyrExpr:
    return AnnplyrExpr(nw.lit(value), frozenset(), 1, "scalar")


obs_names: AnnplyrExpr = col("__annplyr_obs_names__").alias("obs_names")
var_names: AnnplyrExpr = col("__annplyr_var_names__").alias("var_names")


def desc(expr: str | AnnplyrExpr | nw.Expr) -> Desc:
    return Desc(col(expr) if isinstance(expr, str) else expr)


def starts_with(prefix: str) -> AnnplyrSelector:
    return _PatternSelector(f"^{re.escape(prefix)}")


def ends_with(suffix: str) -> AnnplyrSelector:
    return _PatternSelector(f"{re.escape(suffix)}$")


def contains(pattern: str, *, regex: bool = False) -> AnnplyrSelector:
    return _PatternSelector(pattern if regex else re.escape(pattern))


def matches(pattern: str) -> AnnplyrSelector:
    return _PatternSelector(pattern)


def everything() -> AnnplyrSelector:
    return _EverythingSelector()


def all_of(names: str | Iterable[str]) -> AnnplyrSelector:
    if isinstance(names, str):
        names = [names]
    return _NameSelector(tuple(names), strict=True)


def any_of(names: str | Iterable[str]) -> AnnplyrSelector:
    if isinstance(names, str):
        names = [names]
    return _NameSelector(tuple(names), strict=False)


def where(predicate: Any) -> AnnplyrSelector:
    return _WhereSelector(predicate)


def last_col(offset: int = 0) -> AnnplyrSelector:
    return _LastColSelector(offset=offset)


def num_range(prefix: str, range: Iterable[int], *, width: int = 0) -> AnnplyrSelector:
    names = tuple(f"{prefix}{number:0{width}d}" if width else f"{prefix}{number}" for number in range)
    return _NameSelector(names, strict=True)


def pick(selector: Any) -> AnnplyrSelector:
    return _PickSelector(selector)


def across(selector: Any, fns: Any = None, *, names: str | None = None) -> Any:
    return Across(selector=selector, fns=fns, names=names)


def if_any(selector: Any, predicate: Callable[[str], Any]) -> Any:
    return _IfAnyAll(selector=selector, predicate=predicate, how="any")


def if_all(selector: Any, predicate: Callable[[str], Any]) -> Any:
    return _IfAnyAll(selector=selector, predicate=predicate, how="all")


def _expr(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    if isinstance(expr, str):
        return col(expr)
    if isinstance(expr, AnnplyrExpr):
        return expr
    if isinstance(expr, nw.Expr):
        return AnnplyrExpr(expr, None, None, "unknown")
    msg = "expression must be a column name, AnnplyrExpr, or Narwhals expression"
    raise UnknownColumnError(msg)


def n() -> AnnplyrExpr:
    return AnnplyrExpr(nw.len(), frozenset(), 1, "scalar")


def n_distinct(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).n_unique()


def mean(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).mean()


def median(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).median()


def sd(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).std()


def sum(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).sum()


def min(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).min()


def max(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).max()


def first(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).first()


def last(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).last()


def nth(expr: str | AnnplyrExpr | nw.Expr, n: int, *, default: Any = None) -> AnnplyrExpr:
    def _take_nth(series: Any) -> Any:
        index = n if n >= 0 else series.len() + n
        if index < 0 or index >= series.len():
            return default
        return series.item(index)

    return _expr(expr).map_batches(_take_nth, returns_scalar=True)


def lead(expr: str | AnnplyrExpr | nw.Expr, n: int = 1, *, default: Any = None) -> AnnplyrExpr:
    if n < 0:
        return lag(expr, -n, default=default)
    if n == 0:
        return _expr(expr)
    shifted = _expr(expr).shift(-n)
    if default is None:
        return shifted
    length = AnnplyrExpr(nw.len(), frozenset(), 1, "scalar")
    return if_else(col("__annplyr_row_number__") > (length - n), default, shifted)


def lag(expr: str | AnnplyrExpr | nw.Expr, n: int = 1, *, default: Any = None) -> AnnplyrExpr:
    if n < 0:
        return lead(expr, -n, default=default)
    if n == 0:
        return _expr(expr)
    shifted = _expr(expr).shift(n)
    if default is None:
        return shifted
    return if_else(col("__annplyr_row_number__") <= n, default, shifted)


def coalesce(*exprs: str | AnnplyrExpr | nw.Expr | Any) -> AnnplyrExpr:
    values = [
        col(expr) if isinstance(expr, str) else _expr(expr) if isinstance(expr, (AnnplyrExpr, nw.Expr)) else lit(expr)
        for expr in exprs
    ]
    return AnnplyrExpr(
        nw.coalesce(*[_unwrap(value) for value in values]),
        _merge_dependencies(*values),
        1,
        _combined_cardinality(*values),
    )


def na_if(expr: str | AnnplyrExpr | nw.Expr, value: Any) -> AnnplyrExpr:
    base = _expr(expr)
    return if_else(base == value, None, base)


def replace_na(expr: str | AnnplyrExpr | nw.Expr, value: Any) -> AnnplyrExpr:
    return _expr(expr).fill_null(value)


def is_na(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).is_null()


def min_rank(expr: str | AnnplyrExpr | nw.Expr, *, descending: bool = False) -> AnnplyrExpr:
    return _expr(expr).rank("min", descending=descending)


def max_rank(expr: str | AnnplyrExpr | nw.Expr, *, descending: bool = False) -> AnnplyrExpr:
    return _expr(expr).rank("max", descending=descending)


def average_rank(expr: str | AnnplyrExpr | nw.Expr, *, descending: bool = False) -> AnnplyrExpr:
    return _expr(expr).rank("average", descending=descending)


def dense_rank(expr: str | AnnplyrExpr | nw.Expr, *, descending: bool = False) -> AnnplyrExpr:
    return _expr(expr).rank("dense", descending=descending)


def percent_rank(expr: str | AnnplyrExpr | nw.Expr, *, descending: bool = False) -> AnnplyrExpr:
    rank = min_rank(expr, descending=descending) - 1
    denominator = n() - 1
    return if_else(denominator == 0, 0, rank / denominator)


def cume_dist(expr: str | AnnplyrExpr | nw.Expr, *, descending: bool = False) -> AnnplyrExpr:
    rank = _expr(expr).rank("max", descending=descending)
    return rank / n()


def ntile(expr: str | AnnplyrExpr | nw.Expr, buckets: int) -> AnnplyrExpr:
    if buckets < 1:
        msg = "ntile buckets must be a positive integer"
        raise UnknownColumnError(msg)
    return (((_expr(expr).rank("ordinal") - 1) * buckets) / n()).floor() + 1


def cum_sum(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cum_sum()


def cum_min(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cum_min()


def cum_max(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cum_max()


def cum_prod(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cum_prod()


def cummean(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cum_sum() / row_number()


def cumany(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cast(nw.Int64).cum_sum() > 0


def cumall(expr: str | AnnplyrExpr | nw.Expr) -> AnnplyrExpr:
    return _expr(expr).cast(nw.Int64).cum_sum() == row_number()


def near(expr: str | AnnplyrExpr | nw.Expr, other: Any, *, tolerance: float = 1e-8) -> AnnplyrExpr:
    return (_expr(expr) - _expr_or_literal(other)).abs() <= tolerance


def case_match(expr: str | AnnplyrExpr | nw.Expr, *cases: tuple[Any, Any], default: Any = None) -> AnnplyrExpr:
    base = _expr(expr)
    out = _literal_expr(default)
    for values, replacement in reversed(cases):
        condition = base.is_in(list(values)) if isinstance(values, (list, tuple, set, frozenset)) else base == values
        out = if_else(condition, replacement, out)
    return out


def recode(expr: str | AnnplyrExpr | nw.Expr, mapping: Mapping[Any, Any], *, default: Any = _MISSING) -> AnnplyrExpr:
    cases = tuple((value, replacement) for value, replacement in mapping.items())
    return case_match(expr, *cases, default=_expr(expr) if default is _MISSING else default)


def between(
    expr: str | AnnplyrExpr | nw.Expr,
    lower: Any,
    upper: Any,
    *,
    closed: Literal["left", "right", "none", "both"] = "both",
) -> AnnplyrExpr:
    return _expr(expr).is_between(lower, upper, closed=closed)


def _literal_expr(value: Any) -> Any:
    return value if isinstance(value, (AnnplyrExpr, nw.Expr)) else lit(value)


def _expr_or_literal(value: Any) -> Any:
    return _expr(value) if isinstance(value, (str, AnnplyrExpr, nw.Expr)) else lit(value)


def if_else(condition: AnnplyrExpr | nw.Expr, true: Any, false: Any) -> AnnplyrExpr:
    true_expr = _literal_expr(true)
    false_expr = _literal_expr(false)
    result = nw.when(_unwrap(condition)).then(_unwrap(true_expr)).otherwise(_unwrap(false_expr))
    return AnnplyrExpr(
        result,
        _merge_dependencies(condition, true_expr, false_expr),
        1,
        _combined_cardinality(condition, true_expr, false_expr),
    )


def case_when(*cases: tuple[AnnplyrExpr | nw.Expr, Any], default: Any = None) -> AnnplyrExpr:
    expr = _literal_expr(default)
    for condition, value in reversed(cases):
        expr = if_else(condition, value, expr)
    return expr


def row_number() -> AnnplyrExpr:
    return col("__annplyr_row_number__")
