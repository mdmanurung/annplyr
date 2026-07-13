from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, Protocol

import narwhals as nw
import pandas as pd

from annplyr._errors import DuplicateNameError, UnknownColumnError

_MISSING = object()


@dataclass(frozen=True)
class Desc:
    expr: Any


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
            if self.predicate(frame[column]):
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
            return lit(False if self.how == "any" else True)
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


def col(*names: str | Iterable[str]) -> nw.Expr:
    return nw.col(*names)


def lit(value: Any) -> nw.Expr:
    return nw.lit(value)


obs_names: nw.Expr = col("__annplyr_obs_names__").alias("obs_names")
var_names: nw.Expr = col("__annplyr_var_names__").alias("var_names")


def desc(expr: str | nw.Expr) -> Desc:
    return Desc(col(expr) if isinstance(expr, str) else expr)


def starts_with(prefix: str):
    return _PatternSelector(f"^{re.escape(prefix)}")


def ends_with(suffix: str):
    return _PatternSelector(f"{re.escape(suffix)}$")


def contains(pattern: str, *, regex: bool = False):
    return _PatternSelector(pattern if regex else re.escape(pattern))


def matches(pattern: str):
    return _PatternSelector(pattern)


def everything():
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


def _expr(expr: str | nw.Expr) -> nw.Expr:
    return col(expr) if isinstance(expr, str) else expr


def n() -> nw.Expr:
    return nw.len()


def n_distinct(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).n_unique()


def mean(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).mean()


def median(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).median()


def sd(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).std()


def sum(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).sum()


def min(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).min()


def max(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).max()


def first(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).first()


def last(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).last()


def nth(expr: str | nw.Expr, n: int, *, default: Any = None) -> nw.Expr:
    def _take_nth(series: Any) -> Any:
        index = n if n >= 0 else series.len() + n
        if index < 0 or index >= series.len():
            return default
        return series.item(index)

    return _expr(expr).map_batches(_take_nth, returns_scalar=True)


def lead(expr: str | nw.Expr, n: int = 1, *, default: Any = None) -> nw.Expr:
    if n < 0:
        return lag(expr, -n, default=default)
    if n == 0:
        return _expr(expr)
    shifted = _expr(expr).shift(-n)
    if default is None:
        return shifted
    return nw.when(col("__annplyr_row_number__") > (nw.len() - n)).then(_literal_expr(default)).otherwise(shifted)


def lag(expr: str | nw.Expr, n: int = 1, *, default: Any = None) -> nw.Expr:
    if n < 0:
        return lead(expr, -n, default=default)
    if n == 0:
        return _expr(expr)
    shifted = _expr(expr).shift(n)
    if default is None:
        return shifted
    return nw.when(col("__annplyr_row_number__") <= n).then(_literal_expr(default)).otherwise(shifted)


def coalesce(*exprs: str | nw.Expr | Any) -> nw.Expr:
    values = [col(expr) if isinstance(expr, str) else _literal_expr(expr) for expr in exprs]
    return nw.coalesce(*values)


def na_if(expr: str | nw.Expr, value: Any) -> nw.Expr:
    base = _expr(expr)
    return nw.when(base == value).then(None).otherwise(base)


def replace_na(expr: str | nw.Expr, value: Any) -> nw.Expr:
    return _expr(expr).fill_null(value)


def is_na(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).is_null()


def min_rank(expr: str | nw.Expr, *, descending: bool = False) -> nw.Expr:
    return _expr(expr).rank("min", descending=descending)


def max_rank(expr: str | nw.Expr, *, descending: bool = False) -> nw.Expr:
    return _expr(expr).rank("max", descending=descending)


def average_rank(expr: str | nw.Expr, *, descending: bool = False) -> nw.Expr:
    return _expr(expr).rank("average", descending=descending)


def dense_rank(expr: str | nw.Expr, *, descending: bool = False) -> nw.Expr:
    return _expr(expr).rank("dense", descending=descending)


def percent_rank(expr: str | nw.Expr, *, descending: bool = False) -> nw.Expr:
    rank = min_rank(expr, descending=descending) - 1
    denominator = n() - 1
    return nw.when(denominator == 0).then(0).otherwise(rank / denominator)


def cume_dist(expr: str | nw.Expr, *, descending: bool = False) -> nw.Expr:
    rank = _expr(expr).rank("max", descending=descending)
    return rank / n()


def ntile(expr: str | nw.Expr, buckets: int) -> nw.Expr:
    if buckets < 1:
        msg = "ntile buckets must be a positive integer"
        raise UnknownColumnError(msg)
    return (((_expr(expr).rank("ordinal") - 1) * buckets) / n()).floor() + 1


def cum_sum(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cum_sum()


def cum_min(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cum_min()


def cum_max(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cum_max()


def cum_prod(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cum_prod()


def cummean(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cum_sum() / row_number()


def cumany(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cast(nw.Int64).cum_sum() > 0


def cumall(expr: str | nw.Expr) -> nw.Expr:
    return _expr(expr).cast(nw.Int64).cum_sum() == row_number()


def near(expr: str | nw.Expr, other: Any, *, tolerance: float = 1e-8) -> nw.Expr:
    return (_expr(expr) - _expr_or_literal(other)).abs() <= tolerance


def case_match(expr: str | nw.Expr, *cases: tuple[Any, Any], default: Any = None) -> nw.Expr:
    base = _expr(expr)
    out = _literal_expr(default)
    for values, replacement in reversed(cases):
        condition = base.is_in(list(values)) if isinstance(values, (list, tuple, set, frozenset)) else base == values
        out = nw.when(condition).then(_literal_expr(replacement)).otherwise(out)
    return out


def recode(expr: str | nw.Expr, mapping: Mapping[Any, Any], *, default: Any = _MISSING) -> nw.Expr:
    cases = tuple((value, replacement) for value, replacement in mapping.items())
    return case_match(expr, *cases, default=_expr(expr) if default is _MISSING else default)


def between(
    expr: str | nw.Expr,
    lower: Any,
    upper: Any,
    *,
    closed: Literal["left", "right", "none", "both"] = "both",
) -> nw.Expr:
    return _expr(expr).is_between(lower, upper, closed=closed)


def _literal_expr(value: Any) -> Any:
    return value if isinstance(value, nw.Expr) else lit(value)


def _expr_or_literal(value: Any) -> Any:
    return _expr(value) if isinstance(value, (str, nw.Expr)) else lit(value)


def if_else(condition: nw.Expr, true: Any, false: Any) -> nw.Expr:
    return nw.when(condition).then(_literal_expr(true)).otherwise(_literal_expr(false))


def case_when(*cases: tuple[nw.Expr, Any], default: Any = None) -> nw.Expr:
    expr = _literal_expr(default)
    for condition, value in reversed(cases):
        expr = if_else(condition, value, expr)
    return expr


def row_number() -> nw.Expr:
    return col("__annplyr_row_number__")
