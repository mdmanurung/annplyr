from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal, Protocol

import narwhals as nw
import pandas as pd

from annplyr._errors import UnknownColumnError


@dataclass(frozen=True)
class Desc:
    expr: Any


class AnnplyrSelector(Protocol):
    """Selector resolved by annplyr against a pandas frame."""

    def resolve(self, frame: pd.DataFrame, columns: list[str], public_columns: list[str]) -> list[str]:
        """Resolve selected column names."""


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


def col(*names: str | Iterable[str]) -> nw.Expr:
    return nw.col(*names)


def lit(value: Any) -> nw.Expr:
    return nw.lit(value)


obs_names: nw.Expr = col("obs_names")
var_names: nw.Expr = col("var_names")


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


def _expr(expr: str | nw.Expr) -> nw.Expr:
    return col(expr) if isinstance(expr, str) else expr


def n() -> nw.Expr:
    return nw.len()


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


def if_else(condition: nw.Expr, true: Any, false: Any) -> nw.Expr:
    return nw.when(condition).then(_literal_expr(true)).otherwise(_literal_expr(false))


def case_when(*cases: tuple[nw.Expr, Any], default: Any = None) -> nw.Expr:
    expr = _literal_expr(default)
    for condition, value in reversed(cases):
        expr = if_else(condition, value, expr)
    return expr


def row_number() -> nw.Expr:
    return col("__annplyr_row_number__")
