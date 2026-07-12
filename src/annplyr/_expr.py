from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import narwhals as nw
import narwhals.selectors as cs


@dataclass(frozen=True)
class Desc:
    expr: Any


def col(*names: str | Iterable[str]) -> nw.Expr:
    return nw.col(*names)


def lit(value: Any) -> nw.Expr:
    return nw.lit(value)


obs_names: nw.Expr = col("obs_names")
var_names: nw.Expr = col("var_names")


def desc(expr: str | nw.Expr) -> Desc:
    return Desc(col(expr) if isinstance(expr, str) else expr)


def starts_with(prefix: str):
    return cs.matches(f"^{re.escape(prefix)}")


def ends_with(suffix: str):
    return cs.matches(f"{re.escape(suffix)}$")


def contains(pattern: str, *, regex: bool = False):
    return cs.matches(pattern if regex else re.escape(pattern))


def matches(pattern: str):
    return cs.matches(pattern)


def everything():
    return cs.all()


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


def between(expr: str | nw.Expr, lower: Any, upper: Any, *, closed: str = "both") -> nw.Expr:
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
