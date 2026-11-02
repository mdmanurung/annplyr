from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Literal

import narwhals as nw
import pandas as pd

from annplyr._expr import AnnplyrExpr, AnnplyrSelector
from annplyr._expr import Source as _Container
from annplyr._tidy import TidyAnnData

# Verb results are real ``TidyAnnData`` objects, so this alias is the concrete
# class rather than a typing-only stub. It is still an ``AnnData`` subclass, so
# ``isinstance`` checks keep passing for consumers.
AnnDataWithAnnplyr = TidyAnnData


type Axis = Literal["obs", "var"]
type Expression = str | AnnplyrExpr | nw.Expr
type Selector = (
    str
    | AnnplyrSelector
    | AnnplyrExpr
    | nw.Expr
    | _Container
    | Sequence[str | AnnplyrSelector | AnnplyrExpr | nw.Expr | _Container]
)
type Source = Literal["obs", "var", "x", "raw", "obsm", "varm", "obsp", "varp", "uns"]
type SourceSelectors = Mapping[str, Selector]
type JoinInput = pd.DataFrame | Mapping[str, Any]
type JoinBy = str | Sequence[str] | None
type JoinRelationship = Literal["many-to-one", "one-to-one", "one-to-many", "many-to-many"]
type JoinMultiple = Literal["error", "first", "all"]
type JoinUnmatched = Literal["drop", "error"]
type NaMatches = Literal["na", "never"]


__all__ = [
    "AnnDataWithAnnplyr",
    "Axis",
    "Expression",
    "JoinBy",
    "JoinInput",
    "JoinMultiple",
    "JoinRelationship",
    "JoinUnmatched",
    "NaMatches",
    "Selector",
    "Source",
    "SourceSelectors",
    "TidyAnnData",
]
