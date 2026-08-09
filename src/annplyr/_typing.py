from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, Literal

import narwhals as nw
import pandas as pd
from anndata import AnnData

from annplyr._expr import AnnplyrExpr, AnnplyrSelector

if TYPE_CHECKING:
    from annplyr._accessor import AnnplyrAccessor

    class AnnDataWithAnnplyr(AnnData):
        """Static AnnData view whose registered ``ap`` namespace is visible."""

        @property
        def ap(self) -> AnnplyrAccessor: ...

else:
    # Runtime registration modifies AnnData itself. Keeping this as an alias makes
    # casts free and preserves isinstance behavior; type checkers see the class above.
    AnnDataWithAnnplyr = AnnData


type Axis = Literal["obs", "var"]
type Expression = str | AnnplyrExpr | nw.Expr
type Selector = str | AnnplyrSelector | AnnplyrExpr | nw.Expr | Sequence[str | AnnplyrSelector | AnnplyrExpr | nw.Expr]
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
]
