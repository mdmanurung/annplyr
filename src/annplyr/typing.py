"""Consumer-facing type aliases for annplyr pipelines."""

from __future__ import annotations

from annplyr._grouped import GroupedAnnData
from annplyr._typing import (
    AnnDataWithAnnplyr,
    Axis,
    Expression,
    JoinBy,
    JoinInput,
    JoinMultiple,
    JoinRelationship,
    JoinUnmatched,
    NaMatches,
    Selector,
    Source,
    SourceSelectors,
)

type GroupedReturn = AnnDataWithAnnplyr | GroupedAnnData

__all__ = [
    "AnnDataWithAnnplyr",
    "Axis",
    "Expression",
    "GroupedReturn",
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
