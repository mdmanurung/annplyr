from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import get_overloads

from anndata import AnnData

import annplyr as ap
import annplyr.typing as apt

ROOT = Path(__file__).parents[1]


def test_public_runtime_types_and_typing_aliases_are_complete(dense_adata: AnnData) -> None:
    assert ap.AnnplyrAccessor is type(dense_adata.ap)
    assert isinstance(ap.col("score"), ap.AnnplyrExpr)
    assert isinstance(dense_adata.ap.group_by(obs="batch"), ap.GroupedAnnData)
    assert dense_adata.ap.group_by() is dense_adata
    assert apt.AnnDataWithAnnplyr is AnnData
    assert set(apt.__all__) == {
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
    }


def test_group_by_overloads_distinguish_grouped_and_ungrouped_returns() -> None:
    overloads = get_overloads(ap.AnnplyrAccessor.group_by)
    assert len(overloads) == 3
    assert [overload.__annotations__["return"] for overload in overloads] == [
        "AnnDataWithAnnplyr",
        "GroupedAnnData",
        "GroupedAnnData",
    ]


def test_namespace_registration_is_import_order_independent() -> None:
    cases = [
        "import anndata; import annplyr; adata = anndata.AnnData()",
        "import annplyr; import anndata; adata = anndata.AnnData()",
    ]
    for setup in cases:
        code = f"{setup}; assert isinstance(adata.ap, annplyr.AnnplyrAccessor); assert adata.ap is adata.ap"
        result = subprocess.run(
            [sys.executable, "-c", code],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr


def test_generated_typing_contract_is_current() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "render_typing_contract.py"), "--check"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
