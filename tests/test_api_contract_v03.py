from __future__ import annotations

import inspect
import re
from pathlib import Path

import annplyr
from annplyr._accessor import AnnplyrAccessor
from annplyr._grouped import GroupedAnnData

CONTRACT = Path(__file__).parents[1] / "docs" / "development" / "api-contract-v0.3.md"
ROW = re.compile(r"^\| `(?P<section>accessor|grouped|export)\.(?P<name>[^`]+)` \|(?P<body>.*)\|$")

REQUIRED_GROUPED_ADDITIONS = {
    "anti_join",
    "full_join",
    "inner_join",
    "left_join",
    "pipe",
    "relocate",
    "rename",
    "rename_with",
    "right_join",
    "select",
    "semi_join",
    "slice",
    "transmute",
}


def _public_methods(cls: type) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and (not name.startswith("_") or name == "__iter__")
    }


def _inventory() -> dict[str, dict[str, list[str]]]:
    rows: dict[str, dict[str, list[str]]] = {"accessor": {}, "grouped": {}, "export": {}}
    for line in CONTRACT.read_text().splitlines():
        match = ROW.match(line)
        if match is None:
            continue
        section = match.group("section")
        name = match.group("name")
        assert name not in rows[section], f"duplicate contract row: {section}.{name}"
        cells = [cell.strip() for cell in match.group("body").split("|")]
        rows[section][name] = cells
    return rows


def test_contract_inventory_covers_public_api() -> None:
    rows = _inventory()
    assert set(rows["accessor"]) == _public_methods(AnnplyrAccessor)
    expected_grouped = _public_methods(GroupedAnnData) | REQUIRED_GROUPED_ADDITIONS
    assert set(rows["grouped"]) == expected_grouped
    assert set(rows["export"]) == set(annplyr.__all__)


def test_contract_rows_freeze_every_required_dimension() -> None:
    rows = _inventory()
    for section, entries in rows.items():
        assert entries, f"missing {section} inventory"
        for name, cells in entries.items():
            assert len(cells) == 9, f"{section}.{name} must freeze all nine behavior dimensions"
            assert all(cells), f"{section}.{name} has an empty behavior dimension"
            assert not any(re.search(r"\b(?:TBD|TODO)\b", cell, re.IGNORECASE) for cell in cells)
