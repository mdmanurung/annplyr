"""Install concise runtime docstrings for delegated public callables."""

from __future__ import annotations

import inspect
from collections.abc import Callable, Iterable, Mapping
from typing import Any

_AXIS_METHODS = {
    "arrange",
    "distinct",
    "filter",
    "select",
    "slice",
    "slice_head",
    "slice_max",
    "slice_min",
    "slice_sample",
    "slice_tail",
}
_JOIN_METHODS = {"anti_join", "full_join", "inner_join", "left_join", "right_join", "semi_join"}
_INPLACE_METHODS = {"add_count", "add_tally", "mutate", "relocate", "rename", "rename_with"}
_EXTRACTION_METHODS = {"as_frame", "nest_by", "pivot_longer", "pull", "to_df", "to_tidy"}
_TABLE_METHODS = {
    "count",
    "group_data",
    "group_keys",
    "group_vars",
    "summarise",
    "summarize",
    "tally",
}


def _parameter_lines(obj: Callable[..., Any]) -> list[str]:
    lines: list[str] = []
    for name, parameter in inspect.signature(obj).parameters.items():
        if name == "self":
            continue
        display = name
        if parameter.kind is inspect.Parameter.VAR_POSITIONAL:
            display = f"*{name}"
        elif parameter.kind is inspect.Parameter.VAR_KEYWORD:
            display = f"**{name}"
        lines.extend([display, "    See the annotated public signature."])
    return lines or ["None", "    This callable takes no public parameters."]


def _docstring(
    obj: Callable[..., Any],
    *,
    summary: str,
    returns: str,
    ownership: str,
    failures: tuple[str, ...],
) -> str:
    lines = [summary, "", "Parameters", "----------", *_parameter_lines(obj)]
    lines.extend(
        [
            "",
            "Returns",
            "-------",
            returns,
            "",
            "Ownership",
            "---------",
            ownership,
            "",
            "Raises",
            "------",
        ]
    )
    for failure in failures:
        lines.extend([failure, "    If the request violates the corresponding public contract."])
    return "\n".join(lines)


def _fill(
    obj: Callable[..., Any],
    *,
    summary: str,
    returns: str,
    ownership: str,
    failures: tuple[str, ...],
) -> None:
    if inspect.getdoc(obj):
        return
    obj.__doc__ = _docstring(
        obj,
        summary=summary,
        returns=returns,
        ownership=ownership,
        failures=failures,
    )


def _method_policy(name: str, *, grouped: bool) -> tuple[str, str, tuple[str, ...]]:
    if name in _JOIN_METHODS:
        result = "GroupedAnnData" if grouped else "AnnData"
        return (
            result,
            "Returns an independent result by default; ``copy=False`` may be a view or materialized.",
            (
                "SelectionError",
                "UnknownColumnError",
                "JoinRelationshipError",
                "DuplicateNameError",
                "AnnplyrError",
            ),
        )
    if name in _AXIS_METHODS:
        result = "GroupedAnnData" if grouped else "AnnData"
        return (
            result,
            "Returns an independent result by default; ``copy=False`` may be a view or materialized.",
            ("SelectionError", "UnknownColumnError", "UnknownSourceError", "IncompatibleAxisError", "AnnplyrError"),
        )
    if name in _INPLACE_METHODS:
        result = "GroupedAnnData" if grouped else "AnnData"
        return (
            result,
            "Returns a new result unless ``inplace=True`` requests exact input identity.",
            ("SelectionError", "DuplicateNameError", "SizeMismatchError", "AnnplyrError"),
        )
    if name == "transmute":
        result = "GroupedAnnData" if grouped else "AnnData"
        return (
            result,
            "Always returns an independent result and never mutates the input.",
            ("SelectionError", "UnknownSourceError", "SizeMismatchError", "AnnplyrError"),
        )
    if name in _EXTRACTION_METHODS or name in _TABLE_METHODS:
        return (
            "A pandas or extraction result described by the annotated signature.",
            "Leaves the input unchanged; matrix-backed requests may materialize projected values.",
            ("SelectionError", "UnknownSourceError", "NameRepairError", "SizeMismatchError", "AnnplyrError"),
        )
    if name == "ungroup":
        return (
            "AnnData",
            "Returns the wrapped AnnData and removes the grouping boundary.",
            ("AnnplyrError",),
        )
    if name == "pipe":
        return (
            "The callable-defined result.",
            "Passes the exact wrapped object; mutation and ownership are callable-defined.",
            ("AnnplyrError",),
        )
    return (
        "The result described by the annotated signature.",
        "Leaves the input unchanged unless the callable exposes and receives ``inplace=True``.",
        ("AnnplyrError",),
    )


def _install_methods(cls: type[Any], *, grouped: bool) -> None:
    for name, obj in inspect.getmembers(cls, predicate=callable):
        if name.startswith("_") and name != "__iter__":
            continue
        returns, ownership, failures = _method_policy(name, grouped=grouped)
        qualifier = "grouped " if grouped else ""
        _fill(
            obj,
            summary=f"Apply the {qualifier}``{name}`` operation.",
            returns=returns,
            ownership=ownership,
            failures=failures,
        )


def _export_policy(obj: Callable[..., Any]) -> tuple[str, str, str, tuple[str, ...]]:
    module = getattr(obj, "__module__", "")
    name = getattr(obj, "__name__", type(obj).__name__)
    if module.endswith("._expr"):
        return (
            f"Build the ``{name}`` expression or selector.",
            "An AnnplyrExpr, selector, or expression expansion.",
            "Builds lazy expression metadata and does not read or mutate AnnData.",
            ("SelectionError", "UnknownColumnError", "UnknownSourceError"),
        )
    if module.endswith("._utils"):
        return (
            f"Apply the ``{name}`` single-cell metadata utility.",
            "The AnnData, table, mapping, or diagnostic described by the signature.",
            "Readers leave input unchanged; writers are independent unless ``inplace=True``.",
            ("SelectionError", "DuplicateNameError", "SizeMismatchError", "JoinRelationshipError", "AnnplyrError"),
        )
    return (
        f"Apply the ``{name}`` dataframe rectangling helper.",
        "A new pandas object described by the signature.",
        "Leaves the input unchanged and operates after explicit tabular materialization.",
        ("SelectionError", "DuplicateNameError", "NameRepairError", "SizeMismatchError"),
    )


def install_public_docstrings(
    accessor_cls: type[Any],
    grouped_cls: type[Any],
    namespace: Mapping[str, Any],
    exports: Iterable[str],
) -> None:
    """Fill delegated public callables that do not define their own docstring."""
    _install_methods(accessor_cls, grouped=False)
    _install_methods(grouped_cls, grouped=True)
    for name in exports:
        obj = namespace[name]
        if not callable(obj) or inspect.isclass(obj):
            continue
        summary, returns, ownership, failures = _export_policy(obj)
        _fill(obj, summary=summary, returns=returns, ownership=ownership, failures=failures)
