from __future__ import annotations

import inspect
from pathlib import Path

import annplyr as ap
from annplyr import AnnplyrAccessor, GroupedAnnData


def _parameters(obj: object) -> dict[str, inspect.Parameter]:
    return dict(inspect.signature(obj).parameters)


def test_no_verb_takes_copy_or_inplace() -> None:
    callables = [
        value
        for cls in [AnnplyrAccessor, GroupedAnnData]
        for name, value in inspect.getmembers(cls, predicate=callable)
        if not name.startswith("_")
    ]
    callables += [getattr(ap, name) for name in ap.__all__ if callable(getattr(ap, name))]
    for obj in callables:
        if inspect.isclass(obj):
            continue
        assert not {"copy", "inplace"} & set(_parameters(obj)), obj

    for name in ["left_join", "inner_join", "right_join", "full_join", "semi_join", "anti_join"]:
        assert (
            _parameters(getattr(AnnplyrAccessor, name))["other"].annotation
            == _parameters(getattr(GroupedAnnData, name))["other"].annotation
        )


def test_every_public_callable_has_a_docstring() -> None:
    methods = [
        value
        for cls in [AnnplyrAccessor, GroupedAnnData]
        for name, value in inspect.getmembers(cls, predicate=callable)
        if not name.startswith("_") or name == "__iter__"
    ]
    exports = [getattr(ap, name) for name in ap.__all__ if callable(getattr(ap, name))]

    assert all(inspect.getdoc(obj) for obj in [*methods, *exports])


def test_every_public_method_has_an_explicit_return_annotation() -> None:
    methods = [
        value
        for cls in [AnnplyrAccessor, GroupedAnnData]
        for name, value in inspect.getmembers(cls, predicate=callable)
        if not name.startswith("_") or name == "__iter__"
    ]
    assert all(inspect.signature(method).return_annotation is not inspect.Signature.empty for method in methods)


def test_grouped_interface_does_not_leak_internal_properties() -> None:
    public_properties = {
        name
        for name, value in inspect.getmembers_static(GroupedAnnData)
        if isinstance(value, property) and not name.startswith("_")
    }
    assert public_properties == set()


def test_accessor_and_grouped_docstrings_cover_contract_sections() -> None:
    required = ["Parameters\n----------", "Returns\n-------", "Ownership\n---------", "Raises\n------"]
    for cls in [AnnplyrAccessor, GroupedAnnData]:
        for name, value in inspect.getmembers(cls, predicate=callable):
            if name.startswith("_") and name != "__iter__":
                continue
            doc = inspect.getdoc(value)
            assert doc is not None
            assert all(section in doc for section in required), (
                f"missing public docstring section: {cls.__name__}.{name}"
            )


def test_join_and_pipe_docstrings_match_typed_failure_and_ownership_contracts() -> None:
    for cls in [AnnplyrAccessor, GroupedAnnData]:
        for name in ["left_join", "inner_join", "right_join", "full_join", "semi_join", "anti_join"]:
            doc = inspect.getdoc(getattr(cls, name)) or ""
            assert "JoinRelationshipError" in doc
            assert "DuplicateNameError" in doc
        assert "callable-defined" in (inspect.getdoc(cls.pipe) or "").lower()


def test_docs_pages_are_executable_notebooks() -> None:
    assert "nb_execution_raise_on_error = True" in Path("docs/conf.py").read_text()
    for page in ["index", "from_scanpy", "from_r", "concepts", "tutorial", "pseudobulk", "api"]:
        assert "format_name: myst" in Path(f"docs/{page}.md").read_text()


def test_flat_verb_docstrings_are_dedented_for_sphinx() -> None:
    """Python 3.12 keeps docstring indentation; an indented numpy section breaks RST."""
    doc = ap.verbs.filter.__doc__ or ""
    assert "\nParameters\n----------\n" in doc
