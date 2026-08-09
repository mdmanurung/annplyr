from __future__ import annotations

import inspect
import tomllib
from pathlib import Path

import annplyr as ap
from annplyr import AnnplyrAccessor, GroupedAnnData


def _parameters(obj: object) -> dict[str, inspect.Parameter]:
    return dict(inspect.signature(obj).parameters)


def test_v03_changed_signatures_are_consistent() -> None:
    for name in [
        "filter",
        "select",
        "arrange",
        "distinct",
        "slice",
        "slice_head",
        "slice_tail",
        "slice_min",
        "slice_max",
        "slice_sample",
        "left_join",
        "inner_join",
        "right_join",
        "full_join",
        "semi_join",
        "anti_join",
    ]:
        assert _parameters(getattr(AnnplyrAccessor, name))["copy"].default is True

    for name in ["left_join", "inner_join", "right_join", "full_join", "semi_join", "anti_join"]:
        assert (
            _parameters(getattr(AnnplyrAccessor, name))["other"].annotation
            == _parameters(getattr(GroupedAnnData, name))["other"].annotation
        )

    for name in ["mutate", "rename", "rename_with", "relocate", "add_count", "add_tally"]:
        parameters = _parameters(getattr(AnnplyrAccessor, name))
        assert "copy" not in parameters
        assert parameters["inplace"].default is False

    transmute = _parameters(AnnplyrAccessor.transmute)
    assert "copy" not in transmute
    assert "inplace" not in transmute

    for name in [
        "add_sample_meta",
        "rename_obs_names",
        "rename_var_names",
        "add_name_prefix",
        "replace_name_suffix",
        "store_palette",
    ]:
        parameters = _parameters(getattr(ap, name))
        assert "copy" not in parameters
        assert parameters["inplace"].default is False


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


def test_sphinx_doctest_configuration_and_examples_are_present() -> None:
    pyproject = tomllib.loads(Path("pyproject.toml").read_text())
    conf = Path("docs/conf.py").read_text()
    quickstart = Path("docs/quickstart.md").read_text()

    assert "sphinx.ext.doctest" in conf
    assert "doctest_global_setup" in conf
    assert "sphinx-build -M doctest" in pyproject["tool"]["hatch"]["envs"]["docs"]["scripts"]["doctest"]
    assert quickstart.count("```{testcode}") >= 5
