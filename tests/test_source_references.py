"""``source:name`` column reference resolution."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from narwhals import col as nw_col

import annplyr as ap
from annplyr import _resolve as R
from annplyr._errors import IncompatibleAxisError, UnknownSourceError


@pytest.fixture
def sourced() -> AnnData:
    rng = np.random.default_rng(0)
    adata = AnnData(
        X=rng.normal(size=(6, 4)).astype("float32"),
        obs=pd.DataFrame({"n_genes": [1, 2, 3, 4, 5, 6]}, index=[f"c{i}" for i in range(6)]),
        var=pd.DataFrame(index=["ENSG001:1", "MS4A1", "CD79A", "NKG7"]),
    )
    adata.layers["counts"] = (adata.X > 0).astype("float32")
    adata.raw = adata
    adata.obsm["X_pca"] = rng.normal(size=(6, 3))
    adata.varm["loadings"] = rng.normal(size=(4, 2))
    return adata


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("n_genes", (None, "n_genes")),
        ("raw:MS4A1", ("raw", "MS4A1")),
        ("x:LYZ", ("x", "LYZ")),
        ("var:highly_variable", ("var", "highly_variable")),
        ("obsm:X_pca", ("obsm", "X_pca")),
        ("uns:colors", ("uns", "colors")),
    ],
)
def test_reserved_qualifiers_split_without_an_adata_object(text: str, expected: tuple) -> None:
    assert R.split_qualifier(text) == expected


def test_colon_in_a_feature_name_is_never_treated_as_a_qualifier(sourced: AnnData) -> None:
    # ENSG-style identifiers are real feature names, not container prefixes.
    assert R.split_qualifier("ENSG001:1") == (None, "ENSG001:1")
    assert R.resolve_source(sourced, "ENSG001:1").kind == "obs"


def test_colon_feature_name_is_preserved_when_it_really_is_a_column(sourced: AnnData) -> None:
    source = R.resolve_source(sourced, "ENSG001:1")

    assert source.kind == "obs"
    assert source.references == (("ENSG001:1", "ENSG001:1"),)


@pytest.mark.parametrize(
    ("text", "kind", "key", "axis"),
    [
        ("raw:MS4A1", "raw", None, "obs"),
        ("x:MS4A1", "x", None, "obs"),
        ("counts:MS4A1", "layer", "counts", "obs"),
        ("X_pca:0", "obsm", "X_pca", "obs"),
        ("loadings:0", "varm", "loadings", "var"),
        ("var:chrom", "var", None, "var"),
        ("n_genes", "obs", None, "obs"),
    ],
)
def test_references_resolve_to_containers_and_axes(
    sourced: AnnData, text: str, kind: str, key: str | None, axis: str
) -> None:
    source = R.resolve_source(sourced, text)

    assert (source.kind, source.key, source.axis) == (kind, key, axis)


def test_reference_keeps_the_written_text_and_the_bare_name(sourced: AnnData) -> None:
    source = R.resolve_source(sourced, "raw:MS4A1")

    assert source.references == (("raw:MS4A1", "MS4A1"),)
    assert source.bare_names == ("MS4A1",)


def test_virtual_axis_names_resolve_to_the_correct_axis(sourced: AnnData) -> None:
    assert R.resolve_source(sourced, "__annplyr_obs_names__").axis == "obs"
    assert R.resolve_source(sourced, "__annplyr_var_names__").axis == "var"
    assert R.resolve_source(sourced, "__annplyr_var_names__").kind == "var"


def test_mistyped_qualifier_is_reported_with_the_known_containers(sourced: AnnData) -> None:
    with pytest.raises(UnknownSourceError) as info:
        R.resolve_source(sourced, "nosuch:MS4A1")

    message = str(info.value)
    assert "'nosuch' is not a known container" in message
    assert "raw:" in message
    assert "counts:" in message
    assert "X_pca:" in message


def test_unknown_layer_raises_unknown_source(sourced: AnnData) -> None:
    with pytest.raises(UnknownSourceError):
        R.resolve_source(sourced, "missing:MS4A1")


def test_unknown_key_raises_unknown_source(sourced: AnnData) -> None:
    with pytest.raises(UnknownSourceError):
        R.resolve_source(sourced, "X_nope:0")


def test_references_partition_into_one_source_per_container(sourced: AnnData) -> None:
    grouped = R.group_by_source(sourced, {"n_genes", "raw:MS4A1", "counts:CD79A", "X_pca:0", "var:chrom"})

    assert {source.context: sorted(refs) for source, refs in grouped.items()} == {
        "obs": ["n_genes"],
        "raw": ["raw:MS4A1"],
        "layer:counts": ["counts:CD79A"],
        "obsm:X_pca": ["X_pca:0"],
        "var": ["var:chrom"],
    }


def test_several_references_into_one_source_are_merged(sourced: AnnData) -> None:
    grouped = R.group_by_source(sourced, {"raw:MS4A1", "raw:CD79A"})

    assert len(grouped) == 1
    source = next(iter(grouped))
    assert source.kind == "raw"
    assert sorted(source.bare_names) == ["CD79A", "MS4A1"]


def test_sources_split_by_axis(sourced: AnnData) -> None:
    grouped = R.group_by_source(sourced, {"n_genes", "var:chrom", "loadings:0"})

    rows, columns = R.axis_sources(list(grouped))

    assert [source.context for source in rows] == ["obs"]
    assert sorted(source.context for source in columns) == ["var", "varm:loadings"]


def test_cross_source_expression_is_rejected_with_a_useful_message(sourced: AnnData) -> None:
    with pytest.raises(IncompatibleAxisError) as info:
        R.assert_single_source(sourced, {"score": ap.col("n_genes") + ap.col("raw:MS4A1")})

    message = str(info.value)
    assert "cannot mix columns from obs, raw" in message
    assert "mutate()" in message


def test_same_source_expression_is_accepted(sourced: AnnData) -> None:
    assert R.assert_single_source(sourced, {"score": ap.col("raw:MS4A1") + ap.col("raw:CD79A")}) is None


def test_literal_only_payload_is_accepted(sourced: AnnData) -> None:
    assert R.assert_single_source(sourced, {"flag": ap.lit(True)}) is None


def test_references_collect_through_nested_payloads(sourced: AnnData) -> None:
    references = R.collect_references({"a": (ap.col("q") > 1, {"r": ap.col("raw:S")})})

    assert references == {"q", "raw:S"}


def test_opaque_expression_yields_unknown_references(sourced: AnnData) -> None:
    assert R.collect_references({"a": nw_col("z") > 1}) is None
    assert R.collect_references((nw_col("z") > 1,)) is None


def test_opaque_expression_is_not_rejected_for_cross_source(sourced: AnnData) -> None:
    # References are unknowable, so the single-source check cannot fire.
    assert R.assert_single_source(sourced, {"a": nw_col("z") > 1}) is None


def test_metadata_frame_labels_columns_by_reference(sourced: AnnData) -> None:
    source = R.ResolvedSource(kind="obs", key=None, axis="obs", references=(("obs:n_genes", "n_genes"),))

    frame = R.metadata_frame(sourced, source)

    assert "obs:n_genes" in frame.columns
    assert "n_genes" not in frame.columns


def test_label_frame_renames_bare_column_to_the_written_reference() -> None:
    frame = pd.DataFrame({"MS4A1": [1, 2]})
    reference = R.ResolvedSource(kind="raw", key=None, axis="obs", references=(("raw:MS4A1", "MS4A1"),))

    labelled = R.label_frame(frame, reference)

    assert list(labelled.columns) == ["raw:MS4A1"]
    assert labelled["raw:MS4A1"].tolist() == [1, 2]
