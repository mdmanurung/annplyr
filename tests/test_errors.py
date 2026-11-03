from __future__ import annotations

import pytest
from anndata import AnnData

import annplyr as ap


def test_invalid_axis_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="axis"):
        dense_adata.ap.slice_head(axis="cells")


def test_summarize_rejects_mixed_axes_with_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="obs-axis or var-axis"):
        dense_adata.ap.summarize(
            {"cells": ap.n(), "genes": ap.n_distinct("var:__annplyr_var_names__")},
        )


def test_summarize_rejects_a_cross_container_aggregate(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="cannot mix columns from"):
        dense_adata.ap.summarize({"bad": ap.mean("score") + ap.mean("x:g0")})


def test_pull_takes_exactly_one_name(dense_adata: AnnData) -> None:
    with pytest.raises(TypeError):
        dense_adata.ap.pull()

    with pytest.raises(TypeError):
        dense_adata.ap.pull("batch", "x:g0")


def test_pull_requires_single_selected_column(dense_adata: AnnData) -> None:
    with pytest.raises(ap.SelectionError, match="exactly one column"):
        dense_adata.ap.pull(ap.all_of(["batch", "score"]))


def test_missing_metadata_column_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownColumnError, match="missing"):
        dense_adata.ap.select("missing")


def test_missing_obsm_key_raises_typed_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.UnknownSourceError, match="missing"):
        dense_adata.ap.pull("missing:0")


def test_group_by_rejects_mixed_axes(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="mix obs and var"):
        dense_adata.ap.group_by("batch", "var:feature_type")


def test_summarize_key_naming_the_other_axis_is_rejected(dense_adata: AnnData) -> None:
    """Resolving ``var:`` against an obs frame would group the wrong column."""
    with pytest.raises(ap.IncompatibleAxisError, match="names the var axis"):
        dense_adata.ap.summarize({"n": ap.n()}, by="var:feature_type")


def test_a_var_axis_reduce_requires_a_qualified_grouping_key(dense_adata: AnnData) -> None:
    with pytest.raises(ap.IncompatibleAxisError, match="bare name is obs metadata"):
        dense_adata.ap.summarize({"m": ap.max("var:feature_type")}, by="feature_type")


def test_pipe_keyword_collision_raises_annplyr_error(dense_adata: AnnData) -> None:
    def identity(data: AnnData) -> AnnData:
        return data

    with pytest.raises(ap.AnnplyrError, match="pipe target"):
        dense_adata.ap.pipe((identity, "data"), data=dense_adata)


@pytest.mark.parametrize(
    "call",
    [
        lambda a: a.ap.arrange("nope"),
        lambda a: a.ap.slice_max(ap.col("x:g0"), n=1),  # slice_* read metadata only
        lambda a: a.ap.slice_min(ap.col("nope"), n=1),
    ],
)
def test_unknown_sort_key_raises_unknown_column_error(dense_adata: AnnData, call) -> None:
    with pytest.raises(ap.UnknownColumnError):
        call(dense_adata)


def test_select_with_a_repeated_name_raises_duplicate_name_error(dense_adata: AnnData) -> None:
    with pytest.raises(ap.DuplicateNameError):
        dense_adata.ap.select("batch", "batch")


def test_mutate_without_a_dict_names_the_fix(dense_adata: AnnData) -> None:
    with pytest.raises(ap.SelectionError, match=r"mutate\(name=expression\)"):
        dense_adata.ap.mutate(ap.across(ap.everything()))


def test_grouped_summarize_has_no_by_because_group_by_sets_the_keys(dense_adata: AnnData) -> None:
    with pytest.raises(TypeError, match="by"):
        dense_adata.ap.group_by("batch").summarize({"n": ap.n()}, by="cell_type")


def test_a_name_given_as_dict_and_keyword_raises(dense_adata: AnnData) -> None:
    with pytest.raises(ap.DuplicateNameError, match="s"):
        dense_adata.ap.mutate({"s": ap.col("score")}, s=ap.col("score") * 2)


def test_error_messages_name_the_fix() -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd

    from annplyr.data import X, obs

    genes = [f"G{i}" for i in range(40)]
    adata = ad.AnnData(
        np.ones((3, 40)), obs=pd.DataFrame({"louvain": list("abc")}, index=list("xyz")), var=pd.DataFrame(index=genes)
    )
    adata.raw = ad.AnnData(np.ones((3, 41)), obs=adata.obs, var=pd.DataFrame(index=[*genes, "LYZ"]))
    for call in (lambda: adata.ap.mutate(l=X.LYZ), lambda: adata.ap.select(X.LYZ), lambda: adata.ap.pull(X.LYZ)):
        with pytest.raises(ap.UnknownColumnError, match=r"write raw\.LYZ"):
            call()
    with pytest.raises(ap.UnknownColumnError) as typo:
        adata.ap.filter(obs.louvian == "a")
    assert "did you mean louvain" in str(typo.value) and "__annplyr_" not in str(typo.value)
    with pytest.raises(ap.UnknownColumnError) as unknown:
        adata.ap.select(X.G99)
    assert len(str(unknown.value)) < 200  # close matches only, not every gene


def test_mutate_takes_a_bare_constant(dense_adata: AnnData) -> None:
    out = dense_adata.ap.mutate(flag=1, ratio=0.5, ok=True)
    assert out.obs["flag"].tolist() == [1] * dense_adata.n_obs
    assert out.obs["ok"].all() and (out.obs["ratio"] == 0.5).all()
