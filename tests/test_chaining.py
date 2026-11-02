"""Chainable results, the free-function surface, and AnnData interoperability."""

from __future__ import annotations

import inspect
import pickle
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

import annplyr as ap
from annplyr import AnnplyrAccessor
from annplyr._tidy import SHADOWED_VERBS, TidyAnnData
from annplyr.data import X


def _pipeline(adata: AnnData) -> TidyAnnData:
    """Enter the accessor once, then chain with methods only."""
    return (
        adata.ap.filter(ap.col("batch") == "A")
        .select(ap.everything(), X.g0, X.g1)
        .mutate({"gene_score": (ap.col("x:g0") + ap.col("x:g1")) / 2})
    )


def test_accessor_result_is_a_tidy_anndata(dense_adata: AnnData) -> None:
    result = dense_adata.ap.filter(ap.col("batch") == "A")

    assert isinstance(result, TidyAnnData)
    assert isinstance(result, AnnData)


def test_chain_runs_without_reentering_the_accessor(dense_adata: AnnData) -> None:
    summary = _pipeline(dense_adata).summarize({"cells": ap.n(), "mean_score": ap.mean("gene_score")}, by="cell_type")

    assert list(summary.columns) == ["cell_type", "cells", "mean_score"]
    assert summary["cells"].sum() == 3


def test_chained_steps_intermediate_to_the_right_shapes(dense_adata: AnnData) -> None:
    filtered = dense_adata.ap.filter(ap.col("batch") == "A")
    selected = filtered.select(ap.everything(), X.g0, X.g1)
    scored = _pipeline(dense_adata)

    assert filtered.n_obs == 3
    assert selected.n_obs == 3
    assert list(selected.var_names) == ["g0", "g1"]
    assert list(scored.var_names) == ["g0", "g1"]
    assert "gene_score" in scored.obs.columns


def test_free_function_surface_matches_the_chained_accessor(dense_adata: AnnData) -> None:
    from_functions = ap.mutate(
        ap.select(
            ap.filter(dense_adata, ap.col("batch") == "A"),
            ap.everything(),
            X.g0,
            X.g1,
        ),
        {"gene_score": (ap.col("x:g0") + ap.col("x:g1")) / 2},
    )
    chained = _pipeline(dense_adata)

    assert from_functions.obs.equals(chained.obs)
    assert from_functions.shape == chained.shape


def test_every_accessor_verb_has_a_matching_free_function() -> None:
    accessor_verbs = {
        name
        for name, _ in inspect.getmembers(AnnplyrAccessor, predicate=inspect.isfunction)
        if not name.startswith("_")
    }

    assert accessor_verbs == set(ap.verbs.FUNCTIONS)
    for name in accessor_verbs:
        assert getattr(ap, name) is ap.verbs.FUNCTIONS[name]
        assert name in ap.__all__


def test_free_function_signature_takes_adata_first() -> None:
    parameters = list(inspect.signature(ap.filter).parameters)

    assert parameters[0] == "adata"
    assert [p for p in parameters if p != "adata"] == [
        p for p in inspect.signature(AnnplyrAccessor.filter).parameters if p != "self"
    ]


def test_free_functions_expose_docstrings_and_annotations() -> None:
    assert ap.filter.__doc__ is not None
    assert "filter" in ap.filter.__doc__
    assert ap.filter.__annotations__["adata"] is AnnData


def test_result_passes_anndata_gates(dense_adata: AnnData) -> None:
    result = dense_adata.ap.filter(ap.col("batch") == "A")

    assert type(result) is not AnnData
    assert isinstance(result, AnnData)
    assert result.n_obs == 3
    assert result.n_vars == dense_adata.n_vars
    assert isinstance(result.obs, pd.DataFrame)
    assert result.shape[0] == result.n_obs


def test_anndata_methods_that_build_new_objects_stay_chainable(dense_adata: AnnData) -> None:
    wrapped = dense_adata.ap.filter(ap.col("batch") == "A")

    assert isinstance(wrapped[:2], TidyAnnData)
    assert isinstance(wrapped.copy(), TidyAnnData)
    assert isinstance(wrapped.to_memory(), TidyAnnData)
    if hasattr(AnnData, "concatenate"):  # removed in anndata 0.13
        assert isinstance(wrapped.concatenate(wrapped), TidyAnnData)


def test_rename_categories_stays_chainable(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"])
    wrapped = adata.ap.filter(ap.col("batch") == "A")
    renamed = [f"renamed_{c}" for c in wrapped.obs["cell_type"].cat.categories]

    # upstream rename_categories mutates in place and returns None
    assert wrapped.rename_categories("cell_type", renamed) is None

    assert all(str(c).startswith("renamed_") for c in wrapped.obs["cell_type"].cat.categories)
    # the mutated object is still a chainable result
    assert isinstance(wrapped, TidyAnnData)
    assert wrapped.count("cell_type")["n"].sum() == 3


def test_slicing_a_result_does_not_break_the_chain(dense_adata: AnnData) -> None:
    result = dense_adata.ap.filter(ap.col("batch") == "A")[:1]

    assert result.n_obs == 1
    assert result.count("cell_type")["n"].tolist() == [1]


def test_to_anndata_returns_a_plain_anndata(dense_adata: AnnData) -> None:
    result = dense_adata.ap.filter(ap.col("batch") == "A")
    plain = result.to_anndata()

    assert type(plain) is AnnData
    assert not isinstance(plain, TidyAnnData)
    assert plain.shape == result.shape


def test_upstream_anndata_to_df_is_not_shadowed(dense_adata: AnnData) -> None:
    wrapped = dense_adata.ap.filter(ap.col("batch") == "A")

    assert list(wrapped.to_df().columns) == list(dense_adata.var_names)
    assert list(wrapped.to_wide("cell_type", X.g0)) == ["cell_type", "g0"]


def test_no_accessor_verb_is_silently_shadowed() -> None:
    assert SHADOWED_VERBS == {}


def test_wrap_preserves_identity_for_inplace_verbs(sparse_adata: AnnData) -> None:
    result = sparse_adata.ap.mutate({"flag": ap.lit(True)}, inplace=True)

    assert result is sparse_adata
    assert isinstance(sparse_adata, TidyAnnData)
    assert "flag" in sparse_adata.obs.columns


def test_wrap_leaves_terminal_pandas_results_untouched(dense_adata: AnnData) -> None:
    counts = dense_adata.ap.count("cell_type")
    series = dense_adata.ap.pull("score")

    assert type(TidyAnnData.wrap(counts)) is pd.DataFrame
    assert type(TidyAnnData.wrap(series)) is pd.Series
    assert TidyAnnData.wrap(None) is None
    assert TidyAnnData.wrap(3) == 3


def test_wrap_is_idempotent(dense_adata: AnnData) -> None:
    once = TidyAnnData.wrap(dense_adata)

    assert TidyAnnData.wrap(once) is once


def test_view_results_remain_views(dense_adata: AnnData) -> None:
    view = dense_adata.ap.filter(ap.col("score") > 1.0, copy=False)

    assert view.is_view
    assert isinstance(view, TidyAnnData)


def test_sparse_format_survives_the_wrap(sparse_adata: AnnData) -> None:
    from scipy import sparse

    result = sparse_adata.ap.filter(ap.col("batch") == "A")

    assert sparse.issparse(result.X)
    assert sparse.issparse(result.layers["counts"])
    assert result.n_obs == 3


def test_backed_objects_stay_chainable(dense_adata: AnnData, tmp_path: Path) -> None:
    import anndata as ad

    path = tmp_path / "backed.h5ad"
    dense_adata.write_h5ad(path)

    backed = ad.read_h5ad(path, backed="r")
    result = backed.ap.filter(ap.col("batch") == "A")

    assert isinstance(result, TidyAnnData)
    assert result.n_obs == 3
    assert backed.isbacked
    assert backed.n_obs == 5
    backed.file.close()


def test_results_are_picklable(sparse_adata: AnnData) -> None:
    result = sparse_adata.ap.filter(ap.col("batch") == "A")

    restored = pickle.loads(pickle.dumps(result))

    assert isinstance(restored, TidyAnnData)
    assert restored.shape == result.shape
    assert list(restored.layers) == list(result.layers)


def test_grouped_chain_still_returns_a_grouped_result(dense_adata: AnnData) -> None:
    grouped = dense_adata.ap.group_by("cell_type")

    assert isinstance(grouped, ap.GroupedAnnData)
    assert isinstance(grouped.mutate({"rank": ap.row_number()}), ap.GroupedAnnData)
    assert type(grouped.ungroup()) is TidyAnnData


def test_grouped_chain_survives_group_boundaries(dense_adata: AnnData) -> None:
    summary = (
        dense_adata.ap.group_by("cell_type")
        .mutate({"rank": ap.row_number()})
        .filter(ap.col("rank") <= 1)
        .ungroup()
        .count("cell_type", sort=True)
    )

    assert list(summary.columns) == ["cell_type", "n"]
    assert summary["n"].tolist() == [1, 1, 1]


def test_accessor_namespace_still_available_on_results(dense_adata: AnnData) -> None:
    result = dense_adata.ap.filter(ap.col("batch") == "A")

    assert isinstance(result.ap, AnnplyrAccessor)
    assert result.ap.to_wide("cell_type").equals(result.to_wide("cell_type"))


def test_pipe_receives_the_wrapped_object(sparse_adata: AnnData) -> None:
    seen: list[object] = []

    def capture(value: AnnData) -> int:
        seen.append(value)
        return value.n_obs

    filtered = sparse_adata.ap.filter(ap.col("batch") == "A")

    assert filtered.pipe(capture) == filtered.n_obs
    assert isinstance(seen[0], TidyAnnData)
    assert seen[0] is filtered


def test_chained_result_still_works_with_anndata_operations(dense_adata: AnnData) -> None:
    result = dense_adata.ap.filter(ap.col("batch") == "A")
    plain = result[:, ["g0", "g1"]].copy()

    assert plain.shape == (3, 2)
    assert isinstance(plain, AnnData)


class _Subclassed(AnnData):
    """Stands in for a user-defined AnnData subclass."""

    def domain_method(self) -> str:
        return "domain"


class _SlottedSubclass(AnnData):
    __slots__ = ("label",)

    def __init__(self, *args: Any, label: str = "unset", **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.label = label

    def domain_method(self) -> str:
        return f"slotted:{self.label}"


def _subclass_instance(cls: type[AnnData], **kwargs: Any) -> Any:
    return cls(
        X=np.zeros((4, 2), dtype="float32"),
        obs=pd.DataFrame({"g": list("abcd")}, index=list("wxyz")),
        var=pd.DataFrame(index=["g0", "g1"]),
        **kwargs,
    )


def test_wrap_preserves_a_user_anndata_subclass_inplace() -> None:
    adata = _subclass_instance(_Subclassed)

    result = adata.ap.mutate({"z": ap.lit(1)}, inplace=True)

    assert result is adata
    assert isinstance(adata, _Subclassed)
    assert isinstance(adata, TidyAnnData)
    assert adata.domain_method() == "domain"
    # and it is still chainable
    assert adata.count("g")["n"].tolist() == [1, 1, 1, 1]


def test_wrap_does_not_crash_on_a_subclass_with_slots() -> None:
    adata = _subclass_instance(_SlottedSubclass, label="keepme")

    result = adata.ap.mutate({"z": ap.lit(1)}, inplace=True)

    assert result is adata
    assert adata.label == "keepme"
    assert isinstance(adata, TidyAnnData)
    assert adata.domain_method() == "slotted:keepme"


def test_composed_subclass_keeps_its_methods_through_copies() -> None:
    adata = _subclass_instance(_Subclassed)
    adata.ap.mutate({"z": ap.lit(1)}, inplace=True)

    copied = adata.copy()
    sliced = adata[:2]

    assert copied.domain_method() == "domain"
    assert sliced.domain_method() == "domain"
    assert isinstance(copied, TidyAnnData)
    assert isinstance(sliced, TidyAnnData)


def test_to_anndata_restores_a_carried_subclass() -> None:
    adata = _subclass_instance(_Subclassed)
    adata.ap.mutate({"z": ap.lit(1)}, inplace=True)

    plain = adata.to_anndata()

    assert type(plain) is _Subclassed
    assert plain.domain_method() == "domain"
    assert not hasattr(plain, "count")


def test_verbs_returning_new_objects_do_not_claim_the_subclass() -> None:
    # AnnData's own copy() and __getitem__ return a plain AnnData for
    # subclasses, so a derived result is a plain TidyAnnData, not a subclass.
    result = _subclass_instance(_Subclassed).ap.filter(ap.col("g") == "a")

    assert type(result) is TidyAnnData
    assert not isinstance(result, _Subclassed)
    assert hasattr(result, "count")


def test_install_verbs_is_idempotent() -> None:
    from annplyr._tidy import install_verbs

    before = dict(SHADOWED_VERBS)
    install_verbs(AnnplyrAccessor)

    assert SHADOWED_VERBS == before == {}
    assert hasattr(TidyAnnData, "filter")
