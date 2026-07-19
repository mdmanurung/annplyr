from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

import annplyr as ap


def test_sample_meta_extracts_sample_constant_metadata(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["site"] = ["blood", "blood", "marrow", "marrow", "blood"]

    sample_meta = ap.sample_meta(adata, "batch", include=["site"])

    pd.testing.assert_frame_equal(
        sample_meta,
        pd.DataFrame({"batch": ["A", "B"], "site": ["blood", "marrow"]}),
    )

    selected = ap.sample_meta(adata, "batch", include=ap.starts_with("site"))
    pd.testing.assert_frame_equal(selected, sample_meta)

    excluded = ap.sample_meta(adata, "batch", include=["site", "score"], exclude=["score"])
    pd.testing.assert_frame_equal(excluded, sample_meta)

    with pytest.raises(ap.AnnplyrError, match="not sample-constant"):
        ap.sample_meta(adata, "batch", include=["cell_type"])

    with pytest.raises(ap.SelectionError, match="existing obs metadata"):
        ap.sample_meta(adata, "batch", include=ap.obs_names)

    first = ap.sample_meta(adata, "batch", include=["cell_type"], on_conflict="first")
    pd.testing.assert_frame_equal(first, pd.DataFrame({"batch": ["A", "B"], "cell_type": ["T", "T"]}))


def test_sample_summary_and_add_sample_meta_use_existing_tidy_verbs(dense_adata: AnnData) -> None:
    summary = ap.sample_summary(
        dense_adata,
        "batch",
        obs={"mean_score": ap.mean("score")},
        x={"mean_g0": ap.mean("g0")},
    )

    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "mean_score": [11 / 6, 7 / 4], "mean_g0": [2.0, 2.5]}),
    )

    sample_table = pd.DataFrame({"sample_id": ["B", "A"], "treatment": ["drug", "vehicle"]})
    joined = ap.add_sample_meta(dense_adata, sample_table, sample="batch", by="sample_id")

    assert joined.obs_names.tolist() == dense_adata.obs_names.tolist()
    assert joined.obs["treatment"].tolist() == ["vehicle", "vehicle", "drug", "drug", "vehicle"]
    assert "treatment" not in dense_adata.obs

    with pytest.raises(ap.AnnplyrError, match="sample metadata"):
        ap.add_sample_meta(dense_adata, {"batch": "A", "treatment": "vehicle"}, sample="batch")


def test_feature_present_reports_missing_and_wrong_case_features(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.var["symbol"] = ["CD3E", "MS4A1", "LYZ", "MALAT1"]

    result = ap.feature_present(adata, ["g0", "G1", "bad"])
    assert result.found_features == ["g0"]
    assert result.missing_features == ["G1", "bad"]
    assert result.wrong_case_found_features == ["g1"]
    assert not result.all_found
    pd.testing.assert_frame_equal(
        result.to_frame(),
        pd.DataFrame(
            {
                "feature": ["g0", "G1", "bad"],
                "found": [True, False, False],
                "wrong_case_found": [pd.NA, "g1", pd.NA],
            }
        ),
    )

    single = ap.feature_present(adata, "g0")
    assert single.found_features == ["g0"]
    assert single.missing_features == []
    assert single.all_found

    symbols = ap.feature_present(adata, ["CD3E", "ms4a1", "missing"], feature_column="symbol")
    assert symbols.found_features == ["CD3E"]
    assert symbols.missing_features == ["ms4a1", "missing"]
    assert symbols.wrong_case_found_features == ["MS4A1"]

    with pytest.raises(ap.UnknownColumnError, match="feature column"):
        ap.feature_present(adata, ["CD3E"], feature_column="missing")


def test_name_utilities_preserve_alignment_and_reject_duplicates(dense_adata: AnnData) -> None:
    renamed_obs = ap.rename_obs_names(dense_adata, lambda name: f"sample_{name}")
    assert renamed_obs.obs_names.tolist() == [f"sample_c{i}" for i in range(5)]
    assert dense_adata.obs_names.tolist() == [f"c{i}" for i in range(5)]
    assert renamed_obs.obs.index.tolist() == renamed_obs.obs_names.tolist()

    renamed_var = ap.rename_var_names(dense_adata, {"gene0": "g0", "gene1": "g1"})
    assert renamed_var.var_names.tolist() == ["gene0", "gene1", "g2", "g3"]
    assert renamed_var.var.index.tolist() == renamed_var.var_names.tolist()
    assert renamed_var[:, "gene0"].X.tolist() == dense_adata[:, "g0"].X.tolist()

    stripped = ap.replace_name_suffix(renamed_obs, current_suffix="_c0", new_suffix="", axis="obs")
    assert stripped.obs_names[0] == "sample"

    prefixed = ap.add_name_prefix(dense_adata, "run1", axis="var", sep=":")
    assert prefixed.var_names.tolist() == ["run1:g0", "run1:g1", "run1:g2", "run1:g3"]

    with pytest.raises(ap.DuplicateNameError, match="Duplicate obs"):
        ap.rename_obs_names(dense_adata, lambda _name: "same")


def test_name_duplicates_reports_counts_and_positions() -> None:
    adata = AnnData(
        X=np.array([[1.0], [2.0], [3.0]]),
        obs=pd.DataFrame(index=["cell", "cell", "other"]),
        var=pd.DataFrame(index=["gene"]),
    )

    duplicates = ap.name_duplicates(adata, axis="obs")

    pd.testing.assert_frame_equal(
        duplicates,
        pd.DataFrame({"name": ["cell"], "count": [2], "positions": [[0, 1]]}),
    )


def test_palette_helpers_use_scanpy_uns_convention(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"], categories=["B", "Mono", "T"], ordered=True)

    colored = ap.store_palette(adata, "cell_type", {"B": "#1f77b4", "Mono": "#ff7f0e", "T": "#2ca02c"})

    assert "cell_type_colors" in colored.uns
    assert colored.uns["cell_type_colors"] == ["#1f77b4", "#ff7f0e", "#2ca02c"]
    assert "cell_type_colors" not in adata.uns
    assert ap.get_palette(colored, "cell_type") == {"B": "#1f77b4", "Mono": "#ff7f0e", "T": "#2ca02c"}

    sequence = ap.store_palette(adata, "batch", ["#000000", "#ffffff"])
    assert ap.get_palette(sequence, "batch") == {"A": "#000000", "B": "#ffffff"}

    with pytest.raises(ap.SizeMismatchError, match="Palette length"):
        ap.store_palette(adata, "batch", ["#000000"])


def test_mutating_utility_helpers_reject_backed_anndata(tmp_path, dense_adata: AnnData) -> None:
    path = tmp_path / "adata.h5ad"
    dense_adata.write_h5ad(path)
    backed = read_h5ad(path, backed="r")

    try:
        with pytest.raises(ap.AnnplyrError, match="backed"):
            ap.rename_obs_names(backed, lambda name: f"x_{name}")
        with pytest.raises(ap.AnnplyrError, match="backed"):
            ap.store_palette(backed, "batch", ["#000000", "#ffffff"])
    finally:
        backed.file.close()
