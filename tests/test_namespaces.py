from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap
from annplyr import IncompatibleAxisError, SelectionError, UnknownSourceError
from annplyr.data import X, X_pca, counts, obs, raw, var


@pytest.fixture
def adata(dense_adata: AnnData) -> AnnData:
    dense_adata.raw = dense_adata.copy()
    return ap.tidy(dense_adata)


def test_namespace_columns_equal_the_string_forms(adata: AnnData) -> None:
    pairs = [
        (obs.score, "score"),
        (counts.g0, "counts:g0"),
        (X_pca[0], "X_pca:0"),
        (raw["g1"], "raw:g1"),
        (X.g2, "x:g2"),
    ]
    for source, text in pairs:
        assert source._source_column.endswith(text) or source._source_column == text
        assert adata.pull(source).equals(adata.pull(text))


def test_namespace_selectors_equal_the_string_forms(adata: AnnData) -> None:
    by_name = adata.to_wide(obs.batch, counts.g0, X_pca, raw.g1)
    by_text = adata.to_wide("batch", "counts:g0", "X_pca:0", "X_pca:1", "raw:g1")
    assert by_name.equals(by_text.rename(columns={"X_pca:0": "X_pca_0", "X_pca:1": "X_pca_1"}))


def test_mistyped_layer_names_the_known_containers(adata: AnnData) -> None:
    with pytest.raises(UnknownSourceError, match="counts"):
        adata.pull(ap.layer("cuonts").g0)


def test_select_negation_drops_and_mutate_negation_is_numeric(adata: AnnData) -> None:
    assert adata.select(-obs.score).obs.columns.tolist() == ["batch", "cell_type"]
    assert adata.select(-X.g0).var_names.tolist() == ["g1", "g2", "g3"]
    mutated = adata.mutate({"neg": -ap.col("score")})
    assert mutated.obs["neg"].tolist() == (-adata.obs["score"]).tolist()


def test_negated_name_is_refused_where_it_is_not_a_drop(adata: AnnData) -> None:
    with pytest.raises(SelectionError, match="mutate"):
        adata.count(-obs.score)


def test_whole_matrix_needs_allow_all_features(adata: AnnData) -> None:
    with pytest.raises(SelectionError, match="allow_all_features"):
        adata.to_wide(raw)
    with pytest.raises(SelectionError, match="allow_all_features"):
        adata.pivot_longer(counts)
    assert adata.to_wide(raw, allow_all_features=True).shape == (5, 4)
    assert adata.to_wide(obs, X_pca).shape == (5, 5)


def test_var_names_infer_the_axis_and_mixing_raises(adata: AnnData) -> None:
    assert adata.count(var.feature_type)["n"].tolist() == [3, 1]
    assert adata.group_by(var.feature_type).group_vars() == ["feature_type"]
    with pytest.raises(IncompatibleAxisError):
        adata.group_by(obs.batch, var.feature_type)
    with pytest.raises(IncompatibleAxisError):
        adata.to_wide(obs.batch, var.chrom)


def test_computed_expression_key_asks_for_mutate(adata: AnnData) -> None:
    with pytest.raises(SelectionError, match="mutate"):
        adata.group_by(obs.score + 1)
    with pytest.raises(SelectionError, match="mutate"):
        adata.summarize({"n": ap.n()}, by=obs.score + 1)


def test_expression_keys_match_strings(adata: AnnData) -> None:
    pd.testing.assert_frame_equal(adata.count(obs.batch, obs.cell_type), adata.count("batch", "cell_type"))
    pd.testing.assert_frame_equal(
        adata.summarize({"n": ap.n()}, by=obs.batch), adata.summarize({"n": ap.n()}, by="batch")
    )
    assert adata.group_by(obs.batch).slice_max(obs.score, n=1).ungroup().obs_names.tolist() == ["c4", "c3"]
    assert adata.slice_min(obs.score, n=1).obs_names.tolist() == ["c2"]


def test_rename_takes_columns_and_mappings(adata: AnnData) -> None:
    assert adata.rename(kind=obs.cell_type).obs.columns.tolist() == ["batch", "score", "kind"]
    assert adata.rename({"G+": X.g0}).var_names.tolist() == ["G+", "g1", "g2", "g3"]
    assert adata.relocate(obs.cell_type, after=obs.batch).obs.columns.tolist() == ["batch", "cell_type", "score"]


def test_data_module_and_source_behaviour() -> None:
    from annplyr.data import anything

    assert anything.g0._source_column == "anything:g0"
    assert repr(raw) == "Source('raw')"
    assert not hasattr(ap, "mena")
    with pytest.raises(TypeError, match=r"counts\.x"):
        counts("x")
    with pytest.raises(AttributeError):
        _ = counts.__wrapped__
    with pytest.raises(TypeError):
        iter(counts)


def test_glimpse_prints_obs_and_container_names(adata: AnnData, capsys: pytest.CaptureFixture[str]) -> None:
    assert adata.glimpse() is None
    out = capsys.readouterr().out
    assert "cell_type" in out and "layers: counts" in out and "obsm: X_pca" in out and "raw: yes" in out


def test_a_layer_with_a_reserved_name_reads_the_layer_not_obs() -> None:
    import anndata as ad

    adata = ad.AnnData(
        np.arange(6.0).reshape(3, 2),
        obs=pd.DataFrame({"g0": [9.0, 9.0, 9.0]}, index=["c0", "c1", "c2"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    adata.layers["obs"] = adata.X * 100
    assert adata.ap.pull(ap.layer("obs").g0).tolist() == [0.0, 200.0, 400.0]
    assert adata.ap.pull("layers:obs:g0").tolist() == [0.0, 200.0, 400.0]
    assert adata.ap.aggregate_cells("g0", source=ap.layer("obs")).X.tolist() == [[600.0, 900.0]]
    with pytest.raises(ap.AnnplyrError):
        adata.ap.pull(ap.layer("cnts").g0)


def test_aggregate_cells_takes_a_source_object(dense_adata) -> None:
    from annplyr.data import X

    by_object = dense_adata.ap.aggregate_cells("batch", source=X)
    by_name = dense_adata.ap.aggregate_cells("batch", source="x")
    np.testing.assert_array_equal(np.asarray(by_object.X), np.asarray(by_name.X))


def test_var_keys_need_no_axis_and_joins_take_columns(dense_adata) -> None:
    top = dense_adata.ap.slice_max(var.length, n=2)
    assert top.n_vars == 2 and top.n_obs == dense_adata.n_obs
    meta = pd.DataFrame({"batch": dense_adata.obs["batch"].unique(), "site": "s"})
    by_string = dense_adata.ap.left_join(meta, by="batch")
    by_column = dense_adata.ap.left_join(meta, by=obs.batch)
    pd.testing.assert_frame_equal(by_column.obs, by_string.obs)
    long = dense_adata.ap.join_features(X.g0)
    assert ".abundance_X" in long.columns
