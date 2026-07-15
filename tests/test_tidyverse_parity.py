from __future__ import annotations

import pandas as pd
import pytest
from anndata import AnnData

import annplyr as ap


def test_across_expands_assignments_in_mutate_and_summarize(dense_adata: AnnData) -> None:
    mutated = dense_adata.ap.mutate(
        obs=ap.across(
            ap.where(pd.api.types.is_numeric_dtype),
            lambda column: ap.col(column) * 10,
            names="{col}_x10",
        )
    )

    assert mutated.obs["score_x10"].tolist() == [10.0, 20.0, 5.0, 30.0, 25.0]

    summary = dense_adata.ap.summarize(
        obs=ap.across(
            ["score"],
            {"mean": ap.mean, "max": ap.max},
            names="{fn}_{col}",
        ),
        by="batch",
    )

    pd.testing.assert_frame_equal(
        summary.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "mean_score": [11 / 6, 7 / 4], "max_score": [2.5, 3.0]}),
    )

    transmuted = dense_adata.ap.transmute(
        obs=ap.across(["score"], lambda column: ap.col(column) + 1, names="{col}_plus1")
    )
    assert transmuted.obs.columns.tolist() == ["score_plus1"]
    assert transmuted.obs["score_plus1"].tolist() == [2.0, 3.0, 1.5, 4.0, 3.5]

    with pytest.raises(ap.DuplicateNameError, match="Duplicate assignment"):
        dense_adata.ap.mutate(obs={"score": ap.lit(0), "also_score": ap.across(["score"])})


def test_pick_if_any_and_if_all_use_tidyselect_columns(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["qc_mito"] = [0.1, 0.4, 0.2, 0.6, 0.3]
    adata.obs["qc_ribo"] = [0.2, 0.1, 0.3, 0.7, 0.2]

    picked = adata.ap.to_df(obs=ap.pick(ap.starts_with("qc_")))
    assert picked.columns.tolist() == ["qc_mito", "qc_ribo"]

    any_high = adata.ap.filter(obs=ap.if_any(ap.starts_with("qc_"), lambda column: ap.col(column) > 0.5))
    assert any_high.obs_names.tolist() == ["c3"]

    all_low = adata.ap.filter(obs=ap.if_all(ap.starts_with("qc_"), lambda column: ap.col(column) < 0.35))
    assert all_low.obs_names.tolist() == ["c0", "c2", "c4"]


def test_ntile_and_additional_rank_helpers(dense_adata: AnnData) -> None:
    mutated = dense_adata.ap.mutate(
        obs={
            "tile": ap.ntile("score", 2),
            "tile3": ap.ntile("score", 3),
            "max_rank": ap.max_rank("score"),
            "average_rank": ap.average_rank("score"),
        }
    )

    assert mutated.obs["tile"].tolist() == [1, 1, 1, 2, 2]
    assert mutated.obs["tile3"].tolist() == [1, 2, 1, 3, 2]
    assert mutated.obs["max_rank"].tolist() == [2, 3, 1, 5, 4]
    assert mutated.obs["average_rank"].tolist() == [2.0, 3.0, 1.0, 5.0, 4.0]


def test_additional_logic_and_recode_helpers(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["flag"] = [True, False, True, True, False]

    mutated = adata.ap.mutate(
        obs={
            "cummean_score": ap.cummean("score"),
            "ever_flag": ap.cumany("flag"),
            "all_flag": ap.cumall("flag"),
            "near_score": ap.near("score", 2.0, tolerance=0.01),
            "batch_label": ap.recode("batch", {"A": "alpha", "B": "beta"}),
            "batch_or_missing": ap.recode("batch", {"A": "alpha"}, default=None),
            "type_group": ap.case_match("cell_type", (["T", "B"], "lymphoid"), ("Mono", "myeloid")),
        }
    )

    assert mutated.obs["cummean_score"].round(3).tolist() == [1.0, 1.5, 1.167, 1.625, 1.8]
    assert mutated.obs["ever_flag"].tolist() == [True, True, True, True, True]
    assert mutated.obs["all_flag"].tolist() == [True, False, False, False, False]
    assert mutated.obs["near_score"].tolist() == [False, True, False, False, False]
    assert mutated.obs["batch_label"].tolist() == ["alpha", "alpha", "beta", "beta", "alpha"]
    assert mutated.obs["batch_or_missing"].isna().tolist() == [False, False, True, True, False]
    assert mutated.obs["type_group"].tolist() == ["lymphoid", "lymphoid", "lymphoid", "myeloid", "lymphoid"]


def test_grouped_matrix_sources_arrange_distinct_and_slices_are_group_local(dense_adata: AnnData) -> None:
    adata = dense_adata.copy()
    adata.obs["label"] = ["x", "x", "x", "x", "y"]
    grouped = dense_adata.ap.group_by(obs="batch")

    assert grouped.slice_min(ap.col("score"), n=1).obs_names.tolist() == ["c0", "c2"]
    assert grouped.slice_max(ap.col("score"), n=1).obs_names.tolist() == ["c4", "c3"]
    assert grouped.slice_tail(n=0).n_obs == 0

    sampled = grouped.slice_sample(n=1, random_state=7)
    assert sampled.n_obs == 2
    assert sampled.obs["batch"].tolist() == ["A", "B"]

    assert grouped.filter(x=ap.row_number() == 1).obs_names.tolist() == ["c0", "c2"]
    group_rows = grouped.mutate(x={"group_row": ap.row_number()})
    assert group_rows.obs["group_row"].tolist() == [1, 2, 1, 2, 3]
    assert grouped.arrange(x=ap.desc("g0")).obs_names.tolist() == ["c4", "c1", "c0", "c3", "c2"]
    assert adata.ap.group_by(obs="batch").distinct(obs="label", keep_all=True).obs_names.tolist() == ["c0", "c4", "c2"]


def test_slice_helpers_validate_edge_cases(dense_adata: AnnData) -> None:
    assert dense_adata.ap.slice_tail(n=0).n_obs == 0
    assert dense_adata.ap.slice_max([ap.col("score"), ap.col("batch")], n=2).obs_names.tolist() == ["c3", "c4"]

    with pytest.raises(ap.SelectionError, match="n and prop"):
        dense_adata.ap.slice_sample(n=1, prop=0.5)

    with pytest.raises(ap.SelectionError, match="larger"):
        dense_adata.ap.slice_sample(n=dense_adata.n_obs + 1, replace=False)

    with pytest.raises(ap.SelectionError, match="n and prop"):
        dense_adata.ap.group_by(obs="batch").slice_sample(n=1, prop=0.5)


def test_weighted_count_tally_and_add_tally(dense_adata: AnnData) -> None:
    counted = dense_adata.ap.count("batch", wt="score", sort=True, name="score_sum")
    pd.testing.assert_frame_equal(
        counted.reset_index(drop=True),
        pd.DataFrame({"batch": ["A", "B"], "score_sum": [5.5, 3.5]}),
    )

    tally = dense_adata.ap.tally(wt="score", name="score_sum")
    pd.testing.assert_frame_equal(tally, pd.DataFrame({"score_sum": [9.0]}))

    added = dense_adata.ap.add_count("batch", wt="score", name="batch_score_sum")
    assert added.obs["batch_score_sum"].tolist() == [5.5, 5.5, 3.5, 3.5, 5.5]

    grouped_added = dense_adata.ap.group_by(obs="batch").add_tally(wt="score", name="batch_score_sum")
    assert grouped_added.obs["batch_score_sum"].tolist() == [5.5, 5.5, 3.5, 3.5, 5.5]


def test_tidyr_style_dataframe_helpers() -> None:
    data = pd.DataFrame(
        {
            "sample": ["A_1", "B_2", None],
            "batch": ["A", None, "B"],
            "replicate": [1, 2, 3],
            "value": [1.0, None, 3.0],
        }
    )

    dropped = ap.drop_na(data, ["sample", "value"])
    assert dropped["sample"].tolist() == ["A_1"]

    filled = ap.fill(data, ["batch"])
    assert filled["batch"].tolist() == ["A", "A", "B"]

    separated = ap.separate(data, "sample", into=["condition", "unit"], sep="_", remove=False)
    assert separated[["condition", "unit"]].iloc[0].tolist() == ["A", "1"]

    united = ap.unite(separated, "sample_id", ["condition", "unit"], sep="-", na_rm=True)
    assert united["sample_id"].tolist() == ["A-1", "B-2", ""]

    extracted = ap.extract(
        pd.DataFrame({"sample": ["A-1"]}), "sample", into=["condition", "unit"], regex=r"([A-Z])-(\d)"
    )
    assert extracted.iloc[0].tolist() == ["A", "1"]

    row_split = ap.separate_rows(pd.DataFrame({"id": [1, 2], "genes": ["g1,g2", "g3"]}), "genes", sep=",")
    assert row_split["genes"].tolist() == ["g1", "g2", "g3"]

    nested = ap.nest(pd.DataFrame({"batch": ["A", "A", "B"], "score": [1, 2, 3]}), by="batch")
    assert nested["data"].map(len).tolist() == [2, 1]
    pd.testing.assert_frame_equal(
        ap.unnest(nested, "data"), pd.DataFrame({"batch": ["A", "A", "B"], "score": [1, 2, 3]})
    )

    chopped = ap.chop(pd.DataFrame({"id": [1, 1, 2], "tag": ["x", "y", "z"]}), "tag", by="id")
    assert chopped["tag"].tolist() == [["x", "y"], ["z"]]
    pd.testing.assert_frame_equal(ap.unchop(chopped, "tag"), pd.DataFrame({"id": [1, 1, 2], "tag": ["x", "y", "z"]}))

    longer = ap.unnest_longer(
        pd.DataFrame({"id": [1, 2], "vals": [[10, 20], []]}), "vals", indices_to="idx", keep_empty=True
    )
    assert longer["vals"].tolist()[:2] == [10, 20]
    assert pd.isna(longer["vals"].iloc[2])

    widened = ap.unnest_wider(pd.DataFrame({"id": [1], "info": [{"a": 10, "b": 20}]}), "info", names_sep="_")
    assert widened.columns.tolist() == ["id", "info_a", "info_b"]

    packed = ap.pack(pd.DataFrame({"id": [1], "a": [10], "b": [20]}), "metrics", ["a", "b"])
    assert packed["metrics"].iloc[0] == {"a": 10, "b": 20}
    unpacked = ap.unpack(packed, "metrics", names_sep="_")
    assert unpacked.columns.tolist() == ["id", "metrics_a", "metrics_b"]

    hoisted = ap.hoist(pd.DataFrame({"info": [{"sample": "s1", "n": 3}]}), "info", sample="sample", n="n")
    assert hoisted[["sample", "n"]].iloc[0].tolist() == ["s1", 3]


def test_unnest_preserves_inner_columns_when_all_nested_frames_are_empty() -> None:
    empty_inner = pd.DataFrame({"x": pd.Series([], dtype=float), "y": pd.Series([], dtype=str)})
    nested = pd.DataFrame({"id": ["a", "b"], "data": [empty_inner, empty_inner.copy()]})
    result = ap.unnest(nested, "data")
    assert result.columns.tolist() == ["id", "x", "y"]
    assert len(result) == 0


def test_separate_treats_sep_as_regex_and_handles_na() -> None:
    data = pd.DataFrame({"s": ["a1b", "c2d", None]})
    result = ap.separate(data, "s", into=["pre", "post"], sep=r"\d")
    assert result["pre"].tolist()[:2] == ["a", "c"]
    assert result["post"].tolist()[:2] == ["b", "d"]
    assert pd.isna(result["pre"].iloc[2])
    assert pd.isna(result["post"].iloc[2])
