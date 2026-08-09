"""Group keys stay attached to their own aggregates across summary sources.

`summarize()` evaluates each requested source separately and then combines the
pieces. Those pieces are keyed by group identity rather than by key values, so
a combining step that reorders rows cannot hand a group another group's
aggregates. Every fixture here deliberately uses groups whose first-seen order
differs from their sorted order, which is the case the bug needed.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap

# First-seen order is zeta, alpha, mu - deliberately not the sorted order.
LABELS = ["zeta", "alpha", "mu", "zeta", "alpha", "zeta"]
VALUES = np.array(
    [
        [10.0, 1.0],
        [20.0, 2.0],
        [30.0, 3.0],
        [40.0, 4.0],
        [50.0, 5.0],
        [60.0, 6.0],
    ],
    dtype=np.float32,
)


def make_adata(matrix: Any = None) -> AnnData:
    adata = AnnData(
        X=VALUES.copy() if matrix is None else matrix,
        obs=pd.DataFrame(
            {
                "label": LABELS,
                "score": np.arange(1.0, 7.0, dtype=np.float32),
            },
            index=[f"cell_{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.raw = adata
    return adata


def expected() -> pd.DataFrame:
    frame = pd.DataFrame({"label": LABELS, "score": np.arange(1.0, 7.0), "g1": VALUES[:, 0]})
    grouped = frame.groupby("label", sort=False)
    return pd.DataFrame(
        {
            "label": list(dict.fromkeys(LABELS)),
            "cells": grouped.size().reindex(dict.fromkeys(LABELS)).to_numpy(),
            "mean_score": grouped["score"].mean().reindex(dict.fromkeys(LABELS)).to_numpy(),
            "mean_g1": grouped["g1"].mean().reindex(dict.fromkeys(LABELS)).to_numpy(),
        }
    )


@pytest.mark.parametrize("source", ["x", "raw"])
def test_metadata_and_matrix_aggregates_stay_on_their_own_group(source: str) -> None:
    adata = make_adata()

    summary = adata.ap.summarize(
        obs={"cells": ap.n(), "mean_score": ap.mean("score")},
        **{source: {"mean_g1": ap.mean("g1")}},
        by="label",
    )

    reference = expected()
    assert summary["label"].tolist() == reference["label"].tolist()
    assert summary["cells"].tolist() == reference["cells"].tolist()
    assert summary["mean_score"].tolist() == pytest.approx(reference["mean_score"].tolist())
    assert summary["mean_g1"].tolist() == pytest.approx(reference["mean_g1"].tolist())


def test_persistent_grouping_agrees_with_one_off_by() -> None:
    adata = make_adata()

    one_off = adata.ap.summarize(obs={"cells": ap.n()}, x={"mean_g1": ap.mean("g1")}, by="label")
    persistent = adata.ap.group_by(obs="label").summarize(obs={"cells": ap.n()}, x={"mean_g1": ap.mean("g1")})

    pd.testing.assert_frame_equal(one_off, persistent)


def test_sparse_sources_keep_their_group_alignment() -> None:
    adata = make_adata(sparse.csr_matrix(VALUES))

    summary = adata.ap.summarize(obs={"cells": ap.n()}, x={"mean_g1": ap.mean("g1")}, by="label")

    reference = expected()
    assert summary["label"].tolist() == reference["label"].tolist()
    assert summary["mean_g1"].tolist() == pytest.approx(reference["mean_g1"].tolist())


def test_missing_group_keys_keep_their_aggregates() -> None:
    adata = make_adata()
    adata.obs["label"] = pd.Series(["zeta", None, "mu", "zeta", None, "zeta"], index=adata.obs_names)

    summary = adata.ap.summarize(obs={"cells": ap.n()}, x={"mean_g1": ap.mean("g1")}, by="label")

    assert summary["label"].tolist()[0] == "zeta"
    assert pd.isna(summary["label"].tolist()[1])
    assert summary.loc[summary["label"] == "zeta", "cells"].tolist() == [3]
    assert summary.loc[summary["label"] == "zeta", "mean_g1"].tolist() == pytest.approx([(10 + 40 + 60) / 3])
    assert summary.loc[summary["label"].isna(), "mean_g1"].tolist() == pytest.approx([(20 + 50) / 2])


def test_multiple_group_keys_stay_aligned() -> None:
    adata = make_adata()
    adata.obs["batch"] = ["b2", "b1", "b2", "b2", "b1", "b1"]

    summary = adata.ap.summarize(
        obs={"cells": ap.n()},
        x={"mean_g1": ap.mean("g1")},
        by=["label", "batch"],
    )

    frame = pd.DataFrame({"label": LABELS, "batch": adata.obs["batch"].tolist(), "g1": VALUES[:, 0]})
    reference = frame.groupby(["label", "batch"], sort=False, observed=True)["g1"].agg(["size", "mean"])
    for _, row in summary.iterrows():
        size, mean = reference.loc[(row["label"], row["batch"])]
        assert row["cells"] == size
        assert row["mean_g1"] == pytest.approx(mean)
