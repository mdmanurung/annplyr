from __future__ import annotations

import pandas as pd
from anndata import AnnData

from annplyr._nest import nest_adata, unnest_adata


def test_nest_then_unnest_round_trips_every_cell(dense_adata: AnnData) -> None:
    nested = nest_adata(dense_adata, "batch")
    assert nested["batch"].tolist() == list(dict.fromkeys(dense_adata.obs["batch"]))
    for batch, sub in zip(nested["batch"], nested["data"], strict=True):
        assert (sub.obs["batch"] == batch).all()
        sub.obs["seen"] = True  # a per-group step

    back = unnest_adata(nested)
    assert sorted(back.obs_names) == sorted(dense_adata.obs_names)
    back = back[dense_adata.obs_names]
    pd.testing.assert_frame_equal(back.obs.drop(columns="seen"), dense_adata.obs, check_dtype=False)
    assert back.obs["seen"].all()
    assert (back.X == dense_adata.X).all()
