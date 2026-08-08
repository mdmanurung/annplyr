"""Persistent grouping benchmark families."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

import annplyr as ap

from .common import SEED, metadata_frame, realize_frame, unwrap_grouped


def _grouped_fixture(n_groups: int) -> AnnData:
    rng = np.random.default_rng(SEED)
    n_obs = 200_000
    obs = metadata_frame(n_obs, 50, rng=rng)
    obs.insert(0, "group", pd.Categorical(np.arange(n_obs) % n_groups))
    obs.insert(1, "score", rng.standard_normal(n_obs, dtype=np.float32))
    obs.index = pd.Index([f"o{i}" for i in range(n_obs)])
    return AnnData(X=np.zeros((n_obs, 1), dtype=np.float32), obs=obs, var=pd.DataFrame(index=["g0"]))


class GroupedTwenty:
    """200,000 observations split across 20 interleaved groups."""

    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 900

    def setup(self) -> None:
        self.adata = _grouped_fixture(20)
        self.grouped = self.adata.ap.group_by(obs="group")
        frame = self.adata.obs
        means = frame.groupby("group", observed=True)["score"].transform("mean")
        self.expected_filtered = int((frame["score"] > means).sum())

    def time_filter(self) -> None:
        out = self.grouped.filter(obs=ap.col("score") > ap.mean("score"), copy=True)
        assert unwrap_grouped(out).n_obs == self.expected_filtered

    def time_mutate(self) -> None:
        out = self.grouped.mutate(obs={"centered": ap.col("score") - ap.mean("score")})
        ungrouped = unwrap_grouped(out)
        assert ungrouped.n_obs == 200_000
        assert "centered" in ungrouped.obs

    def time_summarize(self) -> None:
        out = self.grouped.summarize(obs={"mean_score": ap.mean("score")})
        assert out.shape == (20, 2)
        realize_frame(out)

    def time_arrange(self) -> None:
        out = self.grouped.arrange(obs=ap.col("score"), copy=True)
        assert unwrap_grouped(out).n_obs == 200_000

    def time_slice_head(self) -> None:
        out = self.grouped.slice_head(n=5, copy=True)
        assert unwrap_grouped(out).n_obs == 100


class GroupedHighCardinality:
    """200,000 observations split across 50,000 observed groups."""

    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 900

    def setup(self) -> None:
        self.adata = _grouped_fixture(50_000)
        self.grouped = self.adata.ap.group_by(obs="group")

    def time_group_keys(self) -> None:
        out = self.grouped.group_keys()
        assert out.shape == (50_000, 1)

    def time_count(self) -> None:
        out = self.grouped.count()
        assert out.shape == (50_000, 2)
        assert int(out["n"].sum()) == 200_000
