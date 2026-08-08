"""CI-sized benchmark execution and correctness smoke cases."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

import annplyr as ap

from .common import SEED, dense_adata, feature_names, sparse_adata, unwrap_grouped


class Smoke:
    """Frozen CI sizes: 1,000 dense, 5,000 sparse, and 5,000 grouped rows."""

    number = 1
    repeat = 1
    timeout = 300

    def setup(self) -> None:
        self.dense = dense_adata(1_000, 200)
        self.sparse = sparse_adata(5_000, 1_000, matrix_format="csr")
        rng = np.random.default_rng(SEED)
        obs = pd.DataFrame(
            {
                "group": pd.Categorical(np.arange(5_000) % 20),
                "score": rng.standard_normal(5_000, dtype=np.float32),
            },
            index=pd.Index([f"o{i}" for i in range(5_000)]),
        )
        grouped_adata = AnnData(X=np.zeros((5_000, 1), dtype=np.float32), obs=obs, var=pd.DataFrame(index=["g0"]))
        self.grouped = grouped_adata.ap.group_by(obs="group")
        self.features = feature_names(10)

    def time_dense_projection(self) -> None:
        out = self.dense.ap.to_df(x=self.features)
        assert out.shape == (1_000, 10)

    def time_sparse_projection(self) -> None:
        out = self.sparse.ap.to_df(x=self.features)
        assert out.shape == (5_000, 10)

    def time_grouped_filter(self) -> None:
        out = self.grouped.filter(obs=ap.col("score") > ap.mean("score"), copy=True)
        assert 0 < unwrap_grouped(out).n_obs < 5_000

    def time_grouped_mutate(self) -> None:
        out = self.grouped.mutate(obs={"centered": ap.col("score") - ap.mean("score")})
        assert "centered" in unwrap_grouped(out).obs

    def time_grouped_summarize(self) -> None:
        out = self.grouped.summarize(obs={"mean_score": ap.mean("score")})
        assert out.shape == (20, 2)
