"""Wide metadata evaluation and join benchmarks."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

import annplyr as ap

from .common import SEED, metadata_frame


class WideMetadata:
    """200,000 observations with 200 float32 metadata columns."""

    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 1200

    def setup(self) -> None:
        rng = np.random.default_rng(SEED)
        n_obs = 200_000
        obs = metadata_frame(n_obs, 200, rng=rng)
        obs.insert(0, "id", np.arange(n_obs, dtype=np.int64))
        obs.index = pd.Index([f"o{i}" for i in range(n_obs)])
        self.adata = AnnData(X=np.zeros((n_obs, 1), dtype=np.float32), obs=obs, var=pd.DataFrame(index=["g0"]))
        self.selected = [f"m{i}" for i in range(20)]
        right_ids = np.arange(0, n_obs, 2, dtype=np.int64)
        self.right = pd.DataFrame(
            {"id": right_ids, "right_value": rng.standard_normal(len(right_ids), dtype=np.float32)}
        )
        self.full_right = pd.DataFrame(
            {"id": np.arange(n_obs, dtype=np.int64), "right_value": rng.standard_normal(n_obs, dtype=np.float32)}
        )
        self.assignments = {f"derived_{i}": ap.col(f"m{i}") * 2 + 1 for i in range(10)}

    def time_select_twenty(self) -> None:
        out = self.adata.ap.select(obs=["id", *self.selected], copy=True)
        assert out.obs.shape == (200_000, 21)

    def time_mutate_ten_independent(self) -> None:
        out = self.adata.ap.mutate(obs=self.assignments)
        assert all(name in out.obs for name in self.assignments)

    def time_left_join(self) -> None:
        out = self.adata.ap.left_join(self.full_right, by="id", copy=True)
        assert out.n_obs == 200_000
        assert "right_value" in out.obs

    def time_inner_join(self) -> None:
        out = self.adata.ap.inner_join(self.right, by="id", copy=True)
        assert out.n_obs == 100_000

    def time_semi_join(self) -> None:
        out = self.adata.ap.semi_join(self.right, by="id", copy=True)
        assert out.n_obs == 100_000

    def time_anti_join(self) -> None:
        out = self.adata.ap.anti_join(self.right, by="id", copy=True)
        assert out.n_obs == 100_000
