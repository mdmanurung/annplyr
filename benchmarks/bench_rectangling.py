"""Rectangling and chained extraction benchmarks."""

from __future__ import annotations

import numpy as np
import pandas as pd

import annplyr as ap

from .common import SEED, dense_adata, feature_names, realize_frame


class Rectangling:
    """100,000 table rows and a 20,000 x 50 selected matrix export."""

    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 900

    def setup(self) -> None:
        rng = np.random.default_rng(SEED)
        n_rows = 100_000
        self.long = pd.DataFrame(
            {
                "id": np.repeat(np.arange(n_rows // 2), 2),
                "key": np.tile(["a", "b"], n_rows // 2),
                "value": rng.standard_normal(n_rows, dtype=np.float32),
            }
        )
        self.nested = pd.DataFrame(
            {
                "group": np.arange(n_rows // 4),
                "values": [[i, i + 1, i + 2, i + 3] for i in range(n_rows // 4)],
            }
        )
        self.adata = dense_adata(20_000, 50, metadata_columns=2)
        self.features = feature_names(50)

    def time_pivot_wider(self) -> None:
        out = ap.pivot_wider(self.long, names_from="key", values_from="value", id_cols="id")
        assert out.shape == (50_000, 3)

    def time_nest_unnest(self) -> None:
        out = ap.unnest_longer(self.nested, "values")
        assert out.shape[0] == 100_000

    def time_chained_extraction(self) -> None:
        tidy = self.adata.ap.pivot_longer(obs=["m0", "m1"], x=self.features)
        out = ap.pivot_wider(tidy, names_from="name", values_from="value", id_cols=["obs_name", "m0", "m1"])
        assert out.shape == (20_000, 53)
        realize_frame(out)
