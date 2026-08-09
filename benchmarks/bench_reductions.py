"""Bounded dense matrix-reduction benchmarks."""

from __future__ import annotations

import hashlib

import numpy as np

import annplyr as ap

from .common import dense_adata


class DenseReductions:
    """100,000 x 500 dense means with exact result hashing."""

    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 900

    def setup(self) -> None:
        self.adata = dense_adata(100_000, 500)
        self.assignments = {f"mean_g{position}": ap.mean(f"g{position}") for position in range(500)}
        values = np.asarray(self.adata.X)
        expected = np.asarray([values[:, position].mean() for position in range(values.shape[1])])
        self.expected_hash = hashlib.sha256(expected.tobytes()).hexdigest()
        self.result_hash = ""

    def time_mean_all_features(self) -> None:
        out = self.adata.ap.summarize(x=self.assignments)
        assert out.shape == (1, 500)
        self.result_hash = hashlib.sha256(out.to_numpy(copy=False).tobytes()).hexdigest()
        assert self.result_hash == self.expected_hash
