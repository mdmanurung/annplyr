"""Dense, sparse, and backed projection benchmarks."""

from __future__ import annotations

import tempfile

import numpy as np
from anndata import read_h5ad
from scipy import sparse

import annplyr as ap

from .common import assert_columns, dense_adata, feature_names, realize_frame, sparse_adata


class DenseProjection:
    """20,000 x 2,000 dense projection scenario."""

    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 600

    def setup(self) -> None:
        self.adata = dense_adata(20_000, 2_000)
        self.features = feature_names(10)
        values = np.asarray(self.adata.X)[:, :3]
        self.expected_filtered = int(np.any(values > 0, axis=1).sum())

    def time_filter_three_features(self) -> None:
        predicate = (ap.col("g0") > 0) | (ap.col("g1") > 0) | (ap.col("g2") > 0)
        out = self.adata.ap.filter(x=predicate, copy=True)
        assert out.n_obs == self.expected_filtered
        assert out.n_vars == self.adata.n_vars

    def time_export_ten_features(self) -> None:
        out = self.adata.ap.to_df(x=self.features)
        assert out.shape == (20_000, 10)
        assert_columns(out.columns, self.features)
        realize_frame(out)

    def time_copy_false_call_ten_features(self) -> None:
        out = self.adata.ap.select(x=self.features, copy=False)
        assert out.shape == (20_000, 10)
        assert out.n_vars == 10
        assert_columns(out.var_names, self.features)
        assert self.adata.shape == (20_000, 2_000)

    def time_copy_false_realization_ten_features(self) -> None:
        out = self.adata.ap.select(x=self.features, copy=False)
        values = np.asarray(out.X)
        assert values.shape == (20_000, 10)
        assert_columns(out.var_names, self.features)
        assert self.adata.shape == (20_000, 2_000)
        float(values[0].sum())


class SparseProjection:
    """100,000 x 10,000 CSR/CSC projection and export scenario."""

    params = ["csr", "csc"]
    param_names = ["format"]
    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 900

    def setup(self, matrix_format: str) -> None:
        self.adata = sparse_adata(100_000, 10_000, matrix_format=matrix_format)
        self.features = feature_names(10)
        positive_counts = (self.adata.X[:, :3] > 0).getnnz(axis=1)
        self.expected_filtered = int((positive_counts > 0).sum())

    def time_filter_three_features(self, matrix_format: str) -> None:
        predicate = (ap.col("g0") > 0) | (ap.col("g1") > 0) | (ap.col("g2") > 0)
        out = self.adata.ap.filter(x=predicate, copy=True)
        assert out.n_obs == self.expected_filtered
        assert out.X.format == matrix_format

    def time_pull_one_feature(self, matrix_format: str) -> None:
        out = self.adata.ap.pull(x="g0")
        assert len(out) == 100_000
        float(out.iloc[0])

    def time_to_df_ten_features(self, matrix_format: str) -> None:
        out = self.adata.ap.to_df(x=self.features)
        assert out.shape == (100_000, 10)
        assert_columns(out.columns, self.features)
        realize_frame(out)

    def time_to_tidy_ten_features(self, matrix_format: str) -> None:
        out = self.adata.ap.to_tidy(x=self.features)
        assert out.shape[0] == 1_000_000
        assert set(out["feature"].unique()) == set(self.features)
        realize_frame(out)

    def time_budget_rejection(self, matrix_format: str) -> None:
        try:
            self.adata.ap.to_df(x=self.features, max_matrix_values=1)
        except ap.AnnplyrError:
            return
        raise AssertionError("projected matrix budget was not rejected")


class BackedProjection:
    """Warm-cache selected reads from backed dense, CSR, and CSC H5AD."""

    params = ["dense", "csr", "csc"]
    param_names = ["storage"]
    number = 1
    repeat = 7
    warmup_time = 0.1
    timeout = 1200

    def setup(self, storage: str) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="annplyr-asv-backed-")
        self.path = f"{self._tmp.name}/{storage}.h5ad"
        if storage == "dense":
            source = dense_adata(20_000, 2_000)
        else:
            source = sparse_adata(100_000, 10_000, matrix_format=storage)
        projected = source.X[:, :3]
        if sparse.issparse(projected):
            positive_counts = (projected > 0).getnnz(axis=1)
            mask = positive_counts > 0
        else:
            mask = np.any(np.asarray(projected) > 0, axis=1)
        self.expected_filtered = int(mask.sum())
        source.write_h5ad(self.path)
        self.adata = read_h5ad(self.path, backed="r")
        self.features = feature_names(10)
        self._legacy_sparse = ap.__version__.startswith("0.2.") and storage != "dense"

    def _target(self):
        """Return the timed v0.2 workaround for unsupported backed sparse X."""
        if self._legacy_sparse:
            return self.adata.to_memory()
        return self.adata

    def teardown(self, storage: str) -> None:
        self.adata.file.close()
        self._tmp.cleanup()

    def time_selected_read(self, storage: str) -> None:
        target = self._target()
        out = target.ap.to_df(x=self.features)
        expected_rows = 20_000 if storage == "dense" else 100_000
        assert out.shape == (expected_rows, 10)
        assert_columns(out.columns, self.features)
        realize_frame(out)

    def time_filter_three_features(self, storage: str) -> None:
        target = self._target()
        predicate = (ap.col("g0") > 0) | (ap.col("g1") > 0) | (ap.col("g2") > 0)
        out = target.ap.filter(x=predicate, copy=False)
        assert out.n_obs == self.expected_filtered
        assert out.n_vars == target.n_vars

    def time_export_ten_features(self, storage: str) -> None:
        target = self._target()
        out = target.ap.to_tidy(x=self.features)
        expected_rows = (20_000 if storage == "dense" else 100_000) * 10
        assert out.shape[0] == expected_rows
        realize_frame(out)

    def time_budget_rejection(self, storage: str) -> None:
        target = self._target()
        try:
            target.ap.to_df(x=self.features, max_matrix_values=1)
        except ap.AnnplyrError:
            return
        raise AssertionError("projected matrix budget was not rejected")
