from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse


@pytest.fixture
def dense_adata() -> AnnData:
    x = np.array(
        [
            [1, 0, 5, 2],
            [2, 3, 0, 1],
            [0, 1, 4, 0],
            [5, 2, 3, 8],
            [3, 0, 1, 4],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(
        {
            "batch": ["A", "A", "B", "B", "A"],
            "score": [1.0, 2.0, 0.5, 3.0, 2.5],
            "cell_type": ["T", "B", "T", "Mono", "B"],
        },
        index=pd.Index([f"c{i}" for i in range(5)], name="cell"),
    )
    var = pd.DataFrame(
        {
            "feature_type": ["rna", "rna", "protein", "rna"],
            "chrom": ["chr1", "chr2", "chr1", "chr3"],
            "length": [100, 200, 150, 80],
        },
        index=pd.Index([f"g{i}" for i in range(4)], name="gene"),
    )
    adata = AnnData(X=x, obs=obs, var=var)
    adata.layers["counts"] = x * 10
    adata.obsm["X_pca"] = np.array(
        [
            [0.1, 1.0],
            [-0.2, 0.0],
            [0.3, -1.0],
            [1.0, 2.0],
            [0.5, 1.5],
        ],
        dtype=float,
    )
    adata.varm["loadings"] = np.array(
        [
            [0.5, 2.0],
            [-0.1, 3.0],
            [0.3, 1.0],
            [1.2, 0.0],
        ],
        dtype=float,
    )
    return adata


@pytest.fixture
def sparse_adata(dense_adata: AnnData) -> AnnData:
    adata = dense_adata.copy()
    adata.X = sparse.csr_matrix(adata.X)
    adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])
    return adata
