from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

FIXTURE_VERSION = 1


def make_representative_adata() -> AnnData:
    """Build the deterministic, redistributable annplyr integration fixture."""
    counts = np.array(
        [
            [8, 0, 1, 0, 3, 0],
            [0, 5, 0, 2, 0, 1],
            [4, 0, 3, 0, 1, 0],
            [0, 2, 0, 6, 0, 2],
            [7, 1, 0, 0, 2, 0],
            [3, 0, 5, 1, 0, 0],
            [0, 4, 0, 3, 0, 5],
            [6, 0, 2, 0, 4, 1],
        ],
        dtype=np.int32,
    )
    obs = pd.DataFrame(
        {
            "fixture_row_id": np.arange(8, dtype=np.int16),
            "cell_type": pd.Categorical(
                ["T cell", "B cell", "T cell", "Monocyte", "B cell", "T cell", "Monocyte", "B cell"],
                categories=["B cell", "T cell", "Monocyte"],
                ordered=True,
            ),
            "batch": pd.Categorical(["run-1"] * 4 + ["run-2"] * 4),
            "detected_genes": pd.array([3, 3, 3, pd.NA, 3, 3, 4, 4], dtype="Int64"),
            "qc_pass": pd.array([True, True, pd.NA, False, True, True, False, True], dtype="boolean"),
            "quality_score": np.array([8, 6, 5, 3, 9, 7, 10, 4], dtype=np.float32),
        },
        index=pd.Index(["cell-a", "cell-b", "cell-a", "cell-d", "cell-b", "cell-f", "cell-d", "cell-b"]),
    )
    var = pd.DataFrame(
        {
            "fixture_col_id": np.arange(6, dtype=np.int16),
            "feature_type": pd.Categorical(
                ["gene", "gene", "protein", "gene", "protein", "gene"],
                categories=["gene", "protein"],
            ),
            "highly_variable": pd.array([True, False, True, pd.NA, False, True], dtype="boolean"),
            "feature_rank": pd.array([1, 2, pd.NA, 4, 5, 6], dtype="Int64"),
        },
        index=pd.Index(["GATA3", "MS4A1", "CD3D", "LYZ", "CD79A", "NKG7"]),
    )
    adata = AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    adata.layers["log1p"] = sparse.csr_matrix(np.log1p(counts).astype(np.float32))
    adata.obsm["X_reference"] = np.array(
        [[-1.0, 0.2], [-0.7, -0.1], [-0.4, 0.4], [0.0, -0.5], [0.3, 0.0], [0.6, 0.5], [0.8, -0.3], [1.0, 0.1]],
        dtype=np.float32,
    )
    adata.varm["reference_loadings"] = np.array(
        [[0.4, -0.1], [-0.3, 0.2], [0.5, 0.6], [-0.2, 0.7], [0.1, -0.4], [0.8, 0.3]],
        dtype=np.float32,
    )
    first_block = sparse.csr_matrix(
        np.array(
            [[1.0, 0.5, 0.2, 0.0], [0.5, 1.0, 0.0, 0.3], [0.2, 0.0, 1.0, 0.4], [0.0, 0.3, 0.4, 1.0]],
            dtype=np.float32,
        )
    )
    second_block = sparse.csr_matrix(
        np.array(
            [[1.0, 0.1, 0.0, 0.6], [0.1, 1.0, 0.7, 0.0], [0.0, 0.7, 1.0, 0.2], [0.6, 0.0, 0.2, 1.0]],
            dtype=np.float32,
        )
    )
    adata.obsp["connectivities"] = sparse.block_diag((first_block, second_block), format="csr")
    adata.varp["feature_similarity"] = sparse.csc_matrix(
        np.array(
            [
                [1.0, 0.2, 0.4, 0.0, 0.1, 0.5],
                [0.2, 1.0, 0.0, 0.3, 0.2, 0.0],
                [0.4, 0.0, 1.0, 0.1, 0.6, 0.3],
                [0.0, 0.3, 0.1, 1.0, 0.0, 0.2],
                [0.1, 0.2, 0.6, 0.0, 1.0, 0.4],
                [0.5, 0.0, 0.3, 0.2, 0.4, 1.0],
            ],
            dtype=np.float32,
        )
    )
    adata.uns["fixture"] = {
        "version": FIXTURE_VERSION,
        "origin": "synthetic explicit counts",
        "license": "BSD-3-Clause",
    }
    adata.raw = AnnData(
        X=adata.X.copy(),
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )
    return adata


def write_fixture(output: Path) -> tuple[Path, Path]:
    """Regenerate both serialized fixture forms below ``output``."""
    output.mkdir(parents=True, exist_ok=True)
    adata = make_representative_adata()
    h5ad_path = output / "representative.h5ad"
    zarr_path = output / "representative.zarr"
    adata.write_h5ad(h5ad_path)
    adata.write_zarr(zarr_path)
    return h5ad_path, zarr_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate the annplyr integration fixture")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    for path in write_fixture(args.output):
        print(path)


if __name__ == "__main__":
    main()
