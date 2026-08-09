"""Shared PBMC3K fixture for the executable documentation examples.

Every guide, the quickstart, and the cohort vignette run against the same real
dataset: the 2,638-cell PBMC3K object distributed by Scanpy, already annotated
with Louvain cell types, QC metrics, PCA/UMAP/t-SNE embeddings, and a `.raw`
attribute holding log-normalised counts for all 13,714 detected genes.

The dataset is downloaded once into a gitignored ``data/`` directory and cached
in memory, so a full documentation build fetches it at most one time.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import scanpy as sc
from anndata import AnnData

#: Download location. ``data/`` is gitignored at the repository root.
DATASET_DIR = Path(__file__).resolve().parents[1] / "data"


@lru_cache(maxsize=1)
def _source() -> AnnData:
    sc.settings.datasetdir = DATASET_DIR
    return sc.datasets.pbmc3k_processed()


def pbmc3k() -> AnnData:
    """Return an independent copy of the processed PBMC3K dataset.

    Each documentation page gets its own copy, so an example that mutates
    metadata cannot leak into another page.
    """
    return _source().copy()
