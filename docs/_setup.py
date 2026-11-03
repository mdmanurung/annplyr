"""Shared setup for the executable documentation pages.

Each page runs ``from _setup import *`` in a hidden cell. Every page executes in its
own kernel, so ``adata`` is a fresh copy and one page cannot leak into another.
PBMC3k is downloaded once into the gitignored ``data/`` directory at the repo root.
"""

from pathlib import Path

import ggann as ag
import numpy as np
import pandas as pd
import scanpy as sc

import annplyr as ap
from annplyr.data import X, obs, raw, var

sc.settings.datasetdir = Path(__file__).resolve().parents[1] / "data"
pd.set_option("display.width", 88)
pd.set_option("display.max_columns", 12)

adata = ap.tidy(sc.datasets.pbmc3k_processed())

__all__ = ["X", "ag", "adata", "ap", "np", "obs", "pd", "raw", "sc", "var"]
