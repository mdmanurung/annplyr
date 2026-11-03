# Purpose: run ten tidySCE vignette tasks with annplyr on PBMC3k and write each result.
# Run: uvx hatch -e integration run python scripts/r_parity/01_annplyr.py
# Outputs: out/X.mtx, out/obs.csv, out/var.csv (shared input for R), out/py/*.csv
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

import annplyr as ap
from annplyr.data import X, obs

OUT = Path(__file__).resolve().parent / "out"  # gitignored; rerun to refresh
MYELOID = ["CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic cells"]

(OUT / "py").mkdir(parents=True, exist_ok=True)
(OUT / "r").mkdir(exist_ok=True)

# Shared input ----------------------------------------------------------------
# R reads this file, so both sides see the same log-normalized matrix as X.
pbmc = sc.datasets.pbmc3k_processed()
cells = ad.AnnData(X=pbmc.raw.X.copy(), obs=pbmc.obs.copy(), var=pbmc.raw.var.copy())
cells.obs["donor"] = np.array(["d1", "d2", "d3"])[np.arange(cells.n_obs) % 3]  # synthetic donors
cells.obs["louvain"] = cells.obs["louvain"].astype(str)
# anndataR's HDF5 dependency does not build here, so R reads Matrix Market + CSV
scipy.io.mmwrite(OUT / "X.mtx", cells.X.T.tocsr().astype(np.float64))  # genes x cells, as in SCE
cells.obs.rename_axis("cell").to_csv(OUT / "obs.csv")
pd.DataFrame(index=cells.var_names.rename("gene")).to_csv(OUT / "var.csv")
adata = ap.tidy(cells)

# Tasks -----------------------------------------------------------------------
kept = adata.filter(obs.percent_mito < 0.05, obs.n_genes > 300)
pd.DataFrame({"cell": kept.obs_names}).to_csv(OUT / "py/01_filter.csv", index=False)

adata.count(obs.donor, obs.louvain).to_csv(OUT / "py/02_count.csv", index=False)

scored = adata.mutate(b_score=(X.MS4A1 + X.CD79A) / 2)
scored.obs[["b_score"]].rename_axis("cell").to_csv(OUT / "py/03_mutate.csv")

adata.group_by(obs.louvain).summarize(n=ap.n(), median_genes=ap.median(obs.n_genes), mean_LYZ=ap.mean(X.LYZ)).to_csv(
    OUT / "py/04_summarise.csv", index=False
)

top = adata.slice_max(obs.n_counts, n=3, by=obs.louvain)
pd.DataFrame({"cell": top.obs_names}).to_csv(OUT / "py/05_slice_max.csv", index=False)

long = adata.join_features(X.MS4A1, X.LYZ, X["HLA-DRA"], shape="long")
long[[".cell", ".feature", ".abundance_X"]].assign(
    value=lambda d: d[".abundance_X"].sparse.to_dense() if hasattr(d[".abundance_X"], "sparse") else d[".abundance_X"]
).drop(columns=".abundance_X").to_csv(OUT / "py/06_join_features.csv", index=False)

classed = adata.rename(cell_type=obs.louvain).mutate(
    cell_class=ap.if_else(obs.cell_type.is_in(MYELOID), "myeloid", "lymphoid")
)
classed.count(obs.cell_class).to_csv(OUT / "py/07_if_else.csv", index=False)

sheet = pd.DataFrame({"donor": ["d1", "d2", "d3"], "age": [31, 47, 62]})
joined = adata.left_join(sheet, by=obs.donor)
joined.obs[["age"]].rename_axis("cell").to_csv(OUT / "py/08_left_join.csv")

pb = adata.aggregate_cells(obs.donor, obs.louvain, source=X)
pb_frame = pd.DataFrame(pb.X.toarray(), index=pb.obs_names, columns=pb.var_names)[["MS4A1", "LYZ", "CD3E", "NKG7"]]
pb_frame.assign(n_cells=pb.obs["n_cells"].to_numpy()).rename_axis("group").to_csv(OUT / "py/09_aggregate.csv")

ordered = adata.arrange(ap.desc(obs.n_genes), obs.n_counts)
pd.DataFrame({"cell": ordered.obs_names[:50]}).to_csv(OUT / "py/10_arrange.csv", index=False)
print("annplyr results written")
