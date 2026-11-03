# Purpose: run the same ten tasks with tidySingleCellExperiment and write each result.
# Run: cd scripts/r_parity/out && Rscript ../02_tidysce.R
#      (R 4.5 with tidySingleCellExperiment from Bioconductor 3.22)
# Inputs: X.mtx, obs.csv, var.csv [01_annplyr.py]
# Outputs: out/r/*.csv, out/r/session_info.txt
suppressPackageStartupMessages({
  library(tidySingleCellExperiment)
  library(dplyr)
  library(tidyr)
})

MYELOID <- c("CD14+ Monocytes", "FCGR3A+ Monocytes", "Dendritic cells")

# Input ------------------------------------------------------------------------
# One assay, X: the log-normalized raw matrix that annplyr reads as X.
obs <- read.csv("obs.csv", check.names = FALSE)
genes <- read.csv("var.csv")$gene
x <- as(Matrix::readMM("X.mtx"), "CsparseMatrix")
dimnames(x) <- list(genes, obs$cell)
stopifnot(identical(colnames(x), obs$cell))
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(X = x),
  colData = S4Vectors::DataFrame(obs[, -1], row.names = obs$cell, check.names = FALSE)
)
cat("assays:", SummarizedExperiment::assayNames(sce), "\n")

# Tasks ------------------------------------------------------------------------
kept <- sce |> filter(percent_mito < 0.05, n_genes > 300)
write.csv(data.frame(cell = colnames(kept)), "r/01_filter.csv", row.names = FALSE)

sce |> count(donor, louvain) |> write.csv("r/02_count.csv", row.names = FALSE)

scored <- sce |>
  join_features(c("MS4A1", "CD79A"), shape = "wide") |>
  mutate(b_score = (MS4A1 + CD79A) / 2)
write.csv(data.frame(cell = colnames(scored), b_score = scored$b_score), "r/03_mutate.csv", row.names = FALSE)

sce |>
  join_features("LYZ", shape = "wide") |>
  group_by(louvain) |>
  summarise(n = n(), median_genes = median(n_genes), mean_LYZ = mean(LYZ)) |>
  write.csv("r/04_summarise.csv", row.names = FALSE)

top <- sce |> slice_max(n_counts, n = 3, by = louvain, with_ties = FALSE)
write.csv(data.frame(cell = colnames(top)), "r/05_slice_max.csv", row.names = FALSE)

sce |>
  join_features(c("MS4A1", "LYZ", "HLA-DRA"), shape = "long") |>
  select(.cell, .feature, value = .abundance_X) |>
  write.csv("r/06_join_features.csv", row.names = FALSE)

sce |>
  rename(cell_type = louvain) |>
  mutate(cell_class = if_else(cell_type %in% MYELOID, "myeloid", "lymphoid")) |>
  count(cell_class) |>
  write.csv("r/07_if_else.csv", row.names = FALSE)

sheet <- tibble(donor = c("d1", "d2", "d3"), age = c(31, 47, 62))
joined <- sce |> left_join(sheet, by = "donor")
write.csv(data.frame(cell = colnames(joined), age = joined$age), "r/08_left_join.csv", row.names = FALSE)

pb <- sce |> aggregate_cells(.sample = c(donor, louvain), assays = "X")
pb_x <- as.matrix(SummarizedExperiment::assay(pb, "X"))[c("MS4A1", "LYZ", "CD3E", "NKG7"), ]
pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
write.csv(
  data.frame(group = paste(pb_meta$donor, pb_meta$louvain, sep = "_"), t(pb_x), n_cells = pb_meta$.aggregated_cells, check.names = FALSE),
  "r/09_aggregate.csv", row.names = FALSE
)

ordered <- sce |> arrange(desc(n_genes), n_counts)
write.csv(data.frame(cell = colnames(ordered)[1:50]), "r/10_arrange.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "r/session_info.txt")
cat("tidySCE results written\n")
