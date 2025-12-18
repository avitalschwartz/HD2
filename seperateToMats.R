#!/usr/bin/env Rscript

library(Matrix)
library(readr)
library(dplyr)

base_dir <- "/ems/elsc-labs/meshorer-e/avital.schwartz1/HD/HD2"

# ---- inputs ----
mtx_filt_path <- file.path(base_dir, "GSE152058_human_snRNA_processed_counts_filtered.mtx")

coldata_path  <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")   # original (read-only)
rowdata_path  <- file.path(base_dir, "GSE152058_human_snRNA_processed_rowdata.tsv")   # original (read-only)

seurat_path   <- file.path(base_dir, "HD1_mapped_to_ref.rds")

# ---- output ----
out_dir <- file.path(base_dir, "split_celltypes_counts_rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load filtered matrix ----
mat_filt <- readMM(mtx_filt_path)
mat_filt <- as(mat_filt, "dgCMatrix")

# ---- load original TSVs (read-only) ----
colData <- read_tsv(coldata_path, show_col_types = FALSE)
rowData <- read_tsv(rowdata_path, show_col_types = FALSE)

# ---- re-create keep_idx EXACTLY as your filtering rule (>=1000 per NBB_ID×Region) ----
threshold <- 1000
group_sizes <- colData %>% count(NBB_ID, Region, name = "n_cells")
drop_groups <- group_sizes %>% filter(n_cells < threshold)

drop_idx <- colData %>%
  mutate(.idx = row_number()) %>%
  semi_join(drop_groups, by = c("NBB_ID", "Region")) %>%
  pull(.idx)

keep_idx <- setdiff(seq_len(nrow(colData)), drop_idx)

# Now keep_idx should correspond to the columns kept in the filtered MTX:
stopifnot(length(keep_idx) == ncol(mat_filt))

# ---- attach names (needed for matching to Seurat) ----
stopifnot("Barcode" %in% names(colData))
stopifnot("Gene" %in% names(rowData))

colnames(mat_filt) <- colData$Barcode[keep_idx]
rownames(mat_filt) <- rowData$Gene

# ---- load Seurat object and predicted cell types ----
obj <- readRDS(seurat_path)
meta <- obj@meta.data

celltype_col <- "predicted_celltype"
stopifnot(celltype_col %in% names(meta))

# match barcodes in filtered matrix to Seurat meta rows
idx <- match(colnames(mat_filt), rownames(meta))

if (anyNA(idx)) {
  warning(sum(is.na(idx)), " filtered cells not found in Seurat object; dropping them.")
  keep_cols <- which(!is.na(idx))
  mat_filt <- mat_filt[, keep_cols, drop = FALSE]
  idx <- idx[keep_cols]
}

celltypes <- meta[[celltype_col]][idx]

# ---- split + save per cell type (one-by-one, not all in RAM) ----
cts <- sort(unique(na.omit(celltypes)))

for (ct in cts) {
  cols <- which(celltypes == ct)
  if (!length(cols)) next

  m_ct <- mat_filt[, cols, drop = FALSE]
  saveRDS(m_ct, file.path(out_dir, paste0("counts_", ct, ".rds")))

  rm(m_ct); gc()
  cat("Saved", ct, ":", length(cols), "cells\n")
}

cat("\nDone. Output:", out_dir, "\n")
