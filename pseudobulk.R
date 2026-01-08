#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
})

# =========================
# Paths (EDIT if needed)
# =========================
base_dir   <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2"
input_dir  <- file.path(base_dir, "split_celltypes_fix_filtered")
coldata_tsv <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")

output_dir <- file.path(base_dir, "split_celltypes_pseudobulk_rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load sample metadata
# =========================
colData <- read_tsv(coldata_tsv, show_col_types = FALSE)

required <- c("Barcode", "NBB_ID", "Region")
missing <- setdiff(required, colnames(colData))
if (length(missing) > 0) {
  stop("coldata TSV is missing columns: ", paste(missing, collapse = ", "))
}

# Ensure Barcode is unique (if duplicated, keep first)
colData <- colData %>% distinct(Barcode, .keep_all = TRUE)

# Fast lookup table: Barcode -> (NBB_ID, Region)
# We'll subset this per matrix
colData_small <- colData %>%
  transmute(
    Barcode = as.character(Barcode),
    NBB_ID  = as.character(NBB_ID),
    Region  = as.character(Region)
  )

# =========================
# Pseudobulk function
# mat: genes x cells (dgCMatrix recommended)
# sample_id: length = ncol(mat), e.g. "A47L__Putamen"
# returns: genes x samples summed counts (dgCMatrix)
# =========================
pseudobulk_sum <- function(mat, sample_id) {
  stopifnot(ncol(mat) == length(sample_id))
  
  # keep sample order as it appears
  sample_fac <- factor(sample_id, levels = unique(sample_id))
  S <- nlevels(sample_fac)
  
  # membership matrix: cells x samples
  samp_mat <- Matrix::sparseMatrix(
    i = seq_along(sample_fac),
    j = as.integer(sample_fac),
    x = 1,
    dims = c(length(sample_fac), S)
  )
  
  # genes x samples sums
  pb <- mat %*% samp_mat
  colnames(pb) <- levels(sample_fac)
  pb
}

# =========================
# Loop over matrices
# =========================
files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

summary <- data.frame(
  file = character(),
  n_genes = integer(),
  n_cells_in = integer(),
  n_samples_out = integer(),
  n_cells_dropped_no_meta = integer(),
  stringsAsFactors = FALSE
)

for (f in files) {
  cat("\n=== Processing:", basename(f), "===\n")
  
  mat <- readRDS(f)
  
  # Ensure sparse dgCMatrix
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }
  
  if (is.null(colnames(mat))) {
    stop("Matrix has no colnames (barcodes): ", basename(f))
  }
  
  barcodes <- colnames(mat)
  
  # Match barcodes to colData
  idx <- match(barcodes, colData_small$Barcode)
  
  if (anyNA(idx)) {
    # Drop cells not found in metadata
    keep_cols <- which(!is.na(idx))
    dropped <- sum(is.na(idx))
    
    warning(basename(f), ": dropping ", dropped, " cells not found in coldata.tsv")
    
    mat <- mat[, keep_cols, drop = FALSE]
    idx <- idx[keep_cols]
    barcodes <- colnames(mat)
  } else {
    dropped <- 0
  }
  
  md <- colData_small[idx, , drop = FALSE]
  
  # Build sample IDs: patient__region
  sample_id <- paste(md$NBB_ID, md$Region, sep = "__")
  
  # Pseudobulk sum
  pb <- pseudobulk_sum(mat, sample_id)
  
  # Save output
  out_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_pseudobulk.rds")
  )
  saveRDS(pb, out_file)
  
  cat("Saved pseudobulk:", out_file, "\n")
  cat("Input cells:", ncol(mat), " | Output samples:", ncol(pb), "\n")
  
  summary <- rbind(summary, data.frame(
    file = basename(f),
    n_genes = nrow(mat),
    n_cells_in = ncol(mat),
    n_samples_out = ncol(pb),
    n_cells_dropped_no_meta = dropped,
    stringsAsFactors = FALSE
  ))
  
  rm(mat, pb, md); gc()
}

# Write summary CSV
summary_file <- file.path(output_dir, "pseudobulk_summary.csv")
write.csv(summary, summary_file, row.names = FALSE)

cat("\nDone.\nPseudobulk folder:", output_dir, "\nSummary:", summary_file, "\n")
obj <- readRDS("C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\split_celltypes_pseudobulk_rds\\counts_Astrocyte_filtered_pseudobulk.rds")

View(as.data.frame(as.matrix(obj)))

