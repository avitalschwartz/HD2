#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(readr)
  library(dplyr)
})

# =========================
# Paths
# =========================
base_dir <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2"

in_dir   <- file.path(base_dir, "split_celltypes_pseudobulk_rds_filterByExpr")
out_dir  <- file.path(base_dir, "split_celltypes_pseudobulk_edgeR_TMM")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# (optional) metadata for nicer per-sample tables
coldata_tsv <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")
colData <- read_tsv(coldata_tsv, show_col_types = FALSE) %>%
  transmute(
    NBB_ID = as.character(NBB_ID),
    Region = as.character(Region),
    Condition = as.character(Condition),
    sample_id = paste(NBB_ID, Region, sep="__")
  ) %>% distinct()

# =========================
# Process files
# =========================
files <- list.files(in_dir, pattern = "_filterByExpr\\.rds$", full.names = TRUE)

summary <- data.frame(
  file = character(),
  genes = integer(),
  samples = integer(),
  stringsAsFactors = FALSE
)

for (f in files) {
  cat("\n=== Processing:", basename(f), "===\n")
  
  pb <- readRDS(f)
  if (!inherits(pb, "dgCMatrix")) pb <- as(pb, "dgCMatrix")
  
  counts <- as.matrix(pb)  # edgeR expects a regular matrix
  stopifnot(!is.null(colnames(counts)))
  
  # Build DGEList
  y <- DGEList(counts = counts)
  
  # TMM normalization (writes norm.factors into y$samples)
  y <- calcNormFactors(y, method = "TMM")
  
  # Save the DGEList (recommended for downstream edgeR pipeline)
  out_rds <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_TMM_DGEList.rds")
  )
  saveRDS(y, out_rds)
  cat("Saved DGEList:", out_rds, "\n")
  
  # Save a small CSV of factors + effective library sizes
  nf <- y$samples$norm.factors
  lib <- y$samples$lib.size
  eff_lib <- lib * nf
  
  samp <- data.frame(
    sample_id = rownames(y$samples),
    lib_size = lib,
    norm_factor = nf,
    effective_lib_size = eff_lib,
    stringsAsFactors = FALSE
  )
  
  # Attach metadata if available
  samp <- samp %>%
    left_join(colData, by = "sample_id")
  
  out_csv <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_TMM_factors.csv")
  )
  write.csv(samp, out_csv, row.names = FALSE)
  cat("Saved factors:", out_csv, "\n")
  
  # OPTIONAL: save logCPM matrix (often used for PCA/plots)
  # (This is not used for edgeR DE testing directly; it’s for exploration.)
  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  out_logcpm <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_logCPM_TMM.csv")
  )
  write.csv(logcpm, out_logcpm)
  cat("Saved logCPM:", out_logcpm, "\n")
  
  summary <- rbind(summary, data.frame(
    file = basename(f),
    genes = nrow(counts),
    samples = ncol(counts),
    stringsAsFactors = FALSE
  ))
  
  rm(pb, counts, y, samp, logcpm); gc()
}

write.csv(summary, file.path(out_dir, "TMM_summary.csv"), row.names = FALSE)
cat("\nDone. Output folder:", out_dir, "\n")
