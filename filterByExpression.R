#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(edgeR)
})
# =========================
# Paths
# =========================
base_dir    <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2"
pb_dir      <- file.path(base_dir, "split_celltypes_pseudobulk_rds")
coldata_tsv <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")

out_dir <- file.path(base_dir, "split_celltypes_pseudobulk_rds_filterByExpr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load sample metadata (coldata)
# =========================
colData <- read_tsv(coldata_tsv, show_col_types = FALSE)

required <- c("Barcode", "NBB_ID", "Region", "Condition")
missing <- setdiff(required, colnames(colData))
if (length(missing) > 0) stop("coldata TSV missing: ", paste(missing, collapse = ", "))

# Build ONE row per sample (patient x region) since pseudobulk is per sample
sample_meta <- colData %>%
  transmute(
    NBB_ID    = as.character(NBB_ID),
    Region    = as.character(Region),
    Condition = as.character(Condition)
  ) %>%
  distinct() %>%
  mutate(sample_id = paste(NBB_ID, Region, sep = "__"))

# If the TSV has repeated (NBB_ID, Region) with conflicting Condition, this will catch it:
dup_check <- sample_meta %>% count(sample_id) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  stop("Duplicate sample_id rows in sample_meta (unexpected). Example: ", dup_check$sample_id[1])
}

# =========================
# Process each pseudobulk matrix
# =========================
files <- list.files(pb_dir, pattern = "_pseudobulk\\.rds$", full.names = TRUE)

summary <- data.frame(
  file = character(),
  genes_before = integer(),
  genes_after = integer(),
  genes_removed = integer(),
  stringsAsFactors = FALSE
)

for (f in files) {
  cat("\n=== Processing:", basename(f), "===\n")
  
  pb <- readRDS(f)
  if (!inherits(pb, "dgCMatrix")) pb <- as(pb, "dgCMatrix")
  
  # edgeR prefers a standard matrix; pseudobulk usually small enough in columns
  counts <- as.matrix(pb)
  
  # Build group/design matching the columns
  cols <- colnames(counts)
  md <- sample_meta[match(cols, sample_meta$sample_id), , drop = FALSE]
  
  if (anyNA(md$sample_id)) {
    missing_cols <- cols[is.na(md$sample_id)]
    stop("Some pseudobulk columns not found in sample_meta. Example: ", missing_cols[1])
  }
  
  # Define design (minimum needed for filterByExpr)
  # Here: keep genes expressed enough for comparing Condition (HD vs Control),
  # while allowing Region as a covariate.
  md$Condition <- factor(md$Condition)
  md$Region    <- factor(md$Region)
  
  design <- model.matrix(~ 0 + Condition + Region, data = md)
  
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)  # logical vector per gene
  
  counts_filt <- counts[keep, , drop = FALSE]
  pb_filt <- Matrix(counts_filt, sparse = TRUE)  # back to sparse
  
  out_file <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_filterByExpr.rds")
  )
  saveRDS(pb_filt, out_file)
  
  genes_before <- nrow(counts)
  genes_after  <- nrow(counts_filt)
  
  cat("Kept", genes_after, "genes out of", genes_before, "\n")
  cat("Saved:", out_file, "\n")
  
  summary <- rbind(summary, data.frame(
    file = basename(f),
    genes_before = genes_before,
    genes_after = genes_after,
    genes_removed = genes_before - genes_after,
    stringsAsFactors = FALSE
  ))
  
  rm(pb, counts, counts_filt, pb_filt, y, keep); gc()
}

write.csv(summary, file.path(out_dir, "filterByExpr_summary.csv"), row.names = FALSE)
cat("\nDone. Output folder:", out_dir, "\n")

obj <- readRDS("C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\split_celltypes_pseudobulk_rds_filterByExpr\\counts_Astrocyte_filtered_pseudobulk_filterByExpr.rds")

View(as.data.frame(as.matrix(obj)))