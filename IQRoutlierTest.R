#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
})

# =========================
# PATHS
# =========================
base_dir    <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2"
tmm_dir     <- file.path(base_dir, "split_celltypes_pseudobulk_edgeR_TMM")  # *_TMM_DGEList.rds
coldata_tsv <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")

out_dir <- file.path(base_dir, "split_celltypes_IQR_outliers_PCA_TMM_byRegion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TOP_VAR_GENES_FOR_PCA <- 3000  # use top variable genes for PCA (recommended)

# =========================
# SAMPLE METADATA (one row per patient x region)
# =========================
colData <- read_tsv(coldata_tsv, show_col_types = FALSE)

req <- c("NBB_ID", "Region", "Condition")
miss <- setdiff(req, colnames(colData))
if (length(miss) > 0) stop("coldata.tsv missing columns: ", paste(miss, collapse = ", "))

sample_meta <- colData %>%
  transmute(
    NBB_ID    = as.character(NBB_ID),
    Region    = as.character(Region),
    Condition = as.character(Condition),
    sample_id = paste(NBB_ID, Region, sep = "__")
  ) %>%
  distinct()

# =========================
# Helpers
# =========================
infer_celltype <- function(fname) {
  x <- tools::file_path_sans_ext(basename(fname))
  x <- sub("_TMM_DGEList$", "", x)
  x <- sub("^counts_", "", x)
  x <- sub("_pseudobulk_filterByExpr$", "", x)
  x <- sub("_filterByExpr$", "", x)
  x <- sub("_pseudobulk$", "", x)
  x <- sub("_filtered$", "", x)
  x <- sub("_proteinCoding$", "", x)
  x
}

pick_top_var_genes <- function(mat, top_n) {
  if (is.null(top_n) || nrow(mat) <= top_n) return(rownames(mat))
  v <- apply(mat, 1, var, na.rm = TRUE)
  rownames(mat)[order(v, decreasing = TRUE)[seq_len(top_n)]]
}

# Classic IQR outlier rule
iqr_test <- function(x, k = 1.5) {
  x <- x[is.finite(x)]
  if (length(x) < 4) {
    return(list(
      lower = NA_real_, upper = NA_real_,
      q1 = NA_real_, q3 = NA_real_, iqr = NA_real_
    ))
  }
  q1 <- as.numeric(quantile(x, 0.25, names = FALSE, type = 7))
  q3 <- as.numeric(quantile(x, 0.75, names = FALSE, type = 7))
  iqr <- q3 - q1
  lower <- q1 - k * iqr
  upper <- q3 + k * iqr
  list(lower = lower, upper = upper, q1 = q1, q3 = q3, iqr = iqr)
}

# =========================
# MAIN LOOP
# =========================
files <- list.files(tmm_dir, pattern = "_TMM_DGEList\\.rds$", full.names = TRUE)

all_outliers <- data.frame()

for (f in files) {
  cat("\n=== Processing:", basename(f), "===\n")
  
  y <- readRDS(f)
  ct <- infer_celltype(f)
  
  # TMM-normalized logCPM (genes x samples)
  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  samp_ids <- colnames(logcpm)
  
  md <- sample_meta[match(samp_ids, sample_meta$sample_id), , drop = FALSE]
  if (anyNA(md$sample_id)) {
    bad <- samp_ids[is.na(md$sample_id)]
    stop("Metadata mismatch for samples. Example: ", bad[1])
  }
  
  # region-wise PCA (you asked earlier for separate Caudate/Putamen)
  for (reg in c("Caudate", "Putamen")) {
    
    keep <- which(md$Region == reg)
    if (length(keep) < 4) {
      cat("Skipping ", ct, " / ", reg, " (need >=4 samples for stable IQR; found ",
          length(keep), ")\n", sep = "")
      next
    }
    
    X <- logcpm[, keep, drop = FALSE]
    md_reg <- md[keep, , drop = FALSE]
    
    genes_use <- pick_top_var_genes(X, TOP_VAR_GENES_FOR_PCA)
    X <- X[genes_use, , drop = FALSE]
    
    # PCA on samples (samples x genes)
    pca <- prcomp(t(X), center = TRUE, scale. = TRUE)
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    scores$sample_id <- rownames(scores)
    
    scores <- scores %>%
      left_join(md_reg, by = "sample_id") %>%
      mutate(
        PCdist = sqrt(PC1^2 + PC2^2)
      )
    
    # IQR fences on PC1, PC2, and distance
    f1 <- iqr_test(scores$PC1)
    f2 <- iqr_test(scores$PC2)
    fd <- iqr_test(scores$PCdist)
    
    scores <- scores %>%
      mutate(
        out_PC1   = is.finite(f1$lower) & (PC1 < f1$lower | PC1 > f1$upper),
        out_PC2   = is.finite(f2$lower) & (PC2 < f2$lower | PC2 > f2$upper),
        out_dist  = is.finite(fd$lower) & (PCdist < fd$lower | PCdist > fd$upper),
        cell_type = ct,
        region    = reg
      )
    
    out <- scores %>%
      filter(out_PC1 | out_PC2 | out_dist) %>%
      transmute(
        cell_type, region, sample_id, patient_id = NBB_ID, condition = Condition,
        PC1, PC2, PCdist,
        out_PC1, out_PC2, out_dist,
        PC1_lower = f1$lower, PC1_upper = f1$upper,
        PC2_lower = f2$lower, PC2_upper = f2$upper,
        dist_lower = fd$lower, dist_upper = fd$upper
      )
    
    if (nrow(out) > 0) {
      out_file <- file.path(out_dir, paste0("IQR_outliers_", ct, "_", reg, ".csv"))
      write.csv(out, out_file, row.names = FALSE)
      cat("Wrote outliers:", out_file, " (n=", nrow(out), ")\n", sep = "")
      all_outliers <- rbind(all_outliers, out)
    } else {
      cat("No IQR outliers for ", ct, " / ", reg, "\n", sep = "")
    }
  }
  
  rm(y, logcpm); gc()
}

# combined summary
sum_file <- file.path(out_dir, "IQR_outliers_ALL_celltypes.csv")
write.csv(all_outliers, sum_file, row.names = FALSE)
cat("\nDone. Outlier tables in:\n", out_dir, "\nSummary:\n", sum_file, "\n", sep = "")
