#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================
# EDIT PATHS
# =========================
base_dir    <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2"
tmm_dir     <- file.path(base_dir, "split_celltypes_pseudobulk_edgeR_TMM")  # where *_TMM_DGEList.rds are
coldata_tsv <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")

out_dir <- file.path(base_dir, "split_celltypes_PCA_TMM_byRegion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load sample metadata (one row per sample = patient x region)
# =========================
colData <- read_tsv(coldata_tsv, show_col_types = FALSE)

req <- c("NBB_ID", "Region", "Condition", "Sex", "Age")
miss <- setdiff(req, colnames(colData))
if (length(miss) > 0) stop("coldata.tsv missing columns: ", paste(miss, collapse = ", "))

sample_meta <- colData %>%
  transmute(
    NBB_ID    = as.character(NBB_ID),
    Region    = as.character(Region),
    Condition = as.character(Condition),
    Sex       = as.character(Sex),
    Age       = suppressWarnings(as.numeric(Age))
  ) %>%
  distinct() %>%
  mutate(sample_id = paste(NBB_ID, Region, sep = "__"))

# =========================
# Infer cell type from filename
# e.g. counts_Mural_filtered_pseudobulk_filterByExpr_TMM_DGEList.rds -> Mural
# =========================
infer_celltype <- function(fname) {
  x <- tools::file_path_sans_ext(basename(fname))
  x <- sub("_TMM_DGEList$", "", x)
  x <- sub("^geneStats_", "", x)
  x <- sub("^counts_", "", x)
  x <- sub("_filtered_SCT$", "", x)
  x <- sub("_proteinCoding_SCT$", "", x)
  x <- sub("_filtered$", "", x)
  x <- sub("_proteinCoding$", "", x)
  x <- sub("_pseudobulk$", "", x)
  x <- sub("_filterByExpr$", "", x)
  x <- sub("_pseudobulk_filterByExpr$", "", x)
  x
}

# =========================
# PCA helper: choose genes (optional) + compute PCA
# =========================
run_pca <- function(logcpm_mat, top_var_genes = 3000) {
  # logcpm_mat: genes x samples
  
  # Optional: restrict to top variable genes for cleaner PCA
  if (!is.null(top_var_genes) && nrow(logcpm_mat) > top_var_genes) {
    v <- apply(logcpm_mat, 1, var, na.rm = TRUE)
    keep <- order(v, decreasing = TRUE)[seq_len(top_var_genes)]
    logcpm_mat <- logcpm_mat[keep, , drop = FALSE]
  }
  
  # prcomp expects samples x genes
  pca <- prcomp(t(logcpm_mat), center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  list(pca = pca, var_expl = var_expl)
}

# =========================
# Main loop over DGEList files
# =========================
files <- list.files(tmm_dir, pattern = "_TMM_DGEList\\.rds$", full.names = TRUE)

for (f in files) {
  cat("\n=== Processing:", basename(f), "===\n")
  
  y <- readRDS(f)
  
  # logCPM using TMM normalization factors inside y
  logcpm <- cpm(y, log = TRUE, prior.count = 1)  # genes x samples
  samples <- colnames(logcpm)
  
  # attach metadata
  md <- sample_meta[match(samples, sample_meta$sample_id), , drop = FALSE]
  if (anyNA(md$sample_id)) {
    missing_ids <- samples[is.na(md$sample_id)]
    stop("Some DGEList sample columns not found in sample_meta. Example: ", missing_ids[1])
  }
  
  ct <- infer_celltype(f)
  
  # make two plots: one per region
  for (reg in c("Caudate", "Putamen")) {
    
    keep_samp <- which(md$Region == reg)
    if (length(keep_samp) < 3) {
      cat("Skipping ", ct, " / ", reg, " (need >=3 samples for a meaningful PCA; found ",
          length(keep_samp), ")\n", sep = "")
      next
    }
    
    logcpm_reg <- logcpm[, keep_samp, drop = FALSE]
    md_reg <- md[keep_samp, , drop = FALSE]
    
    # PCA
    res <- run_pca(logcpm_reg, top_var_genes = 3000)
    pca <- res$pca
    var_expl <- res$var_expl
    
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    scores$sample_id <- rownames(scores)
    scores <- scores %>%
      left_join(md_reg, by = "sample_id")
    
    pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
    pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])
    
    # Plot: color by condition, shape by sex (easy to read)
    g <- ggplot(scores, aes(x = PC1, y = PC2, color = Condition, shape = Sex)) +
      geom_point(size = 3, alpha = 0.9) +
      labs(
        title = paste0("PCA (TMM logCPM) - ", ct, " - ", reg),
        x = pc1_lab,
        y = pc2_lab
      ) +
      theme_bw()
    
    out_png <- file.path(out_dir, paste0("PCA_TMM_", ct, "_", reg, ".png"))
    ggsave(out_png, g, width = 7, height = 5, dpi = 200)
    cat("Saved:", out_png, "\n")
  }
}

cat("\nDone. PCA plots in:\n", out_dir, "\n", sep = "")

