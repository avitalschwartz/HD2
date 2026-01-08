#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

base_dir    <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2"
tmm_dir     <- file.path(base_dir, "split_celltypes_pseudobulk_edgeR_TMM")
coldata_tsv <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")
out_dir     <- file.path(base_dir, "split_celltypes_Boxplot_logCPM_byRegion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# plot density control (avoid millions of points)
TOP_VAR_GENES_TO_PLOT <- 3000   # set NULL to plot all genes (can be huge)
POINT_ALPHA <- 0.05
POINT_SIZE  <- 0.20

colData <- read_tsv(coldata_tsv, show_col_types = FALSE)

sample_meta <- colData %>%
  transmute(
    NBB_ID    = as.character(NBB_ID),
    Region    = as.character(Region),
    Condition = as.character(Condition),
    sample_id = paste(NBB_ID, Region, sep="__")
  ) %>%
  distinct()

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

files <- list.files(tmm_dir, pattern = "_TMM_DGEList\\.rds$", full.names = TRUE)

for (f in files) {
  cat("\n=== Processing:", basename(f), "===\n")
  y <- readRDS(f)
  
  # TMM-normalized logCPM
  logcpm <- cpm(y, log = TRUE, prior.count = 1)   # genes x samples
  samples <- colnames(logcpm)
  
  md <- sample_meta[match(samples, sample_meta$sample_id), , drop = FALSE]
  if (anyNA(md$sample_id)) stop("Metadata mismatch for ", basename(f))
  
  ct <- infer_celltype(f)
  
  for (reg in c("Caudate", "Putamen")) {
    
    keep <- which(md$Region == reg)
    if (length(keep) < 2) next
    
    logcpm_reg <- logcpm[, keep, drop = FALSE]
    md_reg <- md[keep, , drop = FALSE]
    
    genes_use <- pick_top_var_genes(logcpm_reg, TOP_VAR_GENES_TO_PLOT)
    logcpm_reg <- logcpm_reg[genes_use, , drop = FALSE]
    
    genes <- rownames(logcpm_reg)
    ng <- length(genes)
    ns <- ncol(logcpm_reg)
    
    df <- data.frame(
      patient_id = rep(md_reg$NBB_ID, each = ng),
      gene_id    = rep(genes, times = ns),
      logCPM     = as.vector(logcpm_reg),
      stringsAsFactors = FALSE
    )
    
    df$patient_id <- factor(df$patient_id, levels = unique(md_reg$NBB_ID))
    
    g <- ggplot(df, aes(x = patient_id, y = logCPM)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.25, height = 0),
                 alpha = POINT_ALPHA, size = POINT_SIZE) +
      labs(
        title = paste0("Per-patient logCPM distribution - ", ct, " - ", reg),
        x = "Patient ID",
        y = "TMM logCPM (per gene)"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    out_png <- file.path(out_dir, paste0("BOX_logCPM_", ct, "_", reg, ".png"))
    ggsave(out_png, g, width = 12, height = 5, dpi = 200)
    cat("Saved:", out_png, " | genes:", ng, " | patients:", ns, "\n")
  }
  
  rm(y, logcpm); gc()
}

cat("\nDone. Output:", out_dir, "\n")
