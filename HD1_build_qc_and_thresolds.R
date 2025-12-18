#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

## ========= PATHS (edit if your paths differ) =========
BASE      <- "/ems/elsc-labs/meshorer-e/rivka.masar/HD/HD1"
IN_DIR    <- file.path(BASE, "HD1_samples_counts")
OUT_QC    <- file.path(BASE, "analysis_out", "QC_prethreshold")
OUT_THR   <- file.path(BASE, "analysis_out", "QC_thresholds")
QC_PERCELL_CSVGZ <- file.path(OUT_QC, "HD1_QC_percell.csv.gz")

dir.create(OUT_QC,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_THR, recursive = TRUE, showWarnings = FALSE)

## ========= HELPERS =========
read_count_csv_gz <- function(fp) {
  # Expect: first col = gene; other cols = cell barcodes
  dt <- fread(fp)
  gene_col <- names(dt)[1]
  genes <- dt[[gene_col]]
  dt[[gene_col]] <- NULL
  m <- as.matrix(dt)
  rownames(m) <- genes
  Matrix(m, sparse = TRUE)
}

## Compute per-cell QC without creating a Seurat object
qc_from_matrix <- function(mat, mt_regex="^MT-") {
  stopifnot(is(mat, "dgCMatrix"))
  nCount  <- Matrix::colSums(mat)
  nFeat   <- Matrix::colSums(mat > 0)
  mt_idx  <- grepl(mt_regex, rownames(mat))
  pctMT   <- if (any(mt_idx)) {
    mt_counts <- Matrix::colSums(mat[mt_idx, , drop=FALSE])
    100 * (mt_counts / pmax(1, nCount))
  } else {
    rep(0, length(nCount))
  }
  data.table(
    barcode        = colnames(mat),
    nCount_RNA     = as.numeric(nCount),
    nFeature_RNA   = as.numeric(nFeat),
    percent.mt     = as.numeric(pctMT)
  )
}

## ========= STEP A: Build per-cell QC file if missing =========
if (!file.exists(QC_PERCELL_CSVGZ)) {
  message("QC per-cell file not found. Building it now…")
  files <- list.files(IN_DIR, pattern="_mat\\.csv\\.gz$", full.names=TRUE)
  if (length(files) == 0) stop("No matrices found in: ", IN_DIR)

  # write header
  fwrite(data.table(sample=character(), barcode=character(),
                    nCount_RNA=integer(), nFeature_RNA=integer(), percent.mt=double()),
         QC_PERCELL_CSVGZ)

  for (i in seq_along(files)) {
    fp <- files[i]
    sample_id <- sub("_mat\\.csv\\.gz$", "", basename(fp))
    message(sprintf("[%d/%d] %s", i, length(files), sample_id))

    mat <- read_count_csv_gz(fp)
    qc  <- qc_from_matrix(mat, mt_regex="^MT-")
    qc[, sample := sample_id]
    setcolorder(qc, c("sample","barcode","nCount_RNA","nFeature_RNA","percent.mt"))

    fwrite(qc, QC_PERCELL_CSVGZ, append=TRUE)
    rm(mat, qc); gc()
  }
  message("Wrote per-cell QC to: ", QC_PERCELL_CSVGZ)
} else {
  message("Using existing per-cell QC: ", QC_PERCELL_CSVGZ)
}

## ========= STEP B: Load per-cell QC for threshold exploration =========
qc <- fread(QC_PERCELL_CSVGZ)
stopifnot(all(c("sample","barcode","nCount_RNA","nFeature_RNA","percent.mt") %in% names(qc)))
N <- nrow(qc)
message("Cells in QC table: ", format(N, big.mark=","))

## ========= STEP C: function that counts removals for given thresholds =========
count_removed <- function(min_feats = 0, max_feats = Inf,
                          min_counts = 0, max_counts = Inf,
                          max_pctMT = Inf) {
  keep <- (qc$nFeature_RNA >= min_feats) &
          (qc$nFeature_RNA <= max_feats) &
          (qc$nCount_RNA   >= min_counts) &
          (qc$nCount_RNA   <= max_counts) &
          (qc$percent.mt   <= max_pctMT)
  removed <- sum(!keep)
  list(
    removed = removed,
    kept    = sum(keep),
    frac_removed = removed / nrow(qc)
  )
}

## Save a tiny helper text with examples
example_txt <- file.path(OUT_THR, "how_to_use_count_removed.txt")
writeLines(c(
  "Examples in an interactive R session:",
  "  count_removed(min_feats=500, min_counts=1000, max_pctMT=5)",
  "  count_removed(min_feats=800, min_counts=1500, max_pctMT=10)",
  "",
  "Returned list: $removed, $kept, $frac_removed"
), con = example_txt)

## ========= STEP D: elbow-style scans & heatmaps =========

# 1) Percent.mt elbow
grid_mt <- seq(0, 20, by=1)           # you can extend to 30 if you want
mt_df <- data.table(
  max_pctMT = grid_mt,
  removed   = sapply(grid_mt, function(x) count_removed(max_pctMT = x)$removed)
)
p_mt <- ggplot(mt_df, aes(max_pctMT, removed/N)) +
  geom_line() + geom_point(size=1.2) +
  scale_y_continuous(labels=percent_format(accuracy = 1)) +
  labs(x="Max percent.mt (%)", y="% cells removed",
       title="Elbow scan: effect of max percent.mt") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_THR, "scan_percentMT_elbow.png"), p_mt, width=6, height=4, dpi=200)

# 2) nFeature_RNA minimum elbow (fix mt at 10%)
grid_feat <- seq(100, 4000, by=100)
feat_df <- data.table(
  min_feats = grid_feat,
  removed   = sapply(grid_feat, function(x) count_removed(min_feats = x, max_pctMT = 10)$removed)
)
p_feat <- ggplot(feat_df, aes(min_feats, removed/N)) +
  geom_line() + geom_point(size=1.2) +
  scale_y_continuous(labels=percent_format(accuracy = 1)) +
  labs(x="Min nFeature_RNA", y="% cells removed",
       title="Elbow scan: effect of min nFeature_RNA (max percent.mt = 10%)") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_THR, "scan_nFeature_min_elbow.png"), p_feat, width=6, height=4, dpi=200)

# 3) nCount_RNA minimum elbow (fix mt at 10%)
grid_counts <- seq(500, 15000, by=500)
count_df <- data.table(
  min_counts = grid_counts,
  removed    = sapply(grid_counts, function(x) count_removed(min_counts = x, max_pctMT = 10)$removed)
)
p_cnt <- ggplot(count_df, aes(min_counts, removed/N)) +
  geom_line() + geom_point(size=1.2) +
  scale_y_continuous(labels=percent_format(accuracy = 1)) +
  labs(x="Min nCount_RNA", y="% cells removed",
       title="Elbow scan: effect of min nCount_RNA (max percent.mt = 10%)") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_THR, "scan_nCount_min_elbow.png"), p_cnt, width=6, height=4, dpi=200)

# 4) 2D heatmap over (min_feats x max_pctMT), fixing min_counts=1000
grid_feat2 <- seq(200, 3000, by=200)
grid_mt2   <- seq(2, 20, by=2)
hm <- CJ(min_feats = grid_feat2, max_pctMT = grid_mt2)
hm[, frac_removed := vapply(
  seq_len(nrow(hm)),
  function(i) count_removed(min_feats = hm$min_feats[i],
                            min_counts = 1000,
                            max_pctMT  = hm$max_pctMT[i])$frac_removed,
  numeric(1)
)]
p_hm <- ggplot(hm, aes(min_feats, max_pctMT, fill = frac_removed)) +
  geom_tile() +
  scale_fill_viridis_c(labels = percent_format(accuracy = 1)) +
  labs(x="Min nFeature_RNA", y="Max percent.mt",
       fill="% removed",
       title="Fraction of cells removed (min nCount_RNA = 1000)") +
  theme_bw(base_size=11)
ggsave(file.path(OUT_THR, "heatmap_minFeats_vs_maxMT.png"), p_hm, width=6.5, height=4.8, dpi=200)

# 5) Print a few candidate settings to log
cands <- list(
  A = list(min_feats=500,  min_counts=1000, max_pctMT=5),
  B = list(min_feats=800,  min_counts=1500, max_pctMT=7),
  C = list(min_feats=1000, min_counts=2000, max_pctMT=10)
)
cat("\n=== Candidate thresholds ===\n")
for (nm in names(cands)) {
  x <- cands[[nm]]
  r <- count_removed(min_feats=x$min_feats, min_counts=x$min_counts, max_pctMT=x$max_pctMT)
  cat(sprintf("[%s] min_feats=%d  min_counts=%d  max_pctMT=%d  -> removed=%s (%.1f%%), kept=%s\n",
              nm, x$min_feats, x$min_counts, x$max_pctMT,
              format(r$removed, big.mark=","), 100*r$frac_removed,
              format(r$kept, big.mark=",")))
}
cat(sprintf("\nOutputs written to:\n - %s\n - %s\n", OUT_QC, OUT_THR))
