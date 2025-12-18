
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
})

# ---------------- PATHS ----------------
HD1_RDS  <- "C:/Users/Avital/Desktop/AvitalLab/HD/HD2/pbmc_sctransform.RData"
REF_RDS  <- "C:/Users/Avital/Desktop/AvitalLab/HD/HD2/Matsushima_ref_SCT.rds"
OUT_RDS  <- "C:/Users/Avital/Desktop/AvitalLab/HD/HD2/HD2_mapped_to_Matsushima.rds"

message("Loading HD1 (query)…")
load("pbmc_sctransform.RData")   # loads object named pbmc
query <- pbmc                    # assign it to a new variable

message("Loading Matsushima reference…")
ref <- readRDS(REF_RDS)

DefaultAssay(query) <- "SCT"
DefaultAssay(ref)   <- "SCT"

# ---------------- PCA (if not done) ----------------
# Check whether PCA is present
if (!"pca" %in% names(query@reductions)) {
  message("Running PCA on HD1…")
  query <- RunPCA(query, npcs = 30, verbose = FALSE)
}

if (!"pca" %in% names(ref@reductions)) {
  message("Running PCA on reference…")
  ref <- RunPCA(ref, npcs = 30, verbose = FALSE)
}

# ---------------- FIND TRANSFER ANCHORS ----------------
message("Finding anchors...")
anchors <- FindTransferAnchors(
  reference = ref,
  query     = query,
  normalization.method = "SCT",
  dims      = 1:30,
  reference.reduction = "pca",
  verbose = TRUE
)

# ---------------- TRANSFER CELL TYPE LABELS ----------------
# Matsushima cell-type metadata column is typically "CellType" or "Subclass"
# Tell me if the metadata column name is different.
meta_col <- "CellType"


if (!(meta_col %in% colnames(ref@meta.data))) {
  stop(paste("Column", meta_col, "not found in reference metadata — check ref@meta.data."))
}

message("Transferring labels...")
pred <- TransferData(
  anchorset = anchors,
  ref@meta.data[, meta_col],
  dims      = 1:30
)

# Store predictions
query$predicted_celltype <- pred$predicted.id
query$prediction_score   <- pred$prediction.score.max

# ---------------- SAVE OUTPUT ----------------
saveRDS(query, OUT_RDS)

message("Finished mapping. Saved:")
message(OUT_RDS)

