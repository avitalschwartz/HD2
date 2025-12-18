#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
})

# ---------------- PATHS ----------------
HD1_RDS  <- "/ems/elsc-labs/meshorer-e/avital.schwartz1/HD/HD2/pbmc_SCT.rds"
REF_RDS  <- "/ems/elsc-labs/meshorer-e/avital.schwartz1/HD/HD2/Matsushima_ref_SCT.rds"
OUT_RDS  <- "/ems/elsc-labs/meshorer-e/avital.schwartz1/HD/HD2/HD1_mapped_to_Matsushima.rds"
# ---------------- LOAD HD1 QUERY ----------------
message("Loading HD1 (query) from .rds…")
query <- readRDS(HD1_RDS)

if (!inherits(query, "Seurat")) {
  stop("ERROR: Query object is not a Seurat object.")
}
message("Loaded HD1 successfully.")

# ---------------- LOAD REFERENCE ----------------
message("Loading Matsushima reference…")
ref <- readRDS(REF_RDS)

head(rownames(query))
head(rownames(ref))

if (!inherits(ref, "Seurat")) {
  stop("ERROR: Reference object is not a Seurat object.")
}
message("Loaded reference successfully.")

# ---------------- SET ASSAYS ---------------
DefaultAssay(query) <- "SCT"
DefaultAssay(ref)   <- "SCT"


# ---------------- RUN PCA IF MISSING ----------------
if (!"pca" %in% names(query@reductions)) {
  message("Running PCA on HD1 query…")
  query <- RunPCA(query, npcs = 30, verbose = FALSE)
}

if (!"pca" %in% names(ref@reductions)) {
  message("Running PCA on reference…")
  ref <- RunPCA(ref, npcs = 30, verbose = FALSE)
}
Assays(ref)
Assays(query)
length(intersect(rownames(ref), rownames(query)))
length(VariableFeatures(ref))
length(VariableFeatures(query))
# ---------------- FIND TRANSFER ANCHORS ----------------
message("Finding transfer anchors (SCT)…")
anchors <- FindTransferAnchors(
  reference = ref,
  query     = query,
  normalization.method = "SCT",
  dims      = 1:30,
  reference.reduction = "pca",
  verbose   = TRUE
)

# ---------------- TRANSFER CELL TYPE LABELS ----------------
meta_col <- "CellType"   # Adjust if your reference uses a different column name!

if (!(meta_col %in% colnames(ref@meta.data))) {
  stop(paste("ERROR: Metadata column", meta_col, "not found in reference."))
}

message("Transferring cell-type annotations…")
pred <- TransferData(
  anchorset = anchors,
  refdata   = ref@meta.data[[meta_col]],
  dims      = 1:30
)

query$predicted_celltype <- pred$predicted.id
query$prediction_score   <- pred$prediction.score.max

# ---------------- SAVE OUTPUT ----------------
message("Saving mapped object to:")
message(OUT_RDS)
saveRDS(query, OUT_RDS)

message("DONE. Mapping completed successfully.")



