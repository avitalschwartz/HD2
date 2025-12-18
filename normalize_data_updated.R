#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(readr)
  library(dplyr)
})

#--------------------------------------------------
# Paths
#--------------------------------------------------
base_dir      <- "/ems/elsc-labs/meshorer-e/avital.schwartz1/HD/HD2"
counts_file   <- file.path(base_dir, "GSE152058_human_snRNA_processed_counts.mtx")
rowdata_file  <- file.path(base_dir, "GSE152058_human_snRNA_processed_rowdata.tsv")
coldata_file  <- file.path(base_dir, "GSE152058_human_snRNA_processed_coldata.tsv")

out_rds       <- file.path(base_dir, "pbmc_SCT.rds")

#--------------------------------------------------
# Load matrix + annotations
#--------------------------------------------------
message("Reading count matrix from: ", counts_file)
counts <- readMM(counts_file)

message("Reading rowdata (genes) from: ", rowdata_file)
rowdata <- read_tsv(rowdata_file, show_col_types = FALSE)

message("Reading coldata (cells) from: ", coldata_file)
coldata <- read_tsv(coldata_file, show_col_types = FALSE)

## Basic checks
message("Counts dim: ", paste(dim(counts), collapse = " x "))
message("Rowdata rows: ", nrow(rowdata))
message("Coldata rows: ", nrow(coldata))

if (nrow(counts) != nrow(rowdata)) {
  stop("Number of rows in counts (genes) does not match rowdata.")
}
if (ncol(counts) != nrow(coldata)) {
  message("WARNING: ncol(counts) != nrow(coldata). ",
          "Will generate generic cell names.")
}

#--------------------------------------------------
# Set feature and cell names
#--------------------------------------------------
# rowdata has columns: ENSEMBL, Gene, Feature  (from your screenshot)
# we will use gene SYMBOLS for rownames
if (!all(c("ENSEMBL", "Gene") %in% colnames(rowdata))) {
  stop("Expected columns ENSEMBL and Gene in rowdata.tsv.")
}

gene_symbols <- rowdata$Gene
gene_symbols[gene_symbols == "" | is.na(gene_symbols)] <- rowdata$ENSEMBL[gene_symbols == "" | is.na(gene_symbols)]

rownames(counts) <- make.unique(gene_symbols)

# cell names: if coldata has an obvious ID column, use it; otherwise make generic
cell_names <- NULL
possible_id_cols <- c("cell", "Cell", "barcode", "Barcode", "CellID", "cell_id")
id_col <- possible_id_cols[possible_id_cols %in% colnames(coldata)][1]

if (!is.null(id_col)) {
  cell_names <- coldata[[id_col]]
  cell_names <- make.unique(as.character(cell_names))
} else {
  cell_names <- paste0("Cell", seq_len(ncol(counts)))
  message("No obvious cell ID column in coldata; using generic names.")
}

colnames(counts) <- cell_names

#--------------------------------------------------
# Create Seurat object
#--------------------------------------------------
message("Creating Seurat object...")
pbmc <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(coldata))
message(sprintf("Object dimensions: %d genes x %d cells",
                nrow(pbmc), ncol(pbmc)))

#--------------------------------------------------
# Run SCTransform
#--------------------------------------------------
DefaultAssay(pbmc) <- "RNA"

message("Running SCTransform…")
pbmc <- SCTransform(
  pbmc,
  assay                 = "RNA",
  new.assay.name        = "SCT",
  variable.features.n   = 3000,
  vst.flavor            = "v2",
  conserve.memory       = TRUE,
  return.only.var.genes = FALSE,
  verbose               = TRUE
)

DefaultAssay(pbmc) <- "SCT"

#--------------------------------------------------
# Save
#--------------------------------------------------
saveRDS(pbmc, file = out_rds)

message("Done.")
message("SCT-normalized object saved to: ", out_rds)
message(sprintf("Assays present: %s", paste(Assays(pbmc), collapse = ", ")))
message("Example gene names: ", paste(head(rownames(pbmc)), collapse = ", "))
