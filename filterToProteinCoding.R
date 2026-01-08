# --- Step 1 --- DO NOT RUN AGAIN

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Download protein-coding genes
protein_genes <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters    = "biotype",
  values     = "protein_coding",
  mart       = mart
)

protein_gene_list <- unique(protein_genes$hgnc_symbol)

length(protein_gene_list)   #  ~19–20 k


# --- Step 2 ---

library(Matrix)

# Input: matrices filtered by ≥10 cells
input_dir  <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\split_celltypes_counts_filtered"

# Output: matrices filtered to protein-coding genes only
output_dir <- "C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\split_celltypes_Counts_proteinCoding"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load all input files
files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Summary table
summary_table <- data.frame(
  file = character(),
  genes_before = integer(),
  genes_after = integer(),
  genes_removed = integer(),
  stringsAsFactors = FALSE
)

# Loop through each matrix
for (f in files) {
  
  mat <- readRDS(f)
  genes_before <- nrow(mat)
  
  # Filter to protein-coding only
  mat_pc <- mat[rownames(mat) %in% protein_gene_list , , drop = FALSE]
  genes_after <- nrow(mat_pc)
  
  # Save output
  out_file <- paste0(output_dir, "/", basename(f))
  saveRDS(mat_pc, out_file)
  
  # Add to summary table
  summary_table <- rbind(
    summary_table,
    data.frame(
      file = basename(f),
      genes_before = genes_before,
      genes_after = genes_after,
      genes_removed = genes_before - genes_after,
      stringsAsFactors = FALSE
    )
  )
  
  cat("Processed:", basename(f),
      "| Before:", genes_before,
      "| After:", genes_after,
      "| Removed:", genes_before - genes_after, "\n")
}

# Print summary
cat("\n==== Protein-Coding Gene Filtering Summary ====\n")
print(summary_table)

# Save summary as CSV
write.csv(summary_table,
          file = paste0(output_dir, "/protein_coding_filter_summary.csv"),
          row.names = FALSE)

cat("\nSummary saved to:", paste0(output_dir, "/protein_coding_filter_summary.csv"), "\n")
