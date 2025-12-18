library(Matrix)
mat <- readMM("C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\GSE152058_human_snRNA_processed_counts.mtx")

library(readr)
colData <- read_tsv("C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\GSE152058_human_snRNA_processed_coldata.tsv")  
rowData <- read_tsv("C:\\Users\\Avital\\Desktop\\AvitalLab\\HD\\HD2\\GSE152058_human_snRNA_processed_rowdata.tsv")  

library(dplyr)

unique(colData$NBB_ID)

samples_to_check <- c("A47L", "3345", "2030", "2665", "4254", "4308", "4621", "4294",
                      "3881", "A39R", "4494", "2952", "3730", "2903", "A58R", "4225")

counts_by_person_region <- colData %>% count(NBB_ID, Region, name = "n_cells")
counts_by_person_region %>% arrange(n_cells)

for (id in samples_to_check) {
  cat("\n---", id, "---\n")
  print(colData %>% filter(NBB_ID == id) %>% count(Region, name = "n_cells") %>% arrange(desc(n_cells)))
}

threshold <- 1000

group_sizes <- colData %>%
  count(NBB_ID, Region, name = "n_cells")

drop_groups <- group_sizes %>%
  filter(n_cells < threshold) %>%
  arrange(n_cells)

print(drop_groups)  # shows which person×region will be removed

# Row indices (cells) to DROP / KEEP (works even if mat has no colnames)
drop_idx <- colData %>%
  mutate(.idx = row_number()) %>%
  semi_join(drop_groups, by = c("NBB_ID", "Region")) %>%
  pull(.idx)

keep_idx <- setdiff(seq_len(nrow(colData)), drop_idx)

length(drop_idx)
length(keep_idx)

rownames(mat) <- rowData$Gene
colnames(mat) <- colData$Barcode