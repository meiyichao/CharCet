#!/usr/bin/env Rscript
# Single Cell RNA-seq Data Processing Pipeline
# Creates Seurat object and generates cell-type specific expression matrix

# Environment Configuration ----------------------------------------------------
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())

# Command Line Argument Parsing ------------------------------------------------
library(argparse)
parser <- ArgumentParser(description = "Single Cell RNA-seq Processing Pipeline")
parser$add_argument("-i", "--input", 
                    help = "Directory containing input files", 
                    required = TRUE)
parser$add_argument("-o", "--output", 
                    help = "Output directory for results", 
                    required = TRUE)
parser$add_argument("-a", "--anno",
                    help = "Cell and annotation information file", 
                    required = TRUE)                    
args <- parser$parse_args()

input_dir <- normalizePath(args$input)
output_dir <- normalizePath(args$output)
anno <- normalizePath(args$anno)

dir.create(output_dir, recursive = TRUE)
# Library Import ---------------------------------------------------------------
start_time <- Sys.time()
message(Sys.time(), " Loading required packages...")
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(dplyr)
})

# Data Loading and Quality Checks ----------------------------------------------
message(Sys.time(), " Reading expression matrix...")
metadata <- fread(anno, header = TRUE, sep = "\t")
colnames(metadata) <- c("barcode", "celltype")
rna_matrix <- Read10X(data.dir = input_dir ,gene.column = 1)

# Seurat Object Creation -------------------------------------------------------
message(Sys.time(), " Creating Seurat object...")
rna_seurat <- CreateSeuratObject(
  counts = rna_matrix,
  assay = "RNA",
  project = "scRNA",
  min.cells = 3,
  min.features = 300
)

indice = match(rownames(rna_seurat@meta.data),metadata$barcode)
rna_seurat@meta.data$celltype = metadata$celltype[indice]

# Data Preprocessing -----------------------------------------------------------
message(Sys.time(), " Normalizing data...")
rna_seurat <- NormalizeData(
  rna_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

message(Sys.time(), " Identifying highly variable genes...")
rna_seurat <- FindVariableFeatures(
  rna_seurat,
  selection.method = "vst",
  nfeatures = 2500
)

# Extract and sort highly variable genes
hvg_genes <- VariableFeatures(rna_seurat) %>% sort()

message(Sys.time(), " Scaling data...")
rna_seurat <- ScaleData(rna_seurat, features = rownames(rna_seurat))

# Cell-type Specific Analysis --------------------------------------------------
message(Sys.time(), " Processing cell-type specific expression...")
cell_types <- unique(rna_seurat@meta.data$celltype)

# Create list of expression matrices per cell type
cell_type_matrices <- lapply(cell_types, function(ct) {
  ct_cells <- which(rna_seurat@meta.data$celltype == ct)
  rna_seurat@assays[["RNA"]]@layers[["data"]][, ct_cells, drop = FALSE]
})
names(cell_type_matrices) <- cell_types

# Calculate average expression per cell type
avg_expr <- sapply(cell_type_matrices, function(mat) rowMeans(mat))
rownames(avg_expr) <- rownames(rna_seurat)

# Filter and sort by highly variable genes
final_expr_matrix <- avg_expr[hvg_genes, ] %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  arrange(Gene) %>%
  column_to_rownames("Gene") %>%
  t()

# Output Results ---------------------------------------------------------------
message(Sys.time(), " Writing final expression matrix...")

celltype_id <- as.data.frame(row.names(final_expr_matrix))
write.table(celltype_id, 
            file = file.path(output_dir, "celltype_id.txt"),
            sep = "\t",
            col.names = F,
            quote = FALSE)
rownames(final_expr_matrix) <- 1:nrow(final_expr_matrix)
write.table(final_expr_matrix, 
            file = file.path(output_dir, "final_exp_matrix.txt"),
            sep = "\t",
            col.names = NA,
            quote = FALSE)

# Completion Message -----------------------------------------------------------
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs") %>% round(2)
message("Pipeline completed successfully in ", duration, " seconds")
