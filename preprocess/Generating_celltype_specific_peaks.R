#!/usr/bin/env Rscript
# ATAC-seq Data Processing Pipeline
# Creates Seurat object and generates cell-type specific chromatin accessibility profiles

# Environment Configuration ----------------------------------------------------
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())

# Command Line Argument Parsing ------------------------------------------------
library(argparse)
parser <- ArgumentParser(description = "ATAC-seq Data Processing Pipeline")
parser$add_argument("-i", "--input", 
                    help = "Directory containing input files", 
                    required = TRUE)
parser$add_argument("-o", "--output", 
                    help = "Output directory for results", 
                    required = TRUE)
args <- parser$parse_args()

input_dir <- normalizePath(args$input)
output_dir <- normalizePath(args$output)

# Create output subdirectory for BED files
bed_dir <- file.path(output_dir, "bed_class")
if (!dir.exists(bed_dir)) dir.create(bed_dir, recursive = TRUE)

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
message(Sys.time(), " Reading ATAC matrix...")
atac_matrix <- fread(file.path(input_dir, "Pbmc10k_ATAC_matrix.csv"), 
                     header = TRUE, sep = ',', data.table = FALSE)
rownames(atac_matrix) <- atac_matrix[,1]
atac_matrix[,1] <- NULL
message("Matrix dimensions: ", paste(dim(atac_matrix), collapse = " x "))
message("Matrix preview:")
print(atac_matrix[1:5, 1:5])

message(Sys.time(), " Reading ATAC metadata...")
atac_metadata <- fread(file.path(input_dir, "Pbmc10k_ATAC_metadata.csv"), 
                       header = TRUE, sep = ',', data.table = FALSE)
rownames(atac_metadata) <- atac_metadata[,1]
colnames(atac_metadata)[1] <- "cells"
message("Metadata dimensions: ", paste(dim(atac_metadata), collapse = " x "))
message("Metadata preview:")
print(head(atac_metadata))

# Seurat Object Creation -------------------------------------------------------
message(Sys.time(), " Creating ATAC Seurat object...")
atac_seurat <- CreateSeuratObject(
  counts = t(atac_matrix),
  meta.data = atac_metadata,
  assay = "ATAC",
  project = "Pbmc10k_ATAC",
  min.cells = 3,
  min.features = 300
)

message(Sys.time(), " Saving ATAC Seurat object...")
saveRDS(atac_seurat, file = file.path(output_dir, "ATAC_seurat.rds"))

# Cell-type Specific Analysis --------------------------------------------------
message(Sys.time(), " Processing cell-type specific chromatin accessibility...")

# Extract unique cell types
cell_types <- unique(atac_seurat@meta.data$cell_type)
message("Identified ", length(cell_types), " cell types: ", 
        paste(cell_types, collapse = ", "))

# Create list of accessibility matrices per cell type
message(Sys.time(), " Creating cell-type matrices...")
cell_type_matrices <- lapply(cell_types, function(ct) {
  ct_cells <- which(atac_seurat@meta.data$cell_type == ct)
  atac_seurat@assays$ATAC@data[, ct_cells, drop = FALSE]
})
names(cell_type_matrices) <- cell_types

# Filter peaks with sufficient accessibility per cell type
message(Sys.time(), " Filtering peaks by accessibility threshold...")
peak_celltypes_filtered <- lapply(cell_type_matrices, function(mat) {
  # Keep peaks accessible in at least 20% of cells
  keep_peaks <- rowSums(mat > 0) >= (ncol(mat) * 0.2)
  mat[keep_peaks, , drop = FALSE]
})

# Process each cell type's peaks into BED format
message(Sys.time(), " Converting peaks to BED format...")
for (i in seq_along(peak_celltypes_filtered)) {
  cell_type <- names(peak_celltypes_filtered)[i]
  message("Processing ", cell_type, " (", i, "/", length(peak_celltypes_filtered), ")")
  
  # Extract peak names and split genomic coordinates
  peak_names <- rownames(peak_celltypes_filtered[[cell_type]])
  split_names <- strsplit(peak_names, "[:-]")
  
  # Create BED data frame
  bed_df <- data.frame(
    chromosome = sapply(split_names, `[`, 1),
    start = sapply(split_names, `[`, 2),
    end = sapply(split_names, `[`, 3),
    stringsAsFactors = FALSE
  )
  
  # Write BED file
  output_file <- file.path(bed_dir, paste0(i, ".bed"))
  write.table(bed_df, file = output_file, 
              sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  message("Wrote ", nrow(bed_df), " peaks to ", output_file)
}

# Completion Message -----------------------------------------------------------
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs") %>% round(2)
message("\nPipeline completed successfully in ", duration, " seconds")
message("Output files saved to: ", output_dir)
message("BED files saved to: ", bed_dir)