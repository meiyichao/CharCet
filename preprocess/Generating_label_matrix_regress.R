#!/usr/bin/env Rscript
# ATAC-seq Peak Processing Pipeline
# Processes single-cell ATAC-seq data to generate cell-type specific peak matrices

# Environment Configuration ----------------------------------------------------
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())

# Command Line Argument Parsing ------------------------------------------------
library(argparse)
parser <- ArgumentParser(description = "ATAC-seq Peak Processing Pipeline")
parser$add_argument("-i", "--input_rds", 
                    help = "Path to ATAC_seurat.rds file", 
                    required = TRUE)
parser$add_argument("-g", "--genome_bed", 
                    help = "Path to whole_genome_200bp.bed file", 
                    required = TRUE)
parser$add_argument("-b", "--bed_regress_dir", 
                    help = "Directory for intermediate BED files", 
                    required = TRUE)
parser$add_argument("-o", "--output_dir", 
                    help = "Output directory for results", 
                    required = TRUE)
args <- parser$parse_args()

input_rds <- normalizePath(args$input_rds)
genome_bed <- normalizePath(args$genome_bed)
bed_regress_dir <- normalizePath(args$bed_regress_dir)
output_dir <- normalizePath(args$output_dir)

dir.create(bed_regress_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Library Import ---------------------------------------------------------------
start_time <- Sys.time()
message(Sys.time(), " Loading required packages...")
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(gtools)
  library(Matrix)
})

# Data Loading and Preparation -------------------------------------------------
message(Sys.time(), " Reading ATAC Seurat object...")
ATAC_seurat <- readRDS(input_rds)

# Extract cell type information
message(Sys.time(), " Processing cell types...")
cell_types <- unique(as.character(ATAC_seurat@meta.data$cell_type))

# Create cell-type specific matrices
message(Sys.time(), " Creating cell-type specific matrices...")
cell_type_matrices <- lapply(cell_types, function(ct) {
  ct_cells <- which(ATAC_seurat@meta.data$cell_type == ct)
  ATAC_seurat@assays[["ATAC"]]@data[, ct_cells, drop = FALSE]
})
names(cell_type_matrices) <- cell_types

# Filter Peaks and Generate BED Files ------------------------------------------
message(Sys.time(), " Filtering peaks and generating BED files...")
peak_celltypes_filtered <- list()

for (cell_type in names(cell_type_matrices)) {
  message(Sys.time(), " Processing: ", cell_type)
  
  # Extract matrix for current cell type
  ct_matrix <- cell_type_matrices[[cell_type]]
  
  # Filter peaks present in at least 20% of cells
  row_counts <- rowSums(ct_matrix > 0)
  keep_rows <- row_counts >= (ncol(ct_matrix)) / 5
  filtered_matrix <- ct_matrix[keep_rows, ]
  
  # Calculate mean accessibility
  row_means <- Matrix::rowMeans(filtered_matrix)
  peak_celltypes_filtered[[cell_type]] <- data.frame(
    row_means = row_means,
    peak_id = names(row_means)
  )
}

# Generate BED files for each cell type
message(Sys.time(), " Writing intermediate BED files...")
for (i in seq_along(peak_celltypes_filtered)) {
  cell_type <- names(peak_celltypes_filtered)[i]
  df <- peak_celltypes_filtered[[cell_type]]
  
  # Extract genomic coordinates from peak IDs
  split_names <- strsplit(df$peak_id, "[:-]")
  df$chromosome <- sapply(split_names, `[`, 1)
  df$start <- as.integer(sapply(split_names, `[`, 2))
  df$end <- as.integer(sapply(split_names, `[`, 3))
  
  # Prepare BED format
  bed_df <- df[, c("chromosome", "start", "end", "row_means")]
  
  # Write to file
  write.table(bed_df, 
              file = file.path(bed_regress_dir, paste0(i, ".bed")),
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE)
}

# Create Genome-wide Peak Matrix -----------------------------------------------
message(Sys.time(), " Creating genome-wide peak matrix...")

# Read whole genome BED file
whole_genome_200bp <- fread(
  genome_bed, 
  col.names = c("Chromosome", "Start", "End"),
  header = FALSE,
  sep = "\t",
  data.table = FALSE
)

# Create bin identifiers
row_names <- paste(whole_genome_200bp$Chromosome, ":", 
                   whole_genome_200bp$Start, "-", 
                   whole_genome_200bp$End, sep = "")
whole_genome_200bp$bins <- row_names
rownames(whole_genome_200bp) <- row_names

# Process intermediate BED files
bed_files <- list.files(path = bed_regress_dir, pattern = "\\.bed$", full.names = TRUE)
bed_files <- mixedsort(bed_files)

for (i in seq_along(bed_files)) {
  file <- bed_files[i]
  cell_id <- gsub("\\.bed$", "", basename(file))
  message(Sys.time(), " Processing: ", cell_id)
  
  # Read intermediate BED file
  data2 <- fread(file, header = FALSE, data.table = FALSE)
  
  # Create bin identifiers for intermediate data
  data2_bins <- paste(data2$V1, ":", data2$V2, "-", data2$V3, sep = "")
  
  # Add column to genome dataframe
  col_name <- paste0("V", i)
  whole_genome_200bp[[col_name]] <- NA
  
  # Match and assign values
  match_idx <- match(whole_genome_200bp$bins, data2_bins)
  whole_genome_200bp[[col_name]][!is.na(match_idx)] <- data2$V4[match_idx[!is.na(match_idx)]]
}

# Filter and Finalize Matrix ---------------------------------------------------
message(Sys.time(), " Filtering and finalizing matrix...")

# Keep only autosomes (chr1-chr22)
autosome_rows <- grep("^chr[1-9]|^chr1[0-9]|^chr2[0-2]", whole_genome_200bp$Chromosome)
whole_genome_200bp <- whole_genome_200bp[autosome_rows, ]

# Remove unnecessary columns
whole_genome_200bp <- whole_genome_200bp[, !names(whole_genome_200bp) %in% 
                                           c("Chromosome", "Start", "End", "bins")]

# Replace NA with 0
whole_genome_200bp[is.na(whole_genome_200bp)] <- 0

# Filter rows with at least one non-zero value
row_sum <- rowSums(whole_genome_200bp > 0)
selected_rows <- which(row_sum >= 1)
filtered_matrix <- whole_genome_200bp[selected_rows, , drop = FALSE]

# Adjust genomic coordinates
message(Sys.time(), " Adjusting genomic coordinates...")
split_names <- strsplit(rownames(filtered_matrix), "[:-]")
info_df <- data.frame(
  chromosome = sapply(split_names, `[`, 1),
  start = as.integer(sapply(split_names, `[`, 2)) - 400,
  end = as.integer(sapply(split_names, `[`, 3)) + 400
)

# Create new row names
new_row_names <- paste(info_df$chromosome, ":", 
                       info_df$start, "-", 
                       info_df$end, sep = "")
rownames(filtered_matrix) <- new_row_names
colnames(filtered_matrix) <- seq_len(ncol(filtered_matrix))

# Output Results ---------------------------------------------------------------
message(Sys.time(), " Writing final output...")
write.table(filtered_matrix, 
            file = file.path(output_dir, "union.peaks.labels.regress.txt"),
            sep = "\t",
            col.names = NA,
            row.names = TRUE,
            quote = FALSE)

# Completion Message -----------------------------------------------------------
end_time <- Sys.time()
duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
message(Sys.time(), " Pipeline completed successfully in ", duration, " seconds")