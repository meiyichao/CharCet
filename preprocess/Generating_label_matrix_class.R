#!/usr/bin/env Rscript
# ATAC-seq Open Chromatin Region Processing Pipeline
# Creates sparse matrix of open chromatin regions across single cells

# Environment Configuration ----------------------------------------------------
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())

# Command Line Argument Parsing ------------------------------------------------
library(argparse)
parser <- ArgumentParser(description = "ATAC-seq Open Chromatin Processing Pipeline")
parser$add_argument("-g", "--genome_bed", 
                    help = "Whole genome BED file", 
                    required = TRUE)
parser$add_argument("-i", "--input_dir", 
                    help = "Directory containing cell BED files", 
                    required = TRUE)
parser$add_argument("-o", "--output_dir", 
                    help = "Output directory for results", 
                    required = TRUE)
args <- parser$parse_args()

genome_bed <- normalizePath(args$genome_bed)
input_dir <- normalizePath(args$input_dir)
output_dir <- normalizePath(args$output_dir)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Library Import ---------------------------------------------------------------
start_time <- Sys.time()
message(Sys.time(), " Loading required packages...")
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(dplyr)
  library(gtools)
})

# Data Loading and Processing --------------------------------------------------
message(Sys.time(), " Reading genome BED file...")
whole_genome_200bp <- fread(
  genome_bed, 
  col.names = c("Chromosome", "Start", "End"),
  header = FALSE,
  sep = "\t",
  data.table = FALSE
)

# Expand peak regions to 1kb
message(Sys.time(), " Expanding peak regions...")
whole_genome_200bp <- whole_genome_200bp %>%
  mutate(Start = Start - 400,
         End = End + 400)

# Create row names for genomic regions
row_names <- paste(whole_genome_200bp$Chromosome, ":", 
                   whole_genome_200bp$Start, "-", 
                   whole_genome_200bp$End, sep = "")
rownames(whole_genome_200bp) <- row_names

# Process Cell BED Files -------------------------------------------------------
message(Sys.time(), " Collecting cell BED files...")
bed_files <- list.files(path = input_dir, pattern = "\\.bed$", full.names = TRUE)
bed_files <- mixedsort(bed_files)

# Initialize sparse matrix
message(Sys.time(), " Initializing sparse matrix...")
sparse_matrix <- Matrix(0, 
                        nrow = nrow(whole_genome_200bp), 
                        ncol = 0, 
                        sparse = TRUE)
rownames(sparse_matrix) <- row_names

message(Sys.time(), " Processing cell files...")
for (i in seq_along(bed_files)) {
  file <- bed_files[i]
  cell_id <- gsub("\\.bed$", "", basename(file))
  message(Sys.time(), " Processing cell: ", cell_id)
  
  # Read cell BED file
  cell_data <- fread(file, 
                     header = FALSE, 
                     col.names = c("chromosome", "start", "end", cell_id),
                     select = 1:4)
  
  # Add to sparse matrix
  sparse_matrix <- cbind(sparse_matrix, as(cell_data[[4]], "sparseMatrix"))
}

# Filter Matrix ----------------------------------------------------------------
message(Sys.time(), " Filtering matrix...")

# Keep only autosomes (chr1-chr22)
autosome_rows <- which(whole_genome_200bp$Chromosome %in% paste0("chr", 1:22))
sparse_matrix <- sparse_matrix[autosome_rows, ]

# Filter regions with at least one open cell
row_sum <- rowSums(sparse_matrix > 0)
selected_rows <- which(row_sum >= 1)
filtered_sparse <- sparse_matrix[selected_rows, , drop = FALSE]

# Convert to data frame for output
data_frame <- as.data.frame(as.matrix(filtered_sparse))
colnames(data_frame) <- seq_len(ncol(data_frame))

# Extract region information
region_info <- whole_genome_200bp[autosome_rows, ][selected_rows, ]

# Output Results ---------------------------------------------------------------
message(Sys.time(), " Writing results...")

# Write matrix file
write.table(data_frame, 
            file = file.path(output_dir, "union.peaks.labels.class.txt"),
            sep = "\t",
            col.names = NA,
            row.names = TRUE,
            quote = FALSE)

# Write BED file
write.table(region_info[, c("Chromosome", "Start", "End")],
            file = file.path(output_dir, "union.peaks.bed"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Completion Message -----------------------------------------------------------
end_time <- Sys.time()
duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
message(Sys.time(), " Pipeline completed successfully in ", duration, " seconds")