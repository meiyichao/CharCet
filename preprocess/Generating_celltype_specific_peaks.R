setwd("/home/ycmei/model_demo/data")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE,scipen = 999)
rm(list = ls())

library(Seurat)
library(data.table)
library(tidyverse)
library(patchwork)
library(dplyr)

# Get command-line parameters
args <- commandArgs(trailingOnly = TRUE)

# Check the number of parameters
if (length(args) < 2) {
  stop("Usage: Rscript Generation_expression_matrix.R <input_data_directory> <output_data_directory>")
}

input_dir <- args[1]
output_dir <- args[2]

#################### Preprogressing the Seurat object ####################
ATAC_seurat <- readRDS("input_dir/ATAC_seurat.rds")

#ATAC_seurat@assays[["ATAC"]]@data[1:5,1:5]
cell_types = unique(as.character(ATAC_seurat@meta.data$cell_type))
colnames(ATAC_seurat@assays[["ATAC"]]@data) = ATAC_seurat@meta.data$cell_type
#ATAC_seurat@assays[["ATAC"]]@data[1:5,1:5]

##Create a list to store matrices for each cell type
cell_type_matrices <- list()
for (cell_type in cell_types) {
  cell_type_matrix <- ATAC_seurat@assays[["ATAC"]]@data[, colnames(ATAC_seurat@assays[["ATAC"]]@data) == cell_type]
  cell_type_matrices[[cell_type]] <- cell_type_matrix
}
rm(cell_type,cell_type_matrix)

##Create a list to store peaks reserved for different cell types
peak_celltypes_filtered = list()
for (cell_type in names(cell_type_matrices)) {
  cell_type_matrix <- cell_type_matrices[[cell_type]]
  row_counts= rowSums(cell_type_matrix > 0)
  keep_rows <- row_counts >= (ncol(cell_type_matrix))/5
  filtered_matrix <- cell_type_matrix[keep_rows,1]
  peak_celltypes_filtered[[cell_type]] <- as.data.frame(filtered_matrix)
}

i = 1
for (cell_type in names(peak_celltypes_filtered)) {
  split_names <- strsplit(as.character(row.names(peak_celltypes_filtered[[cell_type]])),"[:-]")
  chromsome = c()
  start = c()
  end = c()
  
  chromsome = sapply(split_names, function(x) chromsome=c(chromsome,x[1]))
  start = c(start,as.integer(sapply(split_names, function(x) start = c(start,as.numeric(x[2])))))
  end = c(end,as.integer(sapply(split_names, function(x) end = c(end,as.numeric(x[3])))))
  peak_celltypes_filtered[[cell_type]]$chromosome <- as.character(chromsome)
  peak_celltypes_filtered[[cell_type]]$start <- as.character(start)
  peak_celltypes_filtered[[cell_type]]$end <- as.character(end)
  peak_celltypes_filtered[[cell_type]][,1] <- NULL
  write.table(peak_celltypes_filtered[[cell_type]], file = paste("./output_dir/class",i,"_class.bed",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)  
  i =i +1
}


##Create a list to store peaks reserved for different cell types
peak_celltypes_filtered = list()
for (cell_type in names(cell_type_matrices)) {
  cell_type_matrix <- cell_type_matrices[[cell_type]]
  row_counts= rowSums(cell_type_matrix > 0)
  keep_rows <- row_counts >= (ncol(cell_type_matrix))/5
  filtered_matrix <- cell_type_matrix[keep_rows,1]
  peak_celltypes_filtered[[cell_type]] <- as.data.frame(filtered_matrix)
}

i = 1
for (cell_type in names(peak_celltypes_filtered)) {
  split_names <- strsplit(as.character(row.names(peak_celltypes_filtered[[cell_type]])),"[:-]")
  chromsome = c()
  start = c()
  end = c()
  
  chromsome = sapply(split_names, function(x) chromsome=c(chromsome,x[1]))
  start = c(start,as.integer(sapply(split_names, function(x) start = c(start,as.numeric(x[2])))))
  end = c(end,as.integer(sapply(split_names, function(x) end = c(end,as.numeric(x[3])))))
  
  peak_celltypes_filtered[[cell_type]]$chromosome <- as.character(chromsome)
  peak_celltypes_filtered[[cell_type]]$start <- as.character(start)
  peak_celltypes_filtered[[cell_type]]$end <- as.character(end)
  peak_celltypes_filtered[[cell_type]]$means <- peak_celltypes_filtered[[cell_type]]$row_means
  peak_celltypes_filtered[[cell_type]][,1] <- NULL
  write.table(peak_celltypes_filtered[[cell_type]], file = paste("./output_dir/regress",i,"_regress.bed",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)  
  i =i +1
}
