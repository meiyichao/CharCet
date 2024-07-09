setwd("/home/ycmei/model_demo/data")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE,scipen = 999)
rm(list = ls())

library(Seurat)
library(data.table)
library(tidyverse)
library(patchwork)
library(dplyr)

#################### Creating a Seurat object ####################

message("Reading counts...")
ATAC_matrix = fread("Pbmc10k_ATAC_matrix.csv",header = T,sep = ',',data.table = F)
rownames(ATAC_matrix) <- ATAC_matrix[,1]
ATAC_matrix[,1] <- NULL
print(dim(ATAC_matrix))
print(ATAC_matrix[1:5,1:5])

message("Reading metadata...")
ATAC_metadata <- fread("Pbmc10k_ATAC_metadata.csv",header=T,sep = ',',data.table = F)
rownames(ATAC_metadata) <- ATAC_metadata[,1]
colnames(ATAC_metadata)[1] <- "cells"
print(dim(ATAC_metadata))
print(head(ATAC_metadata))

message("Writing seurat object...")
ATAC_seurat = CreateSeuratObject(counts=t(ATAC_matrix),meta.data=ATAC_metadata,assay="ATAC",project="Pbmc10k_ATAC",min.cells=3,min.features=300)
saveRDS(ATAC_seurat, file = "ATAC_seurat.rds")

#################### Preprogressing the Seurat object ####################
ATAC_seurat <- readRDS("ATAC_seurat.rds")

ATAC_seurat@assays[["ATAC"]]@data[1:5,1:5]
cell_types = unique(as.character(ATAC_seurat@meta.data$cell_type))
colnames(ATAC_seurat@assays[["ATAC"]]@data) = ATAC_seurat@meta.data$cell_type
ATAC_seurat@assays[["ATAC"]]@data[1:5,1:5]

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
  write.table(peak_celltypes_filtered[[cell_type]], file = paste("./bed_class/",i,".bed",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)  
  i =i +1
}

