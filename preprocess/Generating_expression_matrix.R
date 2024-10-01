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
RNA_seurat <- readRDS("input_dir/RNA_seurat.rds")

#RNA_seurat@assays[["RNA"]]@data[1:5,1:5]
RNA_seurat <- NormalizeData(RNA_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
RNA_seurat <- FindVariableFeatures(RNA_seurat, selection.method = "vst", nfeatures = 2500)

##Save highly variable genes
HVG = VariableFeatures(RNA_seurat)[order(VariableFeatures(RNA_seurat))]

##Data standardization (centralization)
scale.genes <-  rownames(RNA_seurat)
RNA_seurat <- ScaleData(RNA_seurat, features = scale.genes)

#RNA_seurat@assays[["RNA"]]@scale.data[1:5,1:5]
cell_types = unique(as.character(RNA_seurat@meta.data$cell_type))
colnames(RNA_seurat@assays[["RNA"]]@scale.data) = RNA_seurat@meta.data$cell_type
#RNA_seurat@assays[["RNA"]]@scale.data[1:5,1:5]

##Create a list to store matrices for each cell type
cell_type_matrices <- list()
for (cell_type in cell_types) {
  cell_type_matrix <- RNA_seurat@assays[["RNA"]]@scale.data[, colnames(RNA_seurat@assays[["RNA"]]@scale.data) == cell_type]
  cell_type_matrices[[cell_type]] <- cell_type_matrix
}
rm(cell_type,cell_type_matrix)

##Create a data box to store the average vector
average_data <- data.frame(Column1 = 1:nrow(RNA_seurat@assays[["RNA"]]@scale.data))
for (cell_type in names(cell_type_matrices)) {
  cell_type_matrix <- cell_type_matrices[[cell_type]]
  average_data <- cbind(average_data, rowMeans(cell_type_matrix))
}
average_data[,1] <- NULL
colnames(average_data) = cell_types

##Obtain the final expression matrix of highly variable genes and sort them by gene name
final_exp_matrix = average_data[HVG,][order(rownames(average_data[HVG,])),]
colnames(final_exp_matrix) = 1:ncol(final_exp_matrix)
final_exp_matrix = t(as.matrix(final_exp_matrix))

write.table(final_exp_matrix, file = "output_dir/final_exp_matrix.txt",sep = "\t" ,col.names = NA,row.names = T,quote = F)

