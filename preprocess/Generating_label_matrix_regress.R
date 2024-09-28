setwd("/home/ycmei/model_demo/data")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE,scipen = 999)
rm(list = ls())

library(Seurat)
library(data.table)
library(tidyverse)
library(patchwork)
library(dplyr)
library(gtools)
library(Matrix)

ATAC_seurat <- readRDS("ATAC_seurat.rds")
ATAC_seurat@assays[["ATAC"]]@data[1:5,1:5]
cell_types = unique(as.character(ATAC_seurat@meta.data$cell_type))
colnames(ATAC_seurat@assays[["ATAC"]]@data) = ATAC_seurat@meta.data$cell_type
ATAC_seurat@assays[["ATAC"]]@data[1:5,1:5]


cell_type_matrices <- list()
for (cell_type in cell_types) {
  cell_type_matrix <- ATAC_seurat@assays[["ATAC"]]@data[, colnames(ATAC_seurat@assays[["ATAC"]]@data) == cell_type]
  cell_type_matrices[[cell_type]] <- cell_type_matrix
}
rm(cell_type,cell_type_matrix)

peak_celltypes_filtered = list()
for (cell_type in names(cell_type_matrices)) {
  cell_type_matrix <- cell_type_matrices[[cell_type]]
  row_counts= rowSums(cell_type_matrix > 0)
  keep_rows <- row_counts >= (ncol(cell_type_matrix))/5
  filtered_matrix <- cell_type_matrix[keep_rows,]
  row_means <- apply(filtered_matrix, 1, mean)
  peak_celltypes_filtered[[cell_type]] <- as.data.frame(row_means)
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
  write.table(peak_celltypes_filtered[[cell_type]], file = paste("./bed_regress/",i,".bed",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)  
  i =i +1
}

bed_files <- list.files(path = "./output_bed_regress",pattern = "\\.bed$")
bed_files <- mixedsort(bed_files)

whole_genome_200bp <- fread("/home/ycmei/model_demo/data/whole_genome_200bp.bed", col.names = c("Chromosome", "Start", "End"),header = FALSE,sep = "\t",data.table = F)
row_names <- paste(whole_genome_200bp$Chromosome, ":",whole_genome_200bp$Start, "-", whole_genome_200bp$End,sep = "")
whole_genome_200bp$bins = row_names
rownames(whole_genome_200bp) <- row_names


for (i in 1:length(bed_files)) {
  file <- bed_files[i]
  file1_name = paste("./output_bed_regress/",file,sep = "")
  file2_name = paste("./bed_regress/",file,sep = "")
  cell_id <- gsub("\\.bed$", "", file)
  print(cell_id)
  data1 <- fread(file1_name, header = FALSE,data.table = F)
  data2 <- fread(file2_name, header = FALSE,data.table = F)
  data1$means <- NA
  match_indices <- match(data1$V5, data2$V2)
  data1$means[!is.na(match_indices)] <- data2$V4[match_indices[!is.na(match_indices)]]
  data1 = data1[,-c(4,5,6)]
  data1_names <- paste(data1$V1, ":",data1$V2, "-", data1$V3,sep = "")
  data1$bins = data1_names
  data1 <- data1[, 4:5]
  data1 <- data1[, c(ncol(data1), 1:(ncol(data1)-1))]
  new_column_name <- paste0("V",i)
  whole_genome_200bp[[new_column_name]] <- NA
  match_indi <- match(whole_genome_200bp$bins, data1$bins)
  whole_genome_200bp[[new_column_name]][!is.na(match_indi)] <- data1$means[match_indi[!is.na(match_indi)]]
}

##Only obtain the rows containing 22 autosomes
whole_genome_200bp = whole_genome_200bp[1:(which(whole_genome_200bp$Chromosome == "chrX")-1),]
dim(whole_genome_200bp)
whole_genome_200bp = whole_genome_200bp[,-c(1,2,3,4)]
whole_genome_200bp[is.na(whole_genome_200bp)] <- 0
##Calculate the number of non-zero elements per row
row_sum <- rowSums(whole_genome_200bp > 0)
##Find the row index that meets the criteria (at least one is open, so we will keep this bin)
selected_rows <- which(row_sum >= 1)
whole_genome_200bp <- whole_genome_200bp[selected_rows, , drop = FALSE]
dim(whole_genome_200bp)


row_names <- rownames(whole_genome_200bp)
split_names <- strsplit(row_names, "[:-]")
info_df <- do.call(rbind, lapply(split_names, function(x) {
  data.frame(chromosome = x[1], start = as.integer(x[2]), end = as.integer(x[3]))
}))

info_df = info_df %>%
  mutate(start = start - 400,
         end = end + 400)
row_names <- paste(info_df$chromosome, ":",info_df$start, "-", info_df$end,sep = "")
rownames(whole_genome_200bp) <- row_names
colnames(whole_genome_200bp) <- 1:19

write.table(whole_genome_200bp, file = "./union.peaks.labels.regress.txt", sep = "\t",col.names = NA,row.names = TRUE,quote = FALSE)

