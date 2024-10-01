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

# Get command-line parameters
args <- commandArgs(trailingOnly = TRUE)

# Check the number of parameters
if (length(args) < 2) {
  stop("Usage: Rscript Generation_expression_matrix.R <input_data_directory> <output_data_directory>")
}

input_dir <- args[1]
output_dir <- args[2]

bed_files <- list.files(path = "./input_dir",pattern = "\\.bed$")
bed_files <- mixedsort(bed_files)

whole_genome_200bp <- fread("input_dir/whole_genome_200bp.bed", col.names = c("Chromosome", "Start", "End"),header = FALSE,sep = "\t",data.table = F)
row_names <- paste(whole_genome_200bp$Chromosome, ":",whole_genome_200bp$Start, "-", whole_genome_200bp$End,sep = "")
whole_genome_200bp$bins = row_names
rownames(whole_genome_200bp) <- row_names

for (i in 1:length(bed_files)) {
  file <- bed_files[i]
  file1_name = paste("./input_dir/regression",file,sep = "")
  file2_name = paste("./input_dir/regress/",file,sep = "")
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

write.table(whole_genome_200bp, file = "./output_dir/union.peaks.labels.regress.txt", sep = "\t",col.names = NA,row.names = TRUE,quote = FALSE)

