setwd("/home/ycmei/model_demo/data")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE,scipen = 999)
rm(list = ls())

library(data.table)
library(gtools)
library(Matrix)
library(dplyr)

# Get command-line parameters
args <- commandArgs(trailingOnly = TRUE)

# Check the number of parameters
if (length(args) < 2) {
  stop("Usage: Rscript Generation_expression_matrix.R <input_data_directory> <output_data_directory>")
}

input_dir <- args[1]
output_dir <- args[2]

whole_genome_200bp <- fread("input_dir/whole_genome_200bp.bed", col.names = c("Chromosome", "Start", "End"),header = FALSE,sep = "\t",data.table = F)

##Expand the range of peaK regions by 400bp each, making the range of each peak region 1kb
whole_genome_200bp = whole_genome_200bp %>%
  mutate(Start = Start - 400,
         End = End + 400)

row_names <- paste(whole_genome_200bp$Chromosome, ":",whole_genome_200bp$Start, "-", whole_genome_200bp$End,sep = "")
rownames(whole_genome_200bp) <- row_names

##########  Store open genomic information of all cells in the form of sparse matrices  ##########
bed_files <- list.files(path = "./input_dir",pattern = "\\.bed$")
bed_files <- mixedsort(bed_files)
sparse_matrix <- Matrix(0, nrow = nrow(whole_genome_200bp), ncol =0, sparse = TRUE)
rownames(sparse_matrix) = row_names

for (i in 1:length(bed_files)) {
  file <- bed_files[i]
  file_name = paste("./input_dir/",file,sep = "")
  cell_id <- gsub("\\.bed$", "", file)
  print(cell_id)
  data <- fread(file_name, header = FALSE, col.names = c("chromosome", "start", "end", cell_id))
  new_column <- as(data[[4]], "sparseMatrix")
  sparse_matrix <- cbind(sparse_matrix,new_column)
}

##Only obtain the rows containing 22 autosomes
sparse_matrix = sparse_matrix[1:(which(whole_genome_200bp$Chromosome == "chrX")-1),]
dim(sparse_matrix)
##Calculate the number of non-zero elements per row
row_sum <- rowSums(sparse_matrix > 0)
##Find the row index that meets the criteria (at least one is open, so we will keep this bin)
selected_rows <- which(row_sum >= 1)
filtered_sparse <- sparse_matrix[selected_rows, , drop = FALSE]
data_frame = as.data.frame(as.matrix(filtered_sparse))

colnames(data_frame) = 1:ncol(data_frame)
write.table(data_frame, file = "output_dir/union.peaks.labels.class.txt", sep = "\t",col.names = NA,row.names = TRUE,quote = FALSE)


row_names <- rownames(data_frame)
split_names <- strsplit(row_names, "[:-]")
info_df <- do.call(rbind, lapply(split_names, function(x) {
  data.frame(chromosome = x[1], start = as.integer(x[2]), end = as.integer(x[3]))
}))
write.table(info_df, file = "output_dir/union.peaks.bed",sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)  
