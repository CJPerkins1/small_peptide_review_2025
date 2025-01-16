#RNA-seq workshop 11/12/2024

# Load necessary libraries
library(DESeq2)
library(edgeR)
library(limma)
library(readr)
library(vsn)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(dplyr)
library(factoextra)
library(formattable)
library(apeglm)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("apegl")

# Import your read count table from Github

urlfile="https://raw.githubusercontent.com/merlinis12/RNA-Seq-Data-Analysis-in-R/refs/heads/main/data/airway_scaledcounts.csv"
count_data_full <- read_csv(url(urlfile))

# View data dimensions and initial rows
dim(count_data_full)
head(count_data_full)

#gene names as rownames
count_data <- count_data_full[,-1]
rownames(count_data) <- count_data_full[[1]]


#read metadata from github
urlfile="https://raw.githubusercontent.com/merlinis12/RNA-Seq-Data-Analysis-in-R/refs/heads/main/data/airway_metadata.csv"
meta_data <- read_csv(url(urlfile))
meta_data$dex <- as.factor(meta_data$dex)

# View data dimensions and initial rows
dim(meta_data)
head(meta_data)

# Define experimental groups
ids <- colnames(count_data)
group <- meta_data[which(meta_data$id == ids),][['dex']]# Modify as per your setup
design <- model.matrix(~group)
