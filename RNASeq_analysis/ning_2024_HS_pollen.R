#analyzing pollen rnaseq data from SRA Ning et al. 2024 PRJNA1130680
# Install SRAdb if not already installed
if (!requireNamespace("SRAdb", quietly = TRUE)) {
  BiocManager::install("SRAdb")
}

# Load necessary libraries
library(DESeq2)
library(edgeR)
library(readr)
library(GEOquery)
library(vsn)
library(ggplot2)
library(pheatmap)
library(formattable)
library(VennDiagram)
library(viridis)
library(dplyr)
library(stringr)
library(SRAdb)

# Download the SRA metadata database
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(), sqlfile)

# Search for SRR or PRJNA
searchForTermInSRA("PRJNA1130680")  # Replace with your PRJNA ID

# Download specific SRR
getSRAfile("SRR12345678", sra_con, destDir = "data", fileType = "fastq")
