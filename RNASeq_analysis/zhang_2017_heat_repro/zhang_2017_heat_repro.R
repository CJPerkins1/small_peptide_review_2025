#rna-seq analysis of Zhang 2017 hs on floral buds
#DGE analysis of HS on Arabidopsis

# Comparing DGE in late floral buds between control and HS conditions-----------------------------------------

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
library(knitr)

# Set working directory
setwd("~/Desktop/palanivelu_lab/small_peptide_review/R")

# Path to your gzipped file
expression_file <- "data/GSE94015_gene_readCount.txt.gz"

# Read the expression data as tab-delimited (tsv) file
count_data_full <- read_tsv(expression_file)

# View data dimensions and initial rows
dim(count_data_full)
head(count_data_full)

# Set gene names as rownames
count_data <- count_data_full[,-1]  # Exclude the first column with gene names

# Convert to matrix
count_data <- as.matrix(count_data)

# Set the row names from the first column
rownames(count_data) <- count_data_full[[1]]  # Set the first column as row names

# Remove columns with names starting with "bzip"
columns_to_remove2 <- grep("^bzip", colnames(count_data), value = TRUE)
count_data <- count_data[, !colnames(count_data) %in% columns_to_remove2]

# Check if row names are intact
head(count_data)

# Convert column names of count_data to lowercase
colnames(count_data) <- tolower(colnames(count_data))

# Remove underscores and convert column names to lowercase again (if necessary)
colnames(count_data) <- gsub("_", "", colnames(count_data))

# Verify the modified column names
head(colnames(count_data))


#load in metadata-------

# Specify the GSE ID
gse_id <- "GSE94015"

# Fetch the GEO data (expression data and metadata)
gse_data <- getGEO(gse_id, GSEMatrix = TRUE)

# Extract metadata (e.g., sample annotations)
meta_data <- pData(gse_data[[1]])

dim(meta_data)  # Check the metadata dimensions
head(meta_data)  # View first few rows of metadata

# Remove spaces in 'meta_data$title' and replace them with underscores
meta_data$title <- gsub(" ", "_", meta_data$title)

# Convert both 'meta_data$title' and 'ids' to lowercase
meta_data$title <- tolower(meta_data$title)

# Remove underscores and convert both to lowercase
meta_data$title <- tolower(gsub("_", "", meta_data$title))

#-----------------

# Check for missing values
missing_counts <- sum(is.na(count_data))
cat("Number of missing values in the read count table:", missing_counts)

# Replace Missing Values with Zero
# Replace NA with zero
count_data[is.na(count_data)] <- 0

# Filter Genes with Excessive Missing Values
# Define a threshold for missing data (20%)
threshold <- 0.2

# Calculate the proportion of missing values per gene
missing_proportion <- rowMeans(is.na(count_data))

# Filter out genes with missing values above the threshold
filtered_counts <- count_data[missing_proportion <= threshold, ]

# Filter for late floral buds metadata
# Subset columns with "lf" in their names
lf_count_data <- filtered_counts[, grepl("lf", colnames(count_data))]

lf_ids <- colnames(lf_count_data)

lf_group <- meta_data[which(meta_data$title == lf_ids),][['treatment:ch1']]# Modify as per your setup

lf_design <- model.matrix(~lf_group)

# Define the sample information as a data frame
lf_col_data <- data.frame(
  row.names = colnames(lf_count_data),
  condition = lf_group # Modify as per your setup
)

# Create DESeq2 dataset object
lf_dds <- DESeqDataSetFromMatrix(countData = lf_count_data, colData = lf_col_data, design = ~condition)
head(lf_dds)

# Apply DESeq’s normalization
lf_dds <- estimateSizeFactors(lf_dds)

# Inspect normalized counts (optional)
lf_norm_counts <- counts(lf_dds, normalized = TRUE)
formattable(head(lf_norm_counts))

# Run the DESeq2 pipeline
lf_dds <- DESeq(lf_dds)

# Extract results for all genes with default settings
lf_res <- results(lf_dds)

# Order by adjusted p-value
lf_res <- lf_res[order(lf_res$padj),]

#get meaning of the columns
mcols(lf_res, use.names = TRUE)

#other summary
summary(lf_res)

# Filter for significant results (adjusted p-value < 0.05)
lf_deseq2_genes <- subset(lf_res, padj < 0.05)

# Save to file
write.csv(as.data.frame(lf_deseq2_genes), file="output/HS_DESeq2_late_floral_results.csv")

## 1. Venn Diagram

###################################---------------------------cle_gene_list <- 

significant_genes <- list(
  DESeq2 = rownames(lf_deseq2_genes))

# # Create the Venn diagram
# grid.newpage()
# 
# venn.plot <- venn.diagram(
#   x = significant_genes,
#   category.names = c("DESeq2"),
#   filename = NULL,  # Set to NULL to plot directly
#   output = TRUE,
#   fill = c("cornflowerblue"),
#   alpha = 0.5,
#   cex = 1.5,
#   fontface = "bold",
#   fontfamily = "sans"
# )
# 
# # Display the Venn plot
# grid.draw(venn.plot)

# Example Data: Assume 'results' is a data frame containing columns logFC, P.Value (or PValue in edgeR), and adj.P.Val (or FDR)
# Add a column to categorize genes as upregulated, downregulated, or non-significant
results <- as.data.frame(lf_res)
results <- results %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Non-significant"
  ))

# Volcano Plot
ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-log10(P-Value)") +
  theme_minimal() +
  theme(legend.position = "top")

# Prepare Data: Select top 50 differentially expressed genes
results <- results %>%
  mutate(gene = rownames(results)) # Add gene names if rownames exist

top_genes <- results %>%
  filter(padj < 0.05) %>%                     # Filter significant genes
  arrange(desc(abs(log2FoldChange))) %>%      # Sort by absolute log2 fold change
  slice_head(n = 40) %>%                      # Select top 50 genes
  pull(gene)                                  # Extract gene names

#CLE gene family list
cle_family <- c("AT2G27250", "AT1G73165", "AT4G18510", "AT1G06225", "AT2G31081", "AT2G31083", "AT2G31085", "AT2G31082", 
                "AT1G67775", "AT1G26600", "AT1G69320", "AT1G49005", "AT1G68795", "AT1G73965", "AT1G63245", "AT2G01505", 
                "AT1G70895", "AT1G66145", "AT3G24225", "AT1G05065", "AT5G64800", "AT5G12235", "AT3G28455", "AT1G69970", 
                "AT3G25905", "AT5G12990", "AT3G24770", "AT2G34925", "AT1G25425", "AT4G13195", "AT1G69588")

ralf_family <- c(
  "AT1G02900", "AT1G23145", "AT1G23147", "AT1G28270", "AT1G35467", 
  "AT1G60625", "AT1G60815", "AT1G61563", "AT1G61566", "AT2G19020", 
  "AT2G19045", "AT2G20660", "AT2G22055", "AT2G32835", "AT2G32885", 
  "AT2G33130", "AT2G33775", "AT2G34825", "AT3G04735", "AT3G05490", 
  "AT3G16570", "AT3G23805", "AT3G25165", "AT3G25170", "AT3G29780", 
  "AT4G11510", "AT4G11653", "AT4G13075", "AT4G13950", "AT4G14010", 
  "AT4G15800", "AT5G67070", "AT1G60913", "AT2G32785"
)

# Filter DESeq2 results
cle_filtered_genes <- results %>%
  filter(gene %in% cle_family, padj < 0.05) %>%             # Filter for genes in the CLE family with padj < 0.05
  arrange(desc(abs(log2FoldChange)))                         # Arrange by absolute log2FoldChange

ralf_filtered_genes <- results %>%
  filter(gene %in% ralf_family, padj < 0.05) %>%             # Filter for genes in the CLE family with padj < 0.05
  arrange(desc(abs(log2FoldChange)))                         # Arrange by absolute log2FoldChange

# Table 1: Filtered DESeq2 results for CLE family
kable(head(cle_filtered_genes), caption = "Filtered DESeq2 Results for CLE Family")

# Table 2: Filtered DESeq2 results for RALF family
kable(head(ralf_filtered_genes), caption = "Filtered DESeq2 Results for RALF Family")

# Filter DESeq2 results for CLE and RALF genes without p-value filtering
cle_all_genes <- results %>%
  filter(gene %in% cle_family) %>%         # Filter for genes in the CLE family
  arrange(desc(abs(log2FoldChange)))      # Arrange by absolute log2FoldChange

ralf_all_genes <- results %>%
  filter(gene %in% ralf_family) %>%       # Filter for genes in the RALF family
  arrange(desc(abs(log2FoldChange)))      # Arrange by absolute log2FoldChange

# Expression Matrix Preparation
# Assuming 'dds' is the DESeq2 object and you have a normalized counts matrix
norm_counts <- counts(lf_dds, normalized = TRUE) # Get normalized counts
heatmap_data <- norm_counts[top_genes, ]      # Subset normalized counts by top genes

# Log Transform for Better Visualization
heatmap_data <- log2(heatmap_data + 1)        # Add pseudocount and log-transform

# Group Metadata
sample_groups <- as.data.frame(colData(lf_dds)[, "condition", drop = FALSE])
colnames(sample_groups) <- "Group" # Rename for clarity
rownames(sample_groups) <- colnames(heatmap_data)

# Custom Colors for Groups
group_colors <- list(Group = c("control" = "cyan", "treated at 37C for 3h" = "magenta"))

# Heatmap Visualization
pheatmap(heatmap_data,
         color = viridis(40),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",                      # Normalize rows (genes)
         show_rownames = TRUE,               # Show gene names
         show_colnames = TRUE,               # Show sample names
         annotation_col = sample_groups,     # Add group annotations
         annotation_colors = group_colors,   # Define group colors
         main = "Heatmap of Top Differentially Expressed Genes",
         fontsize_row = 7,                   # Adjust row font size
         fontsize_col = 7,                   # Adjust column font size
         angle_col = 45,
         fontsize = 8,
         #cellwidth = 20,                     # Width of each cell
         cellheight = 6) 

#---------------------------

# Normalized counts
norm_counts <- counts(lf_dds, normalized = TRUE)

# Custom Colors for Groups
group_colors <- list(Group = c("control" = "cyan", "treated at 37C for 3h" = "magenta"))

# RALF Genes ------------------------------------------------------------
# Define the mapping of gene IDs to common names for RALF
# Define the mapping of gene IDs to common names for RALF
ralf_gene_name_mapping <- c(
  "AT1G02900" = "RALF1", "AT1G23145" = "RALF2", "AT1G23147" = "RALF3", "AT1G28270" = "RALF4",
  "AT1G35467" = "RALF5", "AT1G60625" = "RALF6", "AT1G60815" = "RALF7", "AT1G61563" = "RALF8",
  "AT1G61566" = "RALF9", "AT2G19020" = "RALF10", "AT2G19045" = "RALF13", "AT2G20660" = "RALF14",
  "AT2G22055" = "RALF15", "AT2G32835" = "RALF16", "AT2G32885" = "RALF17", "AT2G33130" = "RALF18",
  "AT2G33775" = "RALF19", "AT2G34825" = "RALF20", "AT3G04735" = "RALF21", "AT3G05490" = "RALF22",
  "AT3G16570" = "RALF23", "AT3G23805" = "RALF24", "AT3G25165" = "RALF25", "AT3G25170" = "RALF26",
  "AT3G29780" = "RALF27", "AT4G11510" = "RALF28", "AT4G11653" = "RALF29", "AT4G13075" = "RALF30",
  "AT4G13950" = "RALF31", "AT4G14010" = "RALF32", "AT4G15800" = "RALF33", "AT5G67070" = "RALF34",
  "AT1G60913" = "RALF35", "AT2G32785" = "RALF36"
)


# RALF Genes -------------------------------------------------------------
# Filter out NAs and extract gene IDs for RALF
ralf_all_genes <- ralf_all_genes[!is.na(ralf_all_genes$gene), ]
ralf_gene_ids <- ralf_all_genes$gene

# Subset normalized counts for RALF genes
ralf_heatmap_data <- norm_counts[rownames(norm_counts) %in% ralf_gene_ids, ]

# Log-transform the data
ralf_heatmap_data <- log2(ralf_heatmap_data + 1)

# Remove rows with all zeros or constant values
ralf_heatmap_data <- ralf_heatmap_data[rowSums(ralf_heatmap_data > 0) > 0, ]

# Map gene IDs to gene_id(common_name) format
row_names <- rownames(ralf_heatmap_data)
new_row_names <- paste(row_names, "(", ralf_gene_name_mapping[row_names], ")", sep = "")

# Update the row names in the heatmap data
rownames(ralf_heatmap_data) <- new_row_names

# Check the number of valid genes for the heatmap
cat("Number of valid RALF genes for the heatmap:", nrow(ralf_heatmap_data), "\n")

# Generate Heatmap -----------------------------------------------------
# RALF Heatmap


# Capture the heatmap as a ggplot-like object
ralf_heatmap <- pheatmap(
  ralf_heatmap_data,
  color = viridis(40),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_groups,
  annotation_colors = group_colors,
  main = "Heatmap of RALF Genes",
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  fontsize = 8,
  filename = NA # Do not save immediately
)

# Save the heatmap directly as a PDF
ggsave("plots/zhang_2017_ralf_heatmap.pdf", plot = ralf_heatmap[[4]], width = 8, height = 6)

# CLE Genes -------------------------------------------------------------

# Define the mapping of gene IDs to common names
gene_name_mapping <- c(
  "AT2G27250" = "CLV3", "AT1G73165" = "CLE1", "AT4G18510" = "CLE2",
  "AT1G06225" = "CLE3", "AT2G31081" = "CLE4", "AT2G31083" = "CLE5",
  "AT2G31085" = "CLE6", "AT2G31082" = "CLE7", "AT1G67775" = "CLE8",
  "AT1G26600" = "CLE9", "AT1G69320" = "CLE10", "AT1G49005" = "CLE11",
  "AT1G68795" = "CLE12", "AT1G73965" = "CLE13", "AT1G63245" = "CLE14",
  "AT2G01505" = "CLE16", "AT1G70895" = "CLE17", "AT1G66145" = "CLE18",
  "AT3G24225" = "CLE19", "AT1G05065" = "CLE20", "AT5G64800" = "CLE21",
  "AT5G12235" = "CLE22", "AT3G28455" = "CLE25", "AT1G69970" = "CLE26",
  "AT3G25905" = "CLE27", "AT5G12990" = "CLE40", "AT3G24770" = "CLE41",
  "AT2G34925" = "CLE42", "AT1G25425" = "CLE43", "AT4G13195" = "CLE44",
  "AT1G69588" = "CLE45"
)

# CLE Genes -------------------------------------------------------------
# Filter out NAs and extract gene IDs for CLE
cle_all_genes <- cle_all_genes[!is.na(cle_all_genes$gene), ]
cle_gene_ids <- cle_all_genes$gene

# Subset normalized counts for CLE genes
cle_heatmap_data <- norm_counts[rownames(norm_counts) %in% cle_gene_ids, ]

# Log-transform the data
cle_heatmap_data <- log2(cle_heatmap_data + 1)

# Remove rows with all zeros or constant values
cle_heatmap_data <- cle_heatmap_data[rowSums(cle_heatmap_data > 0) > 0, ]

# Map gene IDs to gene_id(common_name) format
row_names <- rownames(cle_heatmap_data)
new_row_names <- paste(row_names, "(", gene_name_mapping[row_names], ")", sep = "")

# Update the row names in the heatmap data
rownames(cle_heatmap_data) <- new_row_names

# Check the number of valid genes for the heatmap
cat("Number of valid CLE genes for the heatmap:", nrow(cle_heatmap_data), "\n")

# Generate Heatmap -----------------------------------------------------

# Capture the heatmap as a ggplot-like object
cle_heatmap <- pheatmap(
  cle_heatmap_data,
  color = viridis(40),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_groups,
  annotation_colors = group_colors,
  main = "Heatmap of CLE Genes",
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  fontsize = 8,
  filename = NA # Do not save immediately
)

# Save the heatmap directly as a PDF
ggsave("plots/zhang_2017_cle_heatmap.pdf", plot = cle_heatmap[[4]], width = 8, height = 6)



#---------------------------

# Add common gene names to CLE genes
cle_all_genes <- cle_all_genes %>%
  mutate(common_name = gene_name_mapping[gene])

# Add common gene names to RALF genes
ralf_all_genes <- ralf_all_genes %>%
  mutate(common_name = ralf_gene_name_mapping[gene])

# Sort CLE genes by descending padj values
cle_all_genes <- cle_all_genes %>%
  arrange(padj)

# Sort RALF genes by descending padj values
ralf_all_genes <- ralf_all_genes %>%
  arrange(padj)

# Save CLE genes with common names to CSV
write.csv(cle_all_genes, "RNASeq_analysis/zhang_2017_heat_repro/cle_all_genes_zhang_2017.csv", row.names = FALSE)

# Save RALF genes with common names to CSV
write.csv(ralf_all_genes, "RNASeq_analysis/zhang_2017_heat_repro/ralf_all_genes_zhang_2017.csv", row.names = FALSE)

# Analyzing the characterized lists of reproduction and abiotic stress peptides -----------------------------

repro_peptides <- c("AT5G61605", "AT1G69588", "AT5G43285", "AT5G50423", 
                    "AT5G18403", "AT5G18407", "AT5G48605", "AT1G28270", 
                    "AT2G33775", "AT5G67070", "AT1G76750"
)

abiotic_peptides <- c("AT3G28455", "AT1G26600", "AT1G68765", "AT1G13590", 
                      "AT3G44735", "AT3G49780", "AT5G65870", "AT5G66815", 
                      "AT4G33720", "AT3G05490", "AT3G16570", "AT5G64905", 
                      "AT3G10930", "AT1G69588", "AT1G47485", "AT2G2344", 
                      "AT1G73165", "AT1G06225", "AT2G31081", "AT2G31082", 
                      "AT1G6324", "AT1G13620", "AT5G60810"
)

repro_filtered_genes <- results %>%
  filter(gene %in% repro_peptides, padj < 0.05) %>%             # Filter for genes in the CLE family with padj < 0.05
  arrange(desc(abs(log2FoldChange)))

abiotic_filtered_genes <- results %>%
  filter(gene %in% abiotic_peptides, padj < 0.05) %>%             # Filter for genes in the CLE family with padj < 0.05
  arrange(desc(abs(log2FoldChange)))

# Table 2: Filtered DESeq2 results for RALF family
kable(head(repro_filtered_genes), caption = "Filtered DESeq2 Results for RALF Family")

# Table 2: Filtered DESeq2 results for RALF family
kable(head(abiotic_filtered_genes), caption = "Filtered DESeq2 Results for RALF Family")

# Filter DESeq2 results for reproduction and abiotic genes without p-value filtering
repro_all_genes <- results %>%
  filter(gene %in% repro_peptides) %>%         # Filter for genes
  arrange(desc(abs(log2FoldChange)))      # Arrange by absolute log2FoldChange

abiotic_all_genes <- results %>%
  filter(gene %in% abiotic_peptides) %>%       # Filter for genes 
  arrange(desc(abs(log2FoldChange)))      # Arrange by absolute log2FoldChange

repro_name_mapping <- c(
  "AT5G61605" = "PCP-B⍺",
  "AT1G69588" = "CLE45",
  "AT5G43285" = "LURE1.1",
  "AT5G50423" = "XIUQIU1",
  "AT5G18403" = "XIUQIU2",
  "AT5G18407" = "XIUQIU3",
  "AT5G48605" = "XIUQIU4",
  "AT1G76750" = "EC1",
  "AT1G28270" = "RALF4",
  "AT2G33775" = "RALF19",
  "AT5G67070" = "RALF34"
)

abiotic_name_mapping <- c(
  "AT3G28455" = "CLE25",
  "AT1G26600" = "CLE9",
  "AT1G68765" = "IDA",
  "AT1G13590" = "PSKa1",
  "AT3G44735" = "PSKa2",
  "AT3G49780" = "PSKa3",
  "AT5G65870" = "PSKa4",
  "AT5G66815" = "CEP5",
  "AT4G33720" = "CAPE1",
  "AT3G05490" = "RALF22",
  "AT3G16570" = "RALF23",
  "AT5G64905" = "PEP3",
  "AT3G10930" = "IDL7",
  "AT1G69588" = "CLE45",
  "AT1G47485" = "CEP1",
  "AT2G2344"  = "CEP3",
  "AT1G73165" = "CLE1",
  "AT1G06225" = "CLE3",
  "AT2G31081" = "CLE4",
  "AT2G31082" = "CLE7",
  "AT1G6324"  = "CLE14",
  "AT1G13620" = "RGF2",
  "AT5G60810" = "RGF1"
)

repro_all_genes <- repro_all_genes[!is.na(repro_all_genes$gene), ]
repro_gene_ids <- repro_all_genes$gene

abiotic_all_genes <-abiotic_all_genes[!is.na(abiotic_all_genes$gene), ]
abiotic_gene_ids <- abiotic_all_genes$gene


# Add common gene names to CLE genes
repro_all_genes <- repro_all_genes %>%
  mutate(common_name = repro_name_mapping[gene])

# Add common gene names to RALF genes
abiotic_all_genes <- abiotic_all_genes %>%
  mutate(common_name = abiotic_name_mapping[gene])

# Sort CLE genes by descending padj values
repro_all_genes <- repro_all_genes %>%
  arrange(padj)

# Sort RALF genes by descending padj values
abiotic_all_genes <- abiotic_all_genes %>%
  arrange(padj)

# # Load necessary libraries
# library(DESeq2)
# library(edgeR)
# library(readr)
# library(GEOquery)
# library(vsn)
# library(ggplot2)
# library(pheatmap)
# library(formattable)
# library(VennDiagram)
# library(viridis)
# library(dplyr)
# library(stringr)
# 
# # Set working directory
# setwd("~/Desktop/palanivelu_lab/small_peptide_review/R")
# 
# # Path to your gzipped file
# expression_file <- "data/GSE94015_gene_readCount.txt.gz"
# 
# # Read the expression data as tab-delimited (tsv) file
# count_data_full <- read_tsv(expression_file)
# 
# # View data dimensions and initial rows
# dim(count_data_full)
# head(count_data_full)
# 
# # Set gene names as rownames
# count_data <- count_data_full[,-1]  # Exclude the first column with gene names
# 
# # Convert to matrix
# count_data <- as.matrix(count_data)
# 
# # Set the row names from the first column
# rownames(count_data) <- count_data_full[[1]]  # Set the first column as row names
# 
# # Remove columns with names starting with "bzip"
# columns_to_remove2 <- grep("^bzip", colnames(count_data), value = TRUE)
# count_data <- count_data[, !colnames(count_data) %in% columns_to_remove2]
# 
# # Check if row names are intact
# head(count_data)
# 
# # Convert column names of count_data to lowercase
# colnames(count_data) <- tolower(colnames(count_data))
# 
# # Remove underscores and convert column names to lowercase again (if necessary)
# colnames(count_data) <- gsub("_", "", colnames(count_data))
# 
# # Verify the modified column names
# head(colnames(count_data))
# 
# 
# #load in metadata-------
# 
# # Specify the GSE ID
# gse_id <- "GSE94015"
# 
# # Fetch the GEO data (expression data and metadata)
# gse_data <- getGEO(gse_id, GSEMatrix = TRUE)
# 
# # Extract metadata (e.g., sample annotations)
# meta_data <- pData(gse_data[[1]])
# 
# dim(meta_data)  # Check the metadata dimensions
# head(meta_data)  # View first few rows of metadata
# 
# 
# 
# 
# 
# 
# # Define experimental groups-------------------------------------------------you can remove ids conversions below if above changes work
# ids <- colnames(count_data)
# 
# # Remove spaces in 'meta_data$title' and replace them with underscores
# meta_data$title <- gsub(" ", "_", meta_data$title)
# 
# # Convert both 'meta_data$title' and 'ids' to lowercase
# meta_data$title <- tolower(meta_data$title)
# 
# # Remove underscores and convert both to lowercase
# meta_data$title <- tolower(gsub("_", "", meta_data$title))
# 
# group <- meta_data[which(meta_data$title == ids),][['treatment:ch1']]# Modify as per your setup
# 
# design <- model.matrix(~group)
# 
# 
# # Define the sample information as a data frame
# col_data <- data.frame(
#   row.names = colnames(count_data),
#   condition = group # Modify as per your setup
# )
# 
# # Create DESeq2 dataset object
# dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~condition)
# head(dds)
# 
# # Apply DESeq’s normalization
# dds <- estimateSizeFactors(dds)
# 
# # Inspect normalized counts (optional)
# norm_counts <- counts(dds, normalized = TRUE)
# formattable(head(norm_counts))
# 
# # Run the DESeq2 pipeline
# dds <- DESeq(dds)
# 
# # Extract results for all genes with default settings
# res <- results(dds)
# 
# # Order by adjusted p-value
# res <- res[order(res$padj),]
# 
# #get meaning of the columns
# mcols(res, use.names = TRUE)
# 
# #other summary
# summary(res)
# 
# # The first column, baseMean, is a just the average of the normalized count values, divided by the size factors, taken over all samples in the DESeqDataSet. The remaining four columns refer to a specific contrast, namely the comparison of the treatment level over the control level for the factor variable dex.
# # 
# # The column log2FoldChange is the effect size estimate. It tells us how much the gene’s expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples. This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene’s expression is increased by a multiplicative factor of 21.5≈2.82.
# # 
# # This estimate has an uncertainty associated with it, which is available in the column lfcSE, the standard error estimate for the log2 fold change estimate. The purpose of a test for differential expression is to test whether the data provides sufficient evidence to conclude that this value is really different from zero. DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability (i.e., the type of variability that you can expect between different samples in the same treatment group). As usual in statistics, the result of this test is reported as a p value, and it is found in the column pvalue.
# 
# # Filter for significant results (adjusted p-value < 0.05)
# deseq2_genes <- subset(res, padj < 0.05)
# 
# # Save to file
# write.csv(as.data.frame(deseq2_genes), file="output/HS_DESeq2_WTresults.csv")
# 
# ## 1. Venn Diagram
# 
# ###################################---------------------------cle_gene_list <- 
# 
# significant_genes <- list(
#   DESeq2 = rownames(deseq2_genes))
# 
# # # Create the Venn diagram
# # grid.newpage()
# # 
# # venn.plot <- venn.diagram(
# #   x = significant_genes,
# #   category.names = c("DESeq2"),
# #   filename = NULL,  # Set to NULL to plot directly
# #   output = TRUE,
# #   fill = c("cornflowerblue"),
# #   alpha = 0.5,
# #   cex = 1.5,
# #   fontface = "bold",
# #   fontfamily = "sans"
# # )
# # 
# # # Display the Venn plot
# # grid.draw(venn.plot)
# 
# # Example Data: Assume 'results' is a data frame containing columns logFC, P.Value (or PValue in edgeR), and adj.P.Val (or FDR)
# # Add a column to categorize genes as upregulated, downregulated, or non-significant
# results <- as.data.frame(res)
# results <- results %>%
#   mutate(Significance = case_when(
#     padj < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
#     padj < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
#     TRUE ~ "Non-significant"
#   ))
# 
# # Volcano Plot
# ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   scale_color_manual(values = c("blue", "grey", "red")) +
#   labs(title = "Volcano Plot of Differentially Expressed Genes",
#        x = "Log2 Fold Change",
#        y = "-log10(P-Value)") +
#   theme_minimal() +
#   theme(legend.position = "top")
# 
# # Prepare Data: Select top 50 differentially expressed genes
# results <- results %>%
#   mutate(gene = rownames(results)) # Add gene names if rownames exist
# 
# top_genes <- results %>%
#   filter(padj < 0.05) %>%                     # Filter significant genes
#   arrange(desc(abs(log2FoldChange))) %>%      # Sort by absolute log2 fold change
#   slice_head(n = 40) %>%                      # Select top 50 genes
#   pull(gene)                                  # Extract gene names
# 
# #CLE gene family list
# cle_family <- c("AT2G27250", "AT1G73165", "AT4G18510", "AT1G06225", "AT2G31081", "AT2G31083", "AT2G31085", "AT2G31082", 
#                 "AT1G67775", "AT1G26600", "AT1G69320", "AT1G49005", "AT1G68795", "AT1G73965", "AT1G63245", "AT2G01505", 
#                 "AT1G70895", "AT1G66145", "AT3G24225", "AT1G05065", "AT5G64800", "AT5G12235", "AT3G28455", "AT1G69970", 
#                 "AT3G25905", "AT5G12990", "AT3G24770", "AT2G34925", "AT1G25425", "AT4G13195", "AT1G69588")
# 
# ralf_family <- c(
#   "AT1G02900", "AT1G23145", "AT1G23147", "AT1G28270", "AT1G35467", 
#   "AT1G60625", "AT1G60815", "AT1G61563", "AT1G61566", "AT2G19020", 
#   "AT2G19045", "AT2G20660", "AT2G22055", "AT2G32835", "AT2G32885", 
#   "AT2G33130", "AT2G33775", "AT2G34825", "AT3G04735", "AT3G05490", 
#   "AT3G16570", "AT3G23805", "AT3G25165", "AT3G25170", "AT3G29780", 
#   "AT4G11510", "AT4G11653", "AT4G13075", "AT4G13950", "AT4G14010", 
#   "AT4G15800", "AT5G67070", "AT1G60913", "AT2G32785"
# )
# 
# repro_peptides <- c("AT5G61605", "AT1G69588", "AT5G43285", "AT5G50423", 
#                     "AT5G18403", "AT5G18407", "AT5G48605", "AT1G28270", 
#                     "AT2G33775", "AT5G67070", "AT1G76750"
# )
# 
# abiotic_peptides <- c("AT3G28455", "AT1G26600", "AT1G68765", "AT1G13590", 
#                       "AT3G44735", "AT3G49780", "AT5G65870", "AT5G66815", 
#                       "AT4G33720", "AT3G05490", "AT3G16570", "AT5G64905", 
#                       "AT3G10930", "AT1G69588", "AT1G47485", "AT2G2344", 
#                       "AT1G73165", "AT1G06225", "AT2G31081", "AT2G31082", 
#                       "AT1G6324", "AT1G13620", "AT5G60810"
# )
# 
# 
# # Filter DESeq2 results
# filtered_genes <- results %>%
#   filter(gene %in% cle_family, padj < 0.05) %>%             # Filter for genes in the CLE family with padj < 0.05
#   arrange(desc(abs(log2FoldChange)))                         # Arrange by absolute log2FoldChange
# 
# # Expression Matrix Preparation
# # Assuming 'dds' is the DESeq2 object and you have a normalized counts matrix
# norm_counts <- counts(dds, normalized = TRUE) # Get normalized counts
# heatmap_data <- norm_counts[top_genes, ]      # Subset normalized counts by top genes
# 
# # Log Transform for Better Visualization
# heatmap_data <- log2(heatmap_data + 1)        # Add pseudocount and log-transform
# 
# # Group Metadata
# sample_groups <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
# colnames(sample_groups) <- "Group" # Rename for clarity
# rownames(sample_groups) <- colnames(heatmap_data)
# 
# # Custom Colors for Groups
# group_colors <- list(Group = c("control" = "green", "treated at 37C for 3h" = "purple"))
# 
# # Heatmap Visualization
# pheatmap(heatmap_data,
#          color = viridis(40),
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          scale = "row",                      # Normalize rows (genes)
#          show_rownames = TRUE,               # Show gene names
#          show_colnames = TRUE,               # Show sample names
#          annotation_col = sample_groups,     # Add group annotations
#          annotation_colors = group_colors,   # Define group colors
#          main = "Heatmap of Top Differentially Expressed Genes",
#          fontsize_row = 7,                   # Adjust row font size
#          fontsize_col = 7,                   # Adjust column font size
#          angle_col = 45,
#          fontsize = 8,
#          #cellwidth = 20,                     # Width of each cell
#          cellheight = 6) 
# 
# 
# 
# # Now that I established that the pipeline works generally across all samples and conditions, let's specify into smaller comparisons




