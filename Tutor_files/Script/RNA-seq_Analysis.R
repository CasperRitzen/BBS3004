# This is script for RNA-seq Analysis

# Load Libraries for RNA-seq Data Analysis
library(DESeq2)
library(dplyr)
library(tidyverse)
library(GEOquery)
library(forcats)

# Set your working Directory
setwd("~/Research/PhD/Education/2025/BBS3004/Data/")

#=======================================#
# Step 1. Import Data & metadata into R #
#=======================================#

# This first part of the analysis uses FPKM normalised data

# Load your Data
FPKM_Data <- read.delim("GSE81089_FPKM_cufflinks.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)

# get metadata 
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)  # Have a glimpse of how metadata looks like

# Rename columns in metadata
colnames(metadata)

metadata <- metadata %>%
  rename(
    Sample = `tumor (t) or normal (n):ch1`,
    Source = source_name_ch1,
    Tumor_stage = `stage tnm:ch1`,
    Age = `age:ch1`,
    Sex = `gender:ch1`,
    Life_Status = `dead:ch1`,
    Smoking_Status = `smoking:ch1`
  )

# Check updated column names
colnames(metadata)

# Select only the renamed columns and overwrite metadata
metadata <- metadata %>%
  select(Sample, Source, Tumor_stage, Age, Sex, Life_Status, Smoking_Status) 

rownames(metadata) <- metadata$Sample

# Check if the changes were applied
head(metadata)

# Remove the last row from Data
head(FPKM_Data)
dim(FPKM_Data)
FPKM_Data <- FPKM_Data[-nrow(FPKM_Data), ]

# Check if the last row is removed
dim(FPKM_Data)  # Check new dimensions

# Explore Visualization of the data
# Perform quick visualization of known lung cancer genes:
# EGFR, KRAS, MET, LKB1, BRAF, PIK3CA, ALK, RET, and ROS1

# Making sure the row names in metadata matches to column names in FPKM_Data
all(colnames(FPKM_Data) %in% rownames(metadata))

# Find Columns in Data That Are Not in metadata
setdiff(colnames(FPKM_Data), rownames(metadata))

# Remove everything after the underscore in Data
colnames(FPKM_Data) <- sub("_.*", "", colnames(FPKM_Data))

# Check again if row names in metadata matches to column names in Data
all(colnames(FPKM_Data) %in% rownames(metadata))

# Replace NA values with "Control" in metadata
metadata <- metadata %>%
  mutate(across(everything(), ~replace_na(.x, "Control")))

# Define the genes of interest
genes_of_interest <- c("ENSG00000146648", "ENSG00000133703", "ENSG00000157764")

# Convert Data to a long format (genes in rows, samples in columns)
expression_long <- FPKM_Data %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% genes_of_interest) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")

# Merge expression data with metadata to include sample information
expression_long <- expression_long %>%
  left_join(metadata, by = "Sample")

# View the structure of the transformed data
head(expression_long)

# Plot expression levels of selected genes
ggplot(expression_long, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels

# Boxplot of gene expression grouped by sex
ggplot(expression_long, aes(x = Sex, y = Expression, fill = Sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Sex",
       x = "Sex",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability


# Alternatively, you can use a violin plot
ggplot(expression_long, aes(x = Sex, y = Expression, fill = Sex)) +
  geom_violin(alpha = 0.6) +  # Shows density of expression levels
  geom_jitter(width = 0.2, alpha = 0.7) +  # Adds individual sample points
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Sex",
       x = "Sex",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#==========================================================#
# Differential Gene Expression Analysis (A.K.A DEG Anaysis)#
#==========================================================#
 # Deseq2
# Load Raw counts file
Counts <- read.delim("Raw_Counts_GSE81089.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)

# making sure the row names in metadata matches to column names in Counts
all(colnames(Counts) %in% rownames(metadata))

# Check if they are in the same order
all(colnames(Counts) == rownames(metadata))

# Reorder metadata rows to match the column order in Counts
metadata <- metadata[match(colnames(Counts), rownames(metadata)), , drop = FALSE]

# Check if they now match
all(colnames(Counts) == rownames(metadata))

# Check the values in the Counts file
summary(Counts)

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = Counts,
  colData = metadata,
  design = ~ Source)

dds

# Quality control
# Remove genes with low counts
keep <- rowMeans(counts(dds)) >=10
dds <- dds[keep,]
dds
# Choose one either rowmeans or rowsum

# set the factor level
dds$Source <- relevel(dds$Source, ref = "Human non-malignant tissue")

# Run DESeq 
dds <- DESeq(dds)
res <- results(dds)

res

summary(res)
plotMA(res)
