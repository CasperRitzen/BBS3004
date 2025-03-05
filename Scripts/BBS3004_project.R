1+1
3+7
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(deseq2)
library(DESeq2)
library(dplyr)

library(GEOquery)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)
library(GEOquery)

#script
libary(tidyverse)
install.packages("tidyverse")
library(tidyverse)
library(DEseq2)
libary(DESeq2)

setwd("downloads/bbs3004/data")
library(DESeq2)
gse <- getGEO(GEO = 'FPKM, GSEMatrix = true')

# Load data
Data <- read.delim("FPKM_cufflinks.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
Data <- read.delim("FPKM_cufflinks.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)

# get Metadata, description of data
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata) # have a glimpse of how metadata looks like


metadata <- metadata %>%
  rename(
    Sample = 'tumor (t) or normal (n):ch1',
    Source = 'source_name_ch1',
    Tumor_stage = 'stage tnm:ch1',
    Age = 'age:ch1',
    Sex = 'gender:ch1',
    Life_Status = 'dead:ch1',
    Smoking_Status = 'smoking:ch1'
  )

# Check updated column names
print(colnames(metadata))

# Select only the renamed columns and overwrite metadata
metadata <- metadata %>%
  select(Sample, Source, Tumor_stage, Age, Sex, Life_Status, Smoking_Status)

# Check if the changes were applied
head(metadata)

# Remove the last row from Data
head(Data)
dim(Data)
Data <- Data[-nrow(Data), ]

# Check if the last row is removed
dim(Data)  # Check new dimensions

# Define the genes of interest
genes_of_interest <- c("ENSG00000146648", "ENSG00000118046", "ENSG00000157764")

# Convert Data to a long format (genes in rows, samples in columns)
expression_long <- Data %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% genes_of_interest) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")

# Merge expression data with metadata to include sample information
expression_long <- expression_long %>%
  left_join(metadata, by = "Sample")

# View the structure of the transformed data
print(head(expression_long))

# Plot expression levels of selected genes
ggplot(expression_long, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels

# Boxplot of gene expression by sex
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
ggplot(expression_long, aes(x = Tumor_stage, y = Expression, fill = Tumor_stage)) +
  geom_violin(alpha = 0.6) +  # Shows density of expression levels
  geom_jitter(width = 0.2, alpha = 0.7) +  # Adds individual sample points
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Tumor stage",
       x = "Tumor_stage",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(dplyr)

# Differential gene expression analysis

# Making sure row names in colData match to column names in counts_data
all(colnames(Data) %in% rownames(metadata))

# Find Columns in Data That Are Not in metadata
setdiff(colnames(Data), rownames(metadata))

# Remove everything after the underscore in Data
colnames(Data) <- sub("_.*", "", colnames(Data))

# Check again if row names in metadata matches to column names in Data
all(colnames(Data) %in% rownames(metadata))

# Check if they are in the same order
all(colnames(Data) == rownames(metadata))

# Reorder metadata rows to match the column order in Data
metadata <- metadata[match(colnames(Data), rownames(metadata)), , drop = FALSE]

# Check if they now match
all(colnames(Data) == rownames(metadata))

# Check the values in the data
summary(Data)

# Convert all data values to Absolute values. (Non-negative)
Data <- abs(Data)

# Round values to integers
Data <- round(Data)

# step 2: construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = Data,
                       coldata = metadata,
                       design = ~ Source)

# keeping rows that have at least 10 reads
keep <- rowSums(counts(dds)) >= 10

# set factor level
dds$Source <- relevel(dds$Source, ref = "Human non-malignant tissue")

# Step 3: run DESeq
dds <- DESeq(dds)
res <- results(dds)

# Explore results
summary (res)
res0.01 <- results(dds, alpha = 0.01)

# Contrasts
resultsNames(dds)
