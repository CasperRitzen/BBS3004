# This is script for RNA-seq Analysis

# Load Libraries for RNA-seq Data Analysis
library(DESeq2)
library(dplyr)
library(tidyverse)
library(GEOquery)

# Set your working Directory
setwd("C:/Users/senag/OneDrive/Documenten/BMS year 3/BBS3004")

#=======================================#
# Step 1. Import Data & metadata into R #
#=======================================#

# Load your Data
Data <- read.delim("FPKM_cufflinks.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)

# get metadata 
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)  # Have a glimpse of how metadata looks like

# Increase buffer size by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

# Rename columns in metadata
colnames(metadata)
#write.table(metadata, file ="metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

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
print(colnames(metadata))

# Select only the renamed columns and overwrite metadata
metadata <- metadata %>%
  select(Sample, Source, Tumor_stage, Age, Sex, Life_Status, Smoking_Status) 

rownames(metadata) <- metadata$Sample

# Check if the changes were applied
head(metadata)
#write.table(metadata, file ="metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Remove the last row from Data
head(Data)
dim(Data)
Data <- Data[-nrow(Data), ]

# Check if the last row is removed
dim(Data)  # Check new dimensions

# Explore Visualization of the data
# Perform quick visualization of known lung cancer genes:
# EGFR, KRAS, MET, LKB1, BRAF, PIK3CA, ALK, RET, and ROS1

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
ggplot(expression_long, aes(x = Sex, y = Expression, fill = Sex)) +
  geom_violin(alpha = 0.6) +  # Shows density of expression levels
  geom_jitter(width = 0.2, alpha = 0.7) +  # Adds individual sample points
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Sex",
       x = "Sex",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#=======================================#
# Differential Gene Expression Analysis #
#=======================================#
# Deseq2
# making sure the row names in metadata matches to column names in Data
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

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = Data,
                              colData = metadata,
                              design = ~ Source)

dds

# Quality control
# Remove genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

keep2 <- rowMeans(counts(dds)) >=10
dds <- dds[keep2,]

# Choose one either rowmeans or rowsum

# set the factor level