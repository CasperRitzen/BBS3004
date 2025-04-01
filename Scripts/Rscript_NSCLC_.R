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

# Replace NA values with "Control" in metadata
metadata <- metadata %>%
  mutate(across(everything(), ~replace_na(.x, "Control")))

# Replace NA values with "Control" in all relevant columns
expression_long <- expression_long %>%
  mutate(across(everything(), ~replace_na(.x, "Control")))

# Ensure "Life_Status" specifically has no NA values
expression_long$Life_Status[is.na(expression_long$Life_Status)] <- "Control"

# Remove the "L790T" sample from expression_long
expression_long <- expression_long %>%
  filter(Sample != "L790T")

# Compute IQR for each gene and filter out outliers
expression_filtered <- expression_long %>%
  group_by(Gene) %>%  # Ensure outliers are removed per gene
  mutate(Q1 = quantile(Expression, 0.25),
         Q3 = quantile(Expression, 0.75),
         IQR = Q3 - Q1) %>%
  filter(Expression >= (Q1 - 1.5 * IQR) & Expression <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)  # Remove extra columns

# Plot expression levels of selected genes
ggplot(expression_filtered, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels

# Compute IQR for each gene and filter out outliers
expression_filtered <- expression_long %>%
  group_by(Gene) %>%  # Ensure outliers are removed per gene
  mutate(Q1 = quantile(Expression, 0.25),
         Q3 = quantile(Expression, 0.75),
         IQR = Q3 - Q1) %>%
  filter(Expression >= (Q1 - 1.5 * IQR) & Expression <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)  # Remove extra columns

# Plot expression levels of selected genes
ggplot(expression_filtered, aes(x = Life_Status, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples (Outliers Removed)",
       x = "Life_Status",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels


# Compute IQR for each gene and filter out outliers
expression_filtered <- expression_long %>%
  group_by(Gene) %>%  # Ensure outliers are removed per gene
  mutate(Q1 = quantile(Expression, 0.25),
         Q3 = quantile(Expression, 0.75),
         IQR = Q3 - Q1) %>%
  filter(Expression >= (Q1 - 1.5 * IQR) & Expression <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)  # Remove extra columns


# Boxplot of gene expression by Life_Status
ggplot(expression_filtered, aes(x = Life_Status, y = Expression, fill = Life_Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Life_Status",
       x = "Life_Status",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability

# Compute IQR for each gene and filter out outliers
expression_filtered <- expression_long %>%
  group_by(Gene) %>%  # Ensure outliers are removed per gene
  mutate(Q1 = quantile(Expression, 0.25),
         Q3 = quantile(Expression, 0.75),
         IQR = Q3 - Q1) %>%
  filter(Expression >= (Q1 - 1.5 * IQR) & Expression <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)  # Remove extra columns


# Alternatively, you can use a violin plot
ggplot(expression_filtered, aes(x = Life_Status, y = Expression, fill = Life_Status)) +
  geom_violin(alpha = 0.6) +  # Shows density of expression levels
  geom_jitter(width = 0.2, alpha = 0.7) +  # Adds individual sample points
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Life_Status",
       x = "Life_Status",
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
keep2 <- rowMeans(counts(dds)) >=10
dds <- dds[keep2,]

dds

# set the factor level
dds$Source <- relevel(dds$Source, ref = "Human non-malignant tissue")

# Run DESeq 
dds <- DESeq(dds) # Deseq runs its normalization, and compares the two specified condition/groups (i.e. Human NSCLC tissue vs Human non-malignant tissue... etc)
res <- results(dds) # This gives you the differentially Expressed Genes (DEGs)

res

summary(res)
plotMA(res)

#=======================================#
# Individual analysis Life Status #
#=======================================#

# Ensure column names in Data match row names in metadata
all(colnames(Data) %in% rownames(metadata))

# Find Columns in Data That Are Not in metadata
setdiff(colnames(Data), rownames(metadata))

# Remove everything after the underscore in Data
colnames(Data) <- sub("_.*", "", colnames(Data))

# Check again if row names in metadata match column names in Data
all(colnames(Data) %in% rownames(metadata))

# Ensure they are in the same order
all(colnames(Data) == rownames(metadata))

# Reorder metadata rows to match the column order in Data
metadata <- metadata[match(colnames(Data), rownames(metadata)), , drop = FALSE]

# Verify alignment
all(colnames(Data) == rownames(metadata))

# Summary of data
summary(Data)

# Convert all values in Data to absolute values (non-negative)
Data <- abs(Data)

# Round values to integers
Data <- round(Data)

# Construct a DESeqDataSet object with Life Status as the design factor
dds <- DESeqDataSetFromMatrix(countData = Data,
                              colData = metadata,
                              design = ~ Life_Status)

dds

# Quality control: Remove genes with low counts
keep2 <- rowMeans(counts(dds)) >= 10
dds <- dds[keep2,]

dds

# Set factor levels for Life_Status
dds$Life_Status <- relevel(dds$Life_Status, ref = "Control)  # Set reference level

# Run DESeq for differential expression analysis
dds <- DESeq(dds)

# Filter significant DEGs based on adjusted p-value < 0.01
res0.01 <- results(dds, alpha =0.01)

# Summary of differentially expressed genes (DEGs)
summary(res0.01)

# Plot MA plot for Life Status differential gene expression
plotMA(res0.01, main = "MA Plot - Differential Expression (Life Status)", ylim = c(-5, 5))

# Export DEG results to a file
res_df <- data.frame(Gene = rownames(res0.01), res0.01, row.names = NULL)  # Move row names to first column
write.table(res_df, file = "DEGs_LifeStatus.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
