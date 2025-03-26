
# Load Libraries for RNA-seq Data Analysis
library(DESeq2)
library(dplyr)
library(tidyverse)
library(GEOquery)

# Set your working Directory
setwd("/Users/sofiabaars/Downloads/BBS3004/Data/")

#=======================================#
# Step 1. Import Data & metadata into R #
#=======================================#

# Load your Data
Data <- read.delim("FPKM_cufflinks.tsv", row.names=1, header=TRUE, sep="\t", check.names = FALSE)

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

# Compute IQR for each gene and filter out outliers
expression_filtered <- expression_long %>%
  group_by(Gene) %>%  # Ensure outliers are removed per gene
  mutate(Q1 = quantile(Expression, 0.25),
         Q3 = quantile(Expression, 0.75),
         IQR = Q3 - Q1) %>%
  filter(Expression >= (Q1 - 1.5 * IQR) & Expression <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)  # Remove extra columns

# Rename NA into control and convert it into a factor
metadata$Tumor_stage[is.na(metadata$Tumor_stage)] <- "Control"
metadata$Tumor_stage <- as.factor(metadata$Tumor_stage)

# Plot expression levels of selected genes
ggplot(expression_filtered, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples (Outliers Removed)",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels

# Plot expression levels of selected genes for Tumor stage
ggplot(expression_filtered, aes(x = Tumor_stage, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples (Outliers Removed)",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))# Rotate sample labels

# Rename NA into control in expression_filtered
expression_filtered$Tumor_stage[is.na(expression_filtered$Tumor_stage)] <- "Control"

# Boxplot of gene expression by tumor stage
ggplot(expression_filtered, aes(x = Tumor_stage, y = Expression, fill = Tumor_stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Tumor stage",
       x = "Tumor stage",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability


# Alternatively, you can use a violin plot
ggplot(expression_filtered, aes(x = Tumor_stage, y = Expression, fill = Tumor_stage)) +
  geom_violin(alpha = 0.6) +  # Shows density of expression levels
  geom_jitter(width = 0.2, alpha = 0.7) +  # Adds individual sample points
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Tumor stage",
       x = "Tumor stage",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#=======================================#
# Differential Gene Expression Analysis #
#=======================================#
 
Counts <- read.delim("Raw_Counts_GSE81089.tsv", row.names=1, header=TRUE, sep="\t", check.names = FALSE)

# Deseq2
# making sure the row names in metadata matches to column names in Data
all(colnames(Counts) %in% rownames(metadata))

# Find Columns in Data That Are Not in metadata
setdiff(colnames(Counts), rownames(metadata))

# Remove everything after the underscore in Data
colnames(Counts) <- sub("_.*", "", colnames(Counts))

# Check again if row names in metadata matches to column names in Data
all(colnames(Counts) %in% rownames(metadata))

# Check if they are in the same order
all(colnames(Counts) == rownames(metadata))

# Reorder metadata rows to match the column order in Data
metadata <- metadata[match(colnames(Counts), rownames(metadata)), , drop = FALSE]

# Check if they now match
all(colnames(Counts) == rownames(metadata))

# Check the values in the data
summary(Counts)

# Convert all data values to Absolute values. (Non-negative)
Data <- abs(Counts)

# Round values to integers
Data <- round(Counts)

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = Counts,
  colData = metadata,
  design = ~ Source)

dds

# Quality control
# Remove genes with low counts
keep2 <- rowMeans(counts(dds)) >=10
dds <- dds[keep2,]

# set the factor level
dds$Source <- relevel(dds$Source, ref = "Human non-malignant tissue")

# Run DESeq
dds <- DESeq(dds)
res <- results(dds) # Differentially expressed genes are given

res 

# Explore results
summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# MA plot
plotMA(res0.01)


#=======================================#
# Individual Subquestion Tumor stage #
#=======================================#

# Redefine the genes of interest for subQ
genes_of_interest_2 <- c("ENSG00000204305", "ENSG00000107159", "ENSG00000129173")

# Reconvert Data to a long format (genes in rows, samples in columns)
expression_long_2 <- Data %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% genes_of_interest_2) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")


# Rename NA into control and convert it into a factor
metadata$Tumor_stage[is.na(metadata$Tumor_stage)] <- "Control"
metadata$Tumor_stage <- as.factor(metadata$Tumor_stage)

# Construct a DESeqDataSet object 
dds2 <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = metadata,
                              design = ~ Tumor_stage)

dds2

# Quality control
# Remove genes with low counts
keep2 <- rowMeans(counts(dds2)) >=10
dds2 <- dds2[keep2,]

# set the factor level
dds2$Tumor_stage <- as.factor(dds2$Tumor_stage)
dds2$Tumor_stage <- relevel(dds2$Tumor_stage, ref = "Control")

# Run DESeq
dds2 <- DESeq(dds2)
res2 <- results(dds2) # Differentially expressed genes are given

res2

# Explore results
summary(res2_0.01)
res2_0.01 <- results(dds2, alpha = 0.01)
summary(res2_0.01)

# MA plot
plotMA(res2_0.01)

# Export DEG's to a tsv file
res2_df <- data.frame(Gene = rownames(res2_0.01), res2_0.01, row.names = NULL)  # Move row names to first column
write.table(res2_df, file = "DEG's_CancervsNon_Cancer_Cells.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
