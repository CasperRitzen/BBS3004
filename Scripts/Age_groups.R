# this is a script for RNA-seq Analysis
# libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(GEOquery)      

# Set your working directory
setwd("/Users/oguzhanbayrak/Desktop/BBS3004/")

# step 1 of data extraction
Data <- read.delim("GSE81089_FPKM_cufflinks.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
colnames(metadata)
metadata <- metadata %>%
  #write.table(metadata, file ="metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
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
genes_of_interest <- c("ENSG00000146648", "ENSG00000133703", "ENSG00000157764")

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
                              design = ~ Age)

# Convert Age to numeric (if it's not already)
metadata$Age_Group <- as.numeric(metadata$Age_Group)

# Define Age Groups
metadata$Age_Group <- cut(metadata$Age,
                          breaks = c(45, 55, 65, 75, 85),
                          labels = c("45-55", "56-65", "66-75", "76-85"),
                          include.lowest = TRUE)

# Check the changes
table(metadata$Age_Group)  # Verify that Age Groups are correctly assigned

# Print the first few rows to confirm
head(metadata)

table(metadata$Age_Group)  # Count number of samples per age group
ggplot(metadata, aes(x = Age_Group, fill = Age_Group)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Distribution of Samples Across Age Groups", x = "Age Group", y = "Count")

# Step 1: Remove samples with missing Age_Group
metadata <- metadata %>% filter(!is.na(Age_Group))

# Step 2: Filter Data to match the samples in metadata
Data_filtered <- Data[, colnames(Data) %in% rownames(metadata)]  # Filter Data accordingly

# Step 3: Check the dimensions to ensure that the number of samples in Data matches the number of samples in metadata
print(dim(Data_filtered))  # Check dimensions of Data
print(dim(metadata))       # Check dimensions of metadata

# Step 4: Ensure that the column names of Data match the row names of metadata
if(!all(colnames(Data_filtered) == rownames(metadata))) {
  # If they don't match, reorder metadata to match the columns of Data
  metadata <- metadata[match(colnames(Data_filtered), rownames(metadata)), , drop = FALSE]
}

# Step 5: Confirm that the row names of metadata and column names of Data match
if(all(colnames(Data_filtered) == rownames(metadata))) {
  cat("Data and metadata match!\n")
} else {
  stop("Data and metadata do not match. Please check the filtering process.")
}

# Step 6: Run DESeq2 to create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = Data_filtered,
                              colData = metadata,
                              design = ~ Age_Group)

# Step 7: Run the DESeq2 analysis
dds <- DESeq(dds)

dds <- DESeq(dds)

# Correspondingly filter the Data matrix to match the filtered metadata
Data_filtered <- Data[, colnames(Data) %in% rownames(metadata)]

# Check dimensions to ensure Data and metadata are aligned
dim(Data_filtered)
dim(metadata)






dds <- DESeqDataSetFromMatrix(countData = Data,
                              colData = metadata,
                              design = ~ Age_Group)
dds <- DESeq(dds)

# Extract results for specific age group comparisons (e.g., 56-65 vs 45-55)
res <- results(dds, contrast = c("Age_Group", "56-65", "45-55"))
summary(res)

ggplot(expression_long, aes(x = Age_Group, y = Expression, fill = Age_Group)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Gene Expression Levels by Age Group",
       x = "Age Group", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##################
#new try
