# this is a script for RNA-seq Analysis
# libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(GEOquery)      
# Set your working directory
setwd("C:/Users/ritze/OneDrive/Documenten/Maastricht_University/BBS3004Data/FPKM")
#=======================================#
# Step 1. Import Data & metadata into R #
#=======================================#

# Load your Data
Data <- read.delim("GSE81089_FPKM_cufflinks.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)

# get metadata 
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)  # Have a glimpse of how metadata looks like

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

genes_of_interest <- c("ENSG00000087494", "ENSG00000133048", "ENSG00000001084")

#PTHLH =ENSG00000087494
#ENSG00000133048=CHI3L1
#ENSG00000001084=GCLC

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

#recoding smoking status
expression_long$Smoking_Status <- recode(expression_long$Smoking_Status, 
                                  "1" = "Current", 
                                  "2" = "ex > 1 year", 
                                  "3" = "Never") 
#turning NA into Control
expression_long$Smoking_Status[is.na(expression_long$Smoking_Status)] <- "Control"


# Compute IQR for each gene and filter out outliers
expression_filtered <- expression_long %>%
  group_by(Gene) %>%  # Ensure outliers are removed per gene
  mutate(Q1 = quantile(Expression, 0.25),
         Q3 = quantile(Expression, 0.75),
         IQR = Q3 - Q1) %>%
  filter(Expression >= (Q1 - 1.5 * IQR) & Expression <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)  # Remove extra columns


# Compute mean expression per Gene and Smoking_Status
expression_avg <- expression_filtered %>%
  group_by(Smoking_Status, Gene) %>%
  summarize(Average_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Plot with averaged values
ggplot(expression_avg, aes(x = Smoking_Status, y = Average_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Average_Expression, 1)),  # Show rounded average values
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 4) +  
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Average Gene Expression Levels Across Samples",
       x = "Smoking Status",
       y = "Average Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels

# Boxplot of gene expression by Smoking Status
ggplot(expression_filtered, aes(x = Smoking_Status, y = Expression, fill = Smoking_Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Smoking Status",
       x = "Smoking_Status",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability


# Alternatively, you can use a violin plot
ggplot(expression_filtered, aes(x = Smoking_Status, y = Expression, fill = Smoking_Status)) +
  geom_violin(alpha = 0.6) +  # Shows density of expression levels
  geom_jitter(width = 0.2, alpha = 0.7) +  # Adds individual sample points
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Smoking Status",
       x = "Smoking status",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#=======================================#
# Differential Gene Expression Analysis #
#=======================================#
# Load your Data
Data <- read.delim("Raw_Counts_GSE81089.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)

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

#recoding smoking status
metadata$Smoking_Status <- recode(metadata$Smoking_Status, 
                                  "1" = "Current", 
                                  "2" = "ex > 1 year", 
                                  "3" = "Never") 

#turning NA into Control
metadata$Smoking_Status[is.na(metadata$Smoking_Status)] <- "Control"

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = Data,
                              colData = metadata,
                              design = ~ Smoking_Status)


# Quality control
# Remove genes with low counts

keep <- rowMeans(counts(dds)) >=10
dds <- dds[keep,]
#set factor level
dds$Smoking_Status <- relevel(dds$Smoking_Status, ref = "Never")

#run dds
dds <- DESeq(dds)
res <- results(dds)
res
#Explore results
summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)
#MA plot
plotMA(res)
plotMA(res0.01)
plotMA(res0.05)
# Export DEG results to a file
res_df <- data.frame(Gene = rownames(res0.01), res0.01, row.names = NULL)  # Move row names to first column
write.table(res_df, file = "DEG's_Smoking_Status.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

res0.01
