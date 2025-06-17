library(readr)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(maftools)
library(readxl)
library(biomaRt)
library(rtracklayer)
library(reshape2)

#Exploratory analysis TCGA 

merged_data_wt_vs_v600e <- read_csv("Datasets/merged_data_wt_vs_v600e.csv") #load the merged_data_wt_vs_v600e if you could execute the previous script. The data is in the github repository folder "Datsets".
counts_wt_vs_V600E <- cbind(filtered_wtBRAF_expression,
                            filtered_V600E_expression,
                            filtered_nonV600E_expression)
merged_data_wt_vs_v600e <- merged_data_wt_vs_v600e[colnames(counts_wt_vs_V600E), , drop=FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts_wt_vs_V600E, 
                              colData = merged_data_wt_vs_v600e, 
                              design = ~ condition)
#vst normalization
vsd <- vst(dds, blind = FALSE)
normalized_counts <- assay(vsd)

variable_genes <- apply(normalized_counts, 1, var) > 0
normalized_counts_filtered <- normalized_counts[variable_genes, ]

normalized_counts_filtered_transposed <- t(normalized_counts_filtered)
pca_result <- prcomp(normalized_counts_filtered_transposed, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(normalized_counts_filtered_transposed),
  PC1 = pca_result$x[,3],#you can put 1 and 2 if you want to check the other principal components
  PC2 = pca_result$x[,4],
  Site = merged_data_wt_vs_v600e$Site,
  MSI_status = merged_data_wt_vs_v600e$MsiStatus,
  Gender = merged_data_wt_vs_v600e$gender,
  TumorPurity = merged_data_wt_vs_v600e$TumorPurity,
  condition = merged_data_wt_vs_v600e$condition,
  Site_2 = merged_data_wt_vs_v600e$anatomic_neoplasm_subdivision
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis of Tumor Site (VST Normalized)",
       x = "PC1", y = "PC2") +
  theme_minimal()

pc1_mean <- mean(pca_df$PC1)
pc1_sd <- sd(pca_df$PC1)
pc2_mean <- mean(pca_df$PC2)
pc2_sd <- sd(pca_df$PC2)

# threshold for outliers 
threshold <- 2

#filter out samples that are extreme outliers
pca_df_filtered <- pca_df[
  abs(pca_df$PC1 - pc1_mean) < threshold * pc1_sd & 
    abs(pca_df$PC2 - pc2_mean) < threshold * pc2_sd, 
]

ggplot(pca_df_filtered, aes(x = PC1, y = PC2, color = MSI_status, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis TCGA",
       x = "PC1", y = "PC2") +
  theme_minimal()

ggplot(pca_df_filtered, aes(x = PC1, y = PC2, color = Site, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis TCGA",
       x = "PC1", y = "PC2") +
  theme_minimal()



ggplot(pca_df_filtered, aes(x = PC1, y = PC2, color = Site_2, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis TCGA",
       x = "PC1", y = "PC2") +
  theme_minimal()

female_TCGA <- pca_df_filtered %>% filter(pca_df_filtered$Gender == "FEMALE")
ggplot(female_TCGA, aes(x = PC1, y = PC2, color = Site_2, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis TCGA FEMALES",
       x = "PC1", y = "PC2") +
  theme_minimal()

male_TCGA <- pca_df_filtered %>% filter(pca_df_filtered$Gender == "MALE")
ggplot(male_TCGA, aes(x = PC1, y = PC2, color = Site_2, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis TCGA MALES",
       x = "PC1", y = "PC2") +
  theme_minimal()




