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


#E-MTAB-12862 loading and exploratory analysis

#loading data, please change paths to your local directories
MTAB_counts <- read.table(gzfile("/path_to_file"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
metadata_MTAB<- read_xlsx("/path_to_file")
MTAB_counts_symbol <- read.table(gzfile("/path_to_file"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

MTAB_counts_filtered <- MTAB_counts
colnames(MTAB_counts_filtered) <- gsub("\\.", "-", colnames(MTAB_counts_filtered))
row.names(MTAB_counts_filtered) <- MTAB_counts_filtered$`ENSG-ID`
colnames_cleaned <- gsub("CRC-SW-(U[0-9]+)-[TN]", "\\1", colnames(MTAB_counts_filtered))
colnames(MTAB_counts_filtered) <- colnames_cleaned


metadata_MTAB_filtered <- metadata_MTAB[metadata_MTAB$`Sample ID` %in% colnames_cleaned, ]
row.names(metadata_MTAB_filtered) <- metadata_MTAB_filtered$`Sample ID`


MTAB_counts_filtered <- MTAB_counts_filtered[, metadata_MTAB_filtered$`Sample ID`]

dds_MTAB <- DESeqDataSetFromMatrix(countData = MTAB_counts_filtered, 
                                   colData = metadata_MTAB_filtered, 
                                   design = ~1)


vsd_MTAB <- vst(dds_MTAB, blind = FALSE)
normalized_counts_MTAB <- assay(vsd_MTAB)

variable_genes_MTAB <- apply(normalized_counts_MTAB, 1, var) > 0
normalized_counts_filtered_MTAB <- normalized_counts_MTAB[variable_genes_MTAB, ]

normalized_counts_filtered_MTAB_TRANSPOSED <- t(normalized_counts_filtered_MTAB)
pca_result_mtab <- prcomp(normalized_counts_filtered_MTAB_TRANSPOSED, center = TRUE, scale. = TRUE)


pca_mtab_df <- data.frame(
  Sample = rownames(normalized_counts_filtered_MTAB_TRANSPOSED),
  PC1 = pca_result_mtab$x[,3],#you can change the PCs
  PC2 = pca_result_mtab$x[,4],
  Site = metadata_MTAB_filtered$`Anatomic Organ Subdivision`,
  MSI_status = metadata_MTAB_filtered$`MSI Status`,
  Gender = metadata_MTAB_filtered$Sex,
  CMS = metadata_MTAB_filtered$`CMS Tumour`,
  CRPS = metadata_MTAB_filtered$`CRPS Tumour`,
  site_3 = metadata_MTAB_filtered$`Tumour Site`,
  stage = metadata_MTAB_filtered$`Tumour Stage`
)
pca_mtab_df <- pca_mtab_df %>% filter(!Site == "Rectum")





ggplot(mapping = aes(x = PC1, y = PC2, color = MSI_status, fill = Site)) +
  geom_point(data = pca_mtab_df[pca_mtab_df$Gender == "Female", ], 
             aes(x = PC1, y = PC2, color = MSI_status, fill = Site), 
             size = 3, alpha = 0.7, shape = 21) + 
  geom_point(data = pca_mtab_df[pca_mtab_df$Gender == "Male", ], 
             aes(x = PC1, y = PC2, color = MSI_status, fill = Site), 
             size = 3, alpha = 0.7, shape = 25) + 
  labs(title = "PCA Analysis of MTAB",
       x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


explained_variance <- summary(pca_result_mtab)$importance[2, 1:4] * 100

pc1_pc2 <- ggplot(pca_mtab_df, aes(x = PC1, y = PC2, color = site_3, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis in E-MTAB-12862",
       x = paste0("PC3 (", round(explained_variance[3],1), "%)"), y = paste0("PC4 (", round(explained_variance[4],1), "%)")) +
  theme_minimal()
# ggsave(filename = "PC3_PC4_MTAB.pdf", pc1_pc2, height = 6, width = 8)

#Site PCA
ggplot(pca_mtab_df, aes(x = PC1, y = PC2, color = site_3, shape = Gender)) +
  geom_point(size = 3) +
  labs(title = "PCA Analysis of MTAB",
       x = "PC3", y = "PC4") +
  theme_minimal()




