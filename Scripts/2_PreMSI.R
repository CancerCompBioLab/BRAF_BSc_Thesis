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
#please execute all the previous code
library(PreMSIm)
#Please load exp_COAD_TPM from 0_TCGAbiolinks.R
expression_unstranded_data_tpm <- exp_COAD_TPM

df_expression_unstranded_data_tpm <- data.frame(row.names = row.names(expression_unstranded_data_tpm),
                                                expression_unstranded_data_tpm)
#I take off the .X after the ensemble IDs
df_expression_unstranded_data_tpm$Gene_ID <- sub("\\..*", "", rownames(df_expression_unstranded_data_tpm))

df_cleaned <- df_expression_unstranded_data_tpm %>% group_by(Gene_ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% ungroup()

#Translating the ENSMBL ID to Gene Symbols
gtf_file <- "/path_to_file"#the file is too big to upload in github, you should download it from "https://www.gencodegenes.org/human/release_47.html"
genes_gtf <- rtracklayer::import(gtf_file)
genes <- genes_gtf[genes_gtf$type == "gene"]

gene_info <- data.frame(
  gene_id = genes$gene_id,
  gene_name = genes$gene_name
)

gene_info$gene_id <- gsub("\\..*", "", gene_info$gene_id)
df_cleaned <- merge(df_cleaned, gene_info, by.x = "Gene_ID", by.y = "gene_id", all.x = TRUE)
df_cleaned$gene_name[is.na(df_cleaned$gene_name)] <- df_cleaned$Gene_ID[is.na(df_cleaned$gene_name)]
df_cleaned$gene_name <- make.names(df_cleaned$gene_name, unique = TRUE)
row.names(df_cleaned) <- df_cleaned$gene_name
data_cleaned_MSI <- df_cleaned %>% dplyr::select(-Gene_ID, -gene_name)

# write.table(data_cleaned_MSI, "df_expression_unstranded_data_tpm.txt", sep = "\t")

path = system.file("extdata", "df_expression_unstranded_data_tpm.txt", package = "PreMSIm", mustWork = TRUE)
input_data = data_pre(path, type = "Symbol")
results_PreMSI <- msi_pre(input_data)
results_PreMSI$Sample <- gsub("\\.", "-", results_PreMSI$Sample)


merged_data_wt_vs_v600e <- merge(results_PreMSI, final_metadata_wt_vs_v600e, 
                                 by.x = "Sample", by.y = "sample", 
                                 all.x = TRUE) 
merged_data_wt_vs_v600e <- merged_data_wt_vs_v600e %>% filter(Sample %in% final_metadata_wt_vs_v600e$sample)


# 
# 1 Indicates MSI-High these samples have high levels of microsatellite
# instability, which is often associated with defective DNA mismatch
# repair (MMR). 
# 0 Indicates MSI-Low these samples have low or no
# microsatellite instability, meaning their DNA mismatch repair system is
# functioning normally.

merged_data_wt_vs_v600e <- merged_data_wt_vs_v600e %>%
  mutate(MsiStatus = ifelse(is.na(MsiStatus), 
                            ifelse(MSI_status == 1, "MSI-H", "MSS"), 
                            MsiStatus)) %>% 
  mutate(MsiStatus = ifelse(MsiStatus %in% c("POLE", "MSI-H/POLE"), 
                            ifelse(MSI_status == 1, "MSI-H", "MSS"), 
                            MsiStatus))
merged_data_wt_vs_v600e <- merged_data_wt_vs_v600e %>% dplyr::select(-MSI_status)
row.names(merged_data_wt_vs_v600e) <- merged_data_wt_vs_v600e$Sample
merged_data_wt_vs_v600e <- merged_data_wt_vs_v600e[, -1]






