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

#PAMr classifier
#please execute all previous scripts to get all the files

library(pamr)
library(caret)


metadata_MTAB_filtered_PAMR <- metadata_MTAB_filtered %>% filter(!metadata_MTAB_filtered$`Tumour Site` == "Rectum")
normalized_counts_filtered_MTAB_PAMR <- normalized_counts_filtered_MTAB[, colnames(normalized_counts_filtered_MTAB) %in% metadata_MTAB_filtered_PAMR$`Sample ID`]
normalized_counts_filtered_MTAB_PAMR_df <- as.data.frame(normalized_counts_filtered_MTAB_PAMR)

#same genes in tcga and mtab
normalized_counts_filtered_MTAB_PAMR_df$genes <- rownames(normalized_counts_filtered_MTAB_PAMR_df)

normalized_PAM_TCGA_classifier <- sub("\\..*", "", rownames(counts_wt_vs_V600E))
common_TCGA_MTAB_genes <- intersect(normalized_counts_filtered_MTAB_PAMR_df$genes, normalized_PAM_TCGA_classifier)
normalized_counts_filtered_MTAB_PAMR_df <- normalized_counts_filtered_MTAB_PAMR_df %>% filter(normalized_counts_filtered_MTAB_PAMR_df$genes %in% common_TCGA_MTAB_genes)

normalized_counts_filtered_MTAB_PAMR_df <- normalized_counts_filtered_MTAB_PAMR_df %>% dplyr::select(-genes)
normalized_counts_filtered_MTAB_PAMR_genes_filtered <- as.matrix(normalized_counts_filtered_MTAB_PAMR_df)


X <- t(normalized_counts_filtered_MTAB_PAMR_genes_filtered)#Transpose to have samples as rows
Y <- as.factor(metadata_MTAB_filtered_PAMR$`Tumour Site`)

rownames(X) <- metadata_MTAB_filtered_PAMR$`Sample ID`
table(Y)

#training model with E-MTAB-12862 expression data
pam_data <- list(
  x = t(X), 
  y = Y, 
  geneid = rownames(normalized_counts_filtered_MTAB_PAMR_df)
)

pam_model <- pamr.train(pam_data)

#cv plots
pam_cv <- pamr.cv(pam_model, pam_data)
plotcv <- pamr.plotcv(pam_cv)

#threshold
optimal_threshold <- pam_cv$threshold[which.min(pam_cv$error)]

important_genes <- pamr.listgenes(pam_model, pam_data, threshold = optimal_threshold)
imporatnt_genes_df <- as.data.frame(important_genes)

# write.csv(imporatnt_genes_df, "important_genes_classifier.csv")
# 
# save(pam_model, optimal_threshold, pam_data, file = "PAM_MTAB_model.RData")

#predicting data in TCGA
dds_PAM_TCGA_classifier <- DESeqDataSetFromMatrix(countData = counts_wt_vs_V600E, 
                                                  colData = merged_data_wt_vs_v600e, 
                                                  design = ~ condition)


vsd_PAM_TCGA_classifier <- vst(dds_PAM_TCGA_classifier, blind = FALSE)
normalized_PAM_TCGA_classifier <- assay(vsd_PAM_TCGA_classifier)

rownames(normalized_PAM_TCGA_classifier) <- gsub("\\..*", "", rownames(normalized_PAM_TCGA_classifier))

common_genes <- intersect(rownames(pam_data$x), rownames(normalized_PAM_TCGA_classifier))

#subset both datasets to only include common genes, both datasets need to have the same genes
new_X_filtered <- normalized_PAM_TCGA_classifier[common_genes, , drop = FALSE]
pam_data$x <- pam_data$x[common_genes, , drop = FALSE]

#predicting data in TCGA
new_X_filtered <- new_X_filtered[rownames(pam_data$x), , drop = FALSE]
print(identical(rownames(pam_data$x), rownames(new_X_filtered)))


predictions <- pamr.predict(pam_model, new_X_filtered, threshold = optimal_threshold, type = "class")

predictions_df <- data.frame(Sample = colnames(new_X_filtered), Tumor_Site_Predicted = predictions)


merged_data_wt_vs_v600e_with_classifier <- merged_data_wt_vs_v600e
merged_data_wt_vs_v600e_with_classifier$SAMPLE <- rownames(merged_data_wt_vs_v600e)
colnames(predictions_df)[1] <- "SAMPLE"
merged_data_wt_vs_v600e_with_classifier <- merged_data_wt_vs_v600e_with_classifier %>% inner_join(predictions_df, by = "SAMPLE")

posterior <- pamr.predict(pam_model, new_X_filtered, threshold = optimal_threshold, type="posterior")

posterior <- as.data.frame(posterior)
mean(posterior$`Right Colon`)
mean(posterior$`Left Colon`) 

merged_data_wt_vs_v600e_with_classifier$Tumor_Site_Predicted <- factor(
  ifelse(merged_data_wt_vs_v600e_with_classifier$Tumor_Site_Predicted == "Right Colon", "right", "left")
)

# both columns are factors with the same levels
merged_data_wt_vs_v600e_with_classifier$Site <- factor(merged_data_wt_vs_v600e_with_classifier$Site)


levels(merged_data_wt_vs_v600e_with_classifier$Tumor_Site_Predicted)
levels(merged_data_wt_vs_v600e_with_classifier$Site)

merged_data_wt_vs_v600e_with_classifier_filtered <- merged_data_wt_vs_v600e_with_classifier %>% filter(merged_data_wt_vs_v600e_with_classifier$Site == "right" | merged_data_wt_vs_v600e_with_classifier$Site == "left")
merged_data_wt_vs_v600e_with_classifier_filtered$Site <- factor(merged_data_wt_vs_v600e_with_classifier_filtered$Site)

merged_data_wt_vs_v600e_with_classifier_filtered$Tumor_Site_Predicted <- factor(merged_data_wt_vs_v600e_with_classifier_filtered$Tumor_Site_Predicted)

# confusionMatrix
caret::confusionMatrix(
  merged_data_wt_vs_v600e_with_classifier_filtered$Tumor_Site_Predicted,
  merged_data_wt_vs_v600e_with_classifier_filtered$Site
)

#preparing the data for differential expression analysis
#final tcga metadata

missing_diff_exp_TCGA_metadata <- merged_data_wt_vs_v600e_with_classifier %>% filter(!(merged_data_wt_vs_v600e_with_classifier$Site == "right" | merged_data_wt_vs_v600e_with_classifier$Site == "left"))
missing_diff_exp_TCGA_metadata <- missing_diff_exp_TCGA_metadata %>% dplyr::select(-Site)
colnames(missing_diff_exp_TCGA_metadata)[9] <- "Site"

diff_exp_TCGA_metadata <- merged_data_wt_vs_v600e_with_classifier %>% dplyr::select(-Tumor_Site_Predicted)
diff_exp_TCGA_metadata <- diff_exp_TCGA_metadata %>% filter(!(diff_exp_TCGA_metadata$SAMPLE %in% missing_diff_exp_TCGA_metadata$SAMPLE))
diff_exp_TCGA_metadata <- rbind(diff_exp_TCGA_metadata, missing_diff_exp_TCGA_metadata)


diff_exp_TCGA_metadata$Site[diff_exp_TCGA_metadata$SAMPLE == "TCGA-5M-AAT5-01A-21R-A41B-07"] <- "right"
final_diff_exp_TCGA_metadata <- diff_exp_TCGA_metadata
rownames(final_diff_exp_TCGA_metadata) <- final_diff_exp_TCGA_metadata$SAMPLE

final_diff_exp_TCGA_metadata <- final_diff_exp_TCGA_metadata %>% dplyr::select(-SAMPLE)

#final counts TCGA
final_diff_exp_TCGA_counts <- counts_wt_vs_V600E




