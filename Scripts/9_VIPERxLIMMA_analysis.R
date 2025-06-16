#In this script I perform the VIPER analysis (with limma voom).

#V600E vs WT/nonV600E

metadata_MTAB_v600e_vs_wt_nonv600e <- metadata_MTAB_mutation %>% filter(!metadata_MTAB_mutation$`Anatomic Organ Subdivision` == "Rectum")
metadata_MTAB_v600e_vs_wt_nonv600e <- metadata_MTAB_v600e_vs_wt_nonv600e[,-1]
vector_mut <- c("p.Val640Glu")
metadata_MTAB_v600e_vs_wt_nonv600e$mutation[!(metadata_MTAB_v600e_vs_wt_nonv600e$mutation %in% vector_mut)] <- "nonV600E_wt"

MTAB_counts_v600e_vs_wt_nonv600e <- MTAB_counts
colnames(MTAB_counts_v600e_vs_wt_nonv600e) <- gsub("\\.", "-", colnames(MTAB_counts_v600e_vs_wt_nonv600e))
rownames(MTAB_counts_v600e_vs_wt_nonv600e) <- MTAB_counts_v600e_vs_wt_nonv600e$`ENSG-ID`
MTAB_counts_v600e_vs_wt_nonv600e <- MTAB_counts_v600e_vs_wt_nonv600e[, -1]

MTAB_counts_v600e_vs_wt_nonv600e <- MTAB_counts_v600e_vs_wt_nonv600e[, colnames(MTAB_counts_v600e_vs_wt_nonv600e) %in% metadata_MTAB_v600e_vs_wt_nonv600e$DNA_Tumor_Sample_Barcode]

metadata_MTAB_v600e_vs_wt_nonv600e <- metadata_MTAB_v600e_vs_wt_nonv600e[
  metadata_MTAB_v600e_vs_wt_nonv600e$DNA_Tumor_Sample_Barcode %in% colnames(MTAB_counts_v600e_vs_wt_nonv600e),]

colnames(metadata_MTAB_v600e_vs_wt_nonv600e)[40] <- "MSI_status"
colnames(metadata_MTAB_v600e_vs_wt_nonv600e)[33] <- "Tumour_Cell_Content_Pathology"
colnames(metadata_MTAB_v600e_vs_wt_nonv600e)[17] <- "Age_at_diagnosis"
colnames(metadata_MTAB_v600e_vs_wt_nonv600e)[21] <- "Primary_Site_Disease"
colnames(metadata_MTAB_v600e_vs_wt_nonv600e)[23] <- "Tumour_Site"

dge_V600E_vs_WT.nonV600E <- DGEList(counts = MTAB_counts_v600e_vs_wt_nonv600e)

dge_V600E_vs_WT.nonV600E <- calcNormFactors(dge_V600E_vs_WT.nonV600E)

cutoff <- 18
keep_genes_V600E_vs_WT.nonV600E <- apply(cpm(dge_V600E_vs_WT.nonV600E), 1, max) >= cutoff
dge_V600E_vs_WT.nonV600E_filtered <- dge_V600E_vs_WT.nonV600E[keep_genes_V600E_vs_WT.nonV600E, ]

metadata_MTAB_v600e_vs_wt_nonv600e$mutation <- relevel(factor(metadata_MTAB_v600e_vs_wt_nonv600e$mutation), ref = "nonV600E_wt")

design_V600E_vs_WT.nonV600E <- model.matrix(~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation, data = metadata_MTAB_v600e_vs_wt_nonv600e)

v_V600E_vs_WT.nonV600E <- voom(dge_V600E_vs_WT.nonV600E_filtered, design_V600E_vs_WT.nonV600E, plot = TRUE)
fit_V600E_vs_WT.nonV600E <- lmFit(v_V600E_vs_WT.nonV600E, design_V600E_vs_WT.nonV600E)
fit_V600E_vs_WT.nonV600E <- eBayes(fit_V600E_vs_WT.nonV600E)

topTable_V600E_vs_WT.nonV600E <- topTable(fit_V600E_vs_WT.nonV600E, coef = "mutationp.Val640Glu", adjust.method = "BH", number = Inf)
topTable_V600E_vs_WT.nonV600E$ENSEMBL <- rownames(topTable_V600E_vs_WT.nonV600E)

# Load GENCODE v36 annotations
gencode_v36 <- rtracklayer::import("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.annotation.gtf")
entrez_annot_v36 <- data.table::fread("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.metadata.EntrezGene.gz",
                                      col.names = c("transcript_id", "entrez_id"))


temp_v36 <- gencode_v36 %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  dplyr::select(gene_id, transcript_id) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  merge(entrez_annot_v36, by = "transcript_id", all.x = TRUE)

gencode_v36_genes <- gencode_v36[gencode_v36$type == "gene"]
gencode_v36_genes$entrez_id <- temp_v36$entrez_id[match(gencode_v36_genes$gene_id, temp_v36$gene_id, nomatch = NA)]

# Merge with DE table to get gene names
annot_df <- as.data.frame(gencode_v36_genes)[, c("gene_id", "entrez_id", "gene_name")]
annot_df$gene_id <- sub("\\..*", "", annot_df$gene_id)
topTable_V600E_vs_WT.nonV600E$gene_name <- annot_df$gene_name[match(topTable_V600E_vs_WT.nonV600E$ENSEMBL, annot_df$gene_id)]

# Filter valid genes
topTable_V600E_vs_WT.nonV600E <- topTable_V600E_vs_WT.nonV600E[!is.na(topTable_V600E_vs_WT.nonV600E$gene_name), ]
topTable_V600E_vs_WT.nonV600E <- topTable_V600E_vs_WT.nonV600E[!duplicated(topTable_V600E_vs_WT.nonV600E$gene_name), ]

topTable_V600E_vs_WT.nonV600E$zscore <- qnorm(topTable_V600E_vs_WT.nonV600E$P.Value / 2, lower.tail = FALSE) * sign(topTable_V600E_vs_WT.nonV600E$t)
signature_V600E_vs_WT.nonV600E <- setNames(topTable_V600E_vs_WT.nonV600E$zscore, topTable_V600E_vs_WT.nonV600E$gene_name)

# set.seed(123)
# n_perm <- 1000
# genes <- rownames(v_V600E_vs_WT.nonV600E$E)
# null_matrix_V600E_vs_WT.nonV600E <- matrix(NA, nrow = length(genes), ncol = n_perm)
# rownames(null_matrix_V600E_vs_WT.nonV600E) <- genes
# 
# for (i in 1:n_perm) {
#   metadata_MTAB_v600e_vs_wt_nonv600e$mutation_perm <- sample(metadata_MTAB_v600e_vs_wt_nonv600e$mutation)
#   design_perm <- model.matrix(~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation_perm, 
#                               data = metadata_MTAB_v600e_vs_wt_nonv600e)
#   
#   fit_perm <- lmFit(v_V600E_vs_WT.nonV600E, design_perm)
#   fit_perm <- eBayes(fit_perm)
#   
#   null_matrix_V600E_vs_WT.nonV600E[, i] <- fit_perm$t[, "mutation_permp.Val640Glu"]
# }

# save(v_V600E_vs_WT.nonV600E, metadata_MTAB_v600e_vs_wt_nonv600e, file = "VIPER_V600E_vs_WT.nonV600E.RData")

load("Null_Matrix_V600E_vs_WT.nonV600E.RData")




gencode_v36 <- rtracklayer::import("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.annotation.gtf")

gencode_v36_genes <- gencode_v36[gencode_v36$type == "gene"]
gencode_df <- as.data.frame(gencode_v36_genes)[, c("gene_id", "gene_name")]
gencode_df <- gencode_df[!duplicated(gencode_df$gene_id), ]

# Map ENSEMBL IDs to gene names in null_matrix
ensembl_ids <- rownames(null_matrix_V600E_vs_WT.nonV600E)
gencode_df$gene_id <- sub("\\..*", "", gencode_df$gene_id)
gene_names <- gencode_df$gene_name[match(ensembl_ids, gencode_df$gene_id)]

#Filter and rename null_matrix
null_matrix_mapped <- null_matrix_V600E_vs_WT.nonV600E[!is.na(gene_names), ]
rownames(null_matrix_mapped) <- gene_names[!is.na(gene_names)]

#Remove duplicates
null_matrix_mapped_V600E_vs_WT.nonV600E <- null_matrix_mapped[!duplicated(rownames(null_matrix_mapped)), ]

df_residual <- min(fit_V600E_vs_WT.nonV600E$df.residual)

# Convert entire matrix of t-stats to two-sided p-values
p_null <- 2 * pt(-abs(null_matrix_mapped_V600E_vs_WT.nonV600E), df = df_residual)

# Convert to z-scores 
null_matrix_mapped_V600E_vs_WT.nonV600E_ZSCORE <- qnorm(p_null / 2, lower.tail = FALSE) * sign(null_matrix_mapped_V600E_vs_WT.nonV600E)


#MTAB regulon
network_MTAB <- read.delim("/CCBdata/projects/BRAF_regulon/Aracne_COAD_output_MTAB_protein_coding/network.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
network_MTAB <- network_MTAB[, 1:3]
write.table(network_MTAB, "network_MTAB.adj", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
Expression_Data_MTAB <- read_tsv("MTAB_ARACNE_samples.tsv")
Expression_Data_MTAB <- as.data.frame(Expression_Data_MTAB)
rownames(Expression_Data_MTAB) <- Expression_Data_MTAB[[1]]
Expression_Data_MTAB[[1]] <- NULL
Expression_Data_MTAB <- as.matrix(Expression_Data_MTAB) 



regulon_object <- aracne2regulon("network_MTAB.adj", Expression_Data_MTAB, verbose = FALSE)

# Run msVIPER
mrs_V600E_vs_WT.nonV600E <- msviper(signature_V600E_vs_WT.nonV600E, regulon_object, nullmodel = null_matrix_mapped_V600E_vs_WT.nonV600E_ZSCORE, verbose = TRUE)

mrs_V600E_vs_WT.nonV600E_df <- summary(mrs_V600E_vs_WT.nonV600E, 40)
mrs_V600E_vs_WT.nonV600E_df$TF <- rownames(mrs_V600E_vs_WT.nonV600E_df)
mrs_V600E_vs_WT.nonV600E_df <- mrs_V600E_vs_WT.nonV600E_df[order(mrs_V600E_vs_WT.nonV600E_df$NES), ]
mrs_V600E_vs_WT.nonV600E_df$TF <- factor(mrs_V600E_vs_WT.nonV600E_df$TF, levels = mrs_V600E_vs_WT.nonV600E_df$TF)

# Plot
VIPER_BRAFV600E_vs_nonV600E.WT <- ggplot(mrs_V600E_vs_WT.nonV600E_df, aes(x = TF, y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "NES") +
  theme_minimal(base_size = 12) +
  labs(title = "VIPER: BRAF V600E vs WT/nonV600E",
       x = "Regulator",
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")
VIPER_BRAFV600E_vs_nonV600E.WT
# ggsave("VIPER_BRAFV600E_vs_nonV600E.WT.pdf", VIPER_BRAFV600E_vs_nonV600E.WT)

#shadow anaylsis
mrshadow_V600E_vs_WT.nonV600E <- shadow(mrs_V600E_vs_WT.nonV600E, regulators = 25, verbose = FALSE)
summary(mrshadow_V600E_vs_WT.nonV600E)

shadowed_TF_V600E_vs_WT.nonV600E <- sapply(strsplit(summary(mrshadow_V600E_vs_WT.nonV600E)$Shadow.pairs, " -> "), `[`, 2)

all_TF_V600E_vs_WT.nonV600E <- summary(mrshadow_V600E_vs_WT.nonV600E)$msviper.results$Regulon

non_shadowed_TF_V600E_vs_WT.nonV600E <- setdiff(all_TF_V600E_vs_WT.nonV600E, shadowed_TF_V600E_vs_WT.nonV600E)

mrs_df_V600E_vs_WT.nonV600E <- summary(mrs_V600E_vs_WT.nonV600E, length(mrs_V600E_vs_WT.nonV600E$es$nes))
mrs_df_V600E_vs_WT.nonV600E$TF <- rownames(mrs_df_V600E_vs_WT.nonV600E)

# Filter only non-shadowed TFs
mrs_V600E_vs_WT.nonV600E_filtered_df <- mrs_df_V600E_vs_WT.nonV600E[mrs_df_V600E_vs_WT.nonV600E$TF %in% non_shadowed_TF_V600E_vs_WT.nonV600E, ]


mrs_V600E_vs_WT.nonV600E_filtered_df <- mrs_V600E_vs_WT.nonV600E_filtered_df[order(mrs_V600E_vs_WT.nonV600E_filtered_df$NES), ]
mrs_V600E_vs_WT.nonV600E_filtered_df$TF <- factor(mrs_V600E_vs_WT.nonV600E_filtered_df$TF, levels = mrs_V600E_vs_WT.nonV600E_filtered_df$TF)

library(ggplot2)

VIPER_BRAFV600E_vs_nonV600E.WT_shadow <- ggplot(mrs_V600E_vs_WT.nonV600E_filtered_df, aes(x = TF, y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "NES") +
  theme_minimal(base_size = 12) +
  labs(title = "VIPER (Non-Shadowed TFs): BRAF V600E vs WT/nonV600E",
       x = "Regulator (Non-Shadowed)",
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")
VIPER_BRAFV600E_vs_nonV600E.WT_shadow
# ggsave("VIPER_BRAFV600E_vs_nonV600E.WT_shadow.pdf", VIPER_BRAFV600E_vs_nonV600E.WT_shadow)



#nonV600E vs WT/V600E

dge_nonV600E_vs_WT.V600E <- DGEList(counts = MTAB_counts_nonv600e_vs_wt_v600e)

dge_nonV600E_vs_WT.V600E <- calcNormFactors(dge_nonV600E_vs_WT.V600E)

cutoff <- 18
keep_genes_nonV600E_vs_WT.V600E <- apply(cpm(dge_nonV600E_vs_WT.V600E), 1, max) >= cutoff
dge_nonV600E_vs_WT.V600E_filtered <- dge_nonV600E_vs_WT.V600E[keep_genes_nonV600E_vs_WT.V600E, ]

metadata_MTAB_nonv600e_vs_wt_v600e$mutation <- relevel(factor(metadata_MTAB_nonv600e_vs_wt_v600e$mutation), ref = "V600E_wt")

design_nonV600E_vs_WT.V600E <- model.matrix(~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation, data = metadata_MTAB_nonv600e_vs_wt_v600e)

v_nonV600E_vs_WT.V600E <- voom(dge_nonV600E_vs_WT.V600E_filtered, design_nonV600E_vs_WT.V600E, plot = TRUE)
fit_nonV600E_vs_WT.V600E <- lmFit(v_nonV600E_vs_WT.V600E, design_nonV600E_vs_WT.V600E)
fit_nonV600E_vs_WT.V600E <- eBayes(fit_nonV600E_vs_WT.V600E)

topTable_nonV600E_vs_WT.V600E <- topTable(fit_nonV600E_vs_WT.V600E, coef = "mutationnonV600E", adjust.method = "BH", number = Inf)
topTable_nonV600E_vs_WT.V600E$ENSEMBL <- rownames(topTable_nonV600E_vs_WT.V600E)

# Load GENCODE v36 annotations
gencode_v36 <- rtracklayer::import("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.annotation.gtf")
entrez_annot_v36 <- data.table::fread("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.metadata.EntrezGene.gz",
                                      col.names = c("transcript_id", "entrez_id"))


temp_v36 <- gencode_v36 %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  dplyr::select(gene_id, transcript_id) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  merge(entrez_annot_v36, by = "transcript_id", all.x = TRUE)

gencode_v36_genes <- gencode_v36[gencode_v36$type == "gene"]
gencode_v36_genes$entrez_id <- temp_v36$entrez_id[match(gencode_v36_genes$gene_id, temp_v36$gene_id, nomatch = NA)]


annot_df <- as.data.frame(gencode_v36_genes)[, c("gene_id", "entrez_id", "gene_name")]
annot_df$gene_id <- sub("\\..*", "", annot_df$gene_id)
topTable_nonV600E_vs_WT.V600E$gene_name <- annot_df$gene_name[match(topTable_nonV600E_vs_WT.V600E$ENSEMBL, annot_df$gene_id)]

# Filter valid genes
topTable_nonV600E_vs_WT.V600E <- topTable_nonV600E_vs_WT.V600E[!is.na(topTable_nonV600E_vs_WT.V600E$gene_name), ]
topTable_nonV600E_vs_WT.V600E <- topTable_nonV600E_vs_WT.V600E[!duplicated(topTable_nonV600E_vs_WT.V600E$gene_name), ]


topTable_nonV600E_vs_WT.V600E$zscore <- qnorm(topTable_nonV600E_vs_WT.V600E$P.Value / 2, lower.tail = FALSE) * sign(topTable_nonV600E_vs_WT.V600E$t)
signature_nonV600E_vs_WT.V600E <- setNames(topTable_nonV600E_vs_WT.V600E$zscore, topTable_nonV600E_vs_WT.V600E$gene_name)

# set.seed(123)
# n_perm <- 1000
# genes <- rownames(v_nonV600E_vs_WT.V600E$E)
# null_matrix_nonV600E_vs_WT.V600E <- matrix(NA, nrow = length(genes), ncol = n_perm)
# rownames(null_matrix_nonV600E_vs_WT.V600E) <- genes
# 
# for (i in 1:n_perm) {
#   metadata_MTAB_nonv600e_vs_wt_v600e$mutation_perm <- sample(metadata_MTAB_nonv600e_vs_wt_v600e$mutation)
#   design_perm <- model.matrix(~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation_perm,
#                               data = metadata_MTAB_nonv600e_vs_wt_v600e)
# 
#   fit_perm <- lmFit(v_nonV600E_vs_WT.V600E, design_perm)
#   fit_perm <- eBayes(fit_perm)
# 
#   null_matrix_nonV600E_vs_WT.V600E[, i] <- fit_perm$t[, "mutation_permnonV600E"]
# }
# 
# save(v_nonV600E_vs_WT.V600E, metadata_MTAB_nonv600e_vs_wt_v600e, file = "VIPER_nonV600E_vs_WT.V600E.RData")

load("Null_Matrix_nonV600E_vs_WT.V600E.RData")


gencode_v36 <- rtracklayer::import("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.annotation.gtf")

gencode_v36_genes <- gencode_v36[gencode_v36$type == "gene"]
gencode_df <- as.data.frame(gencode_v36_genes)[, c("gene_id", "gene_name")]
gencode_df <- gencode_df[!duplicated(gencode_df$gene_id), ]

ensembl_ids <- rownames(null_matrix_nonV600E_vs_WT.V600E)
gencode_df$gene_id <- sub("\\..*", "", gencode_df$gene_id)
gene_names <- gencode_df$gene_name[match(ensembl_ids, gencode_df$gene_id)]


null_matrix_mapped <- null_matrix_nonV600E_vs_WT.V600E[!is.na(gene_names), ]
rownames(null_matrix_mapped) <- gene_names[!is.na(gene_names)]


null_matrix_mapped_nonV600E_vs_WT.V600E <- null_matrix_mapped[!duplicated(rownames(null_matrix_mapped)), ]

df_residual <- min(fit_nonV600E_vs_WT.V600E$df.residual)

#convert entire matrix of t-stats to two-sided p-values
p_null <- 2 * pt(-abs(null_matrix_mapped_nonV600E_vs_WT.V600E), df = df_residual)

#convert to z-scores
null_matrix_mapped_nonV600E_vs_WT.V600E_ZSCORE <- qnorm(p_null / 2, lower.tail = FALSE) * sign(null_matrix_mapped_nonV600E_vs_WT.V600E)

network_MTAB <- read.delim("/CCBdata/projects/BRAF_regulon/Aracne_COAD_output_MTAB_protein_coding/network.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
network_MTAB <- network_MTAB[, 1:3]
write.table(network_MTAB, "network_MTAB.adj", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
Expression_Data_MTAB <- read_tsv("MTAB_ARACNE_samples.tsv")
Expression_Data_MTAB <- as.data.frame(Expression_Data_MTAB)
rownames(Expression_Data_MTAB) <- Expression_Data_MTAB[[1]]
Expression_Data_MTAB[[1]] <- NULL
Expression_Data_MTAB <- as.matrix(Expression_Data_MTAB)

# regulon_object <- aracne2regulon("network_MTAB.adj", Expression_Data_MTAB, verbose = FALSE)

# Run msVIPER
mrs_nonV600E_vs_WT.V600E <- msviper(signature_nonV600E_vs_WT.V600E, regulon_object, nullmodel = null_matrix_mapped_nonV600E_vs_WT.V600E_ZSCORE, verbose = TRUE)


mrs_nonV600E_vs_WT.V600E_df <- summary(mrs_nonV600E_vs_WT.V600E, 50)
mrs_nonV600E_vs_WT.V600E_df$TF <- rownames(mrs_nonV600E_vs_WT.V600E_df)
mrs_nonV600E_vs_WT.V600E_df <- mrs_nonV600E_vs_WT.V600E_df[order(mrs_nonV600E_vs_WT.V600E_df$NES), ]
mrs_nonV600E_vs_WT.V600E_df$TF <- factor(mrs_nonV600E_vs_WT.V600E_df$TF, levels = mrs_nonV600E_vs_WT.V600E_df$TF)


VIPER_nonV600E_vs_V600E.WT <- ggplot(mrs_nonV600E_vs_WT.V600E_df, aes(x = TF, y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "NES") +
  theme_minimal(base_size = 12) +
  labs(title = "VIPER: BRAF nonV600E vs WT/V600E",
       x = "Regulator",
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")
VIPER_nonV600E_vs_V600E.WT
# ggsave("VIPER_nonV600E_vs_V600E.WT.pdf", VIPER_nonV600E_vs_V600E.WT)

#shadow analysis
mrshadow_nonV600E_vs_WT.V600E <- shadow(mrs_nonV600E_vs_WT.V600E, regulators = 50, verbose = FALSE)
summary(mrshadow_nonV600E_vs_WT.V600E)

shadowed_TF_nonV600E_vs_WT.V600E <- sapply(strsplit(summary(mrshadow_nonV600E_vs_WT.V600E)$Shadow.pairs, " -> "), `[`, 2)

all_TF_nonV600E_vs_WT.V600E <- summary(mrshadow_nonV600E_vs_WT.V600E)$msviper.results$Regulon

non_shadowed_TF_nonV600E_vs_WT.V600E <- setdiff(all_TF_nonV600E_vs_WT.V600E, shadowed_TF_nonV600E_vs_WT.V600E)

mrs_df <- summary(mrs_nonV600E_vs_WT.V600E, length(mrs_nonV600E_vs_WT.V600E$es$nes))
mrs_df$TF <- rownames(mrs_df)

# Filter only non-shadowed TFs
filtered_df <- mrs_df[mrs_df$TF %in% non_shadowed_TF_nonV600E_vs_WT.V600E, ]

# Order for plotting
filtered_df <- filtered_df[order(filtered_df$NES), ]
filtered_df$TF <- factor(filtered_df$TF, levels = filtered_df$TF)

library(ggplot2)

VIPER_nonV600E_vs_V600E.WT_shadow <- ggplot(filtered_df, aes(x = TF, y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "NES") +
  theme_minimal(base_size = 12) +
  labs(title = "VIPER (Non-Shadowed TFs): BRAF nonV600E vs WT/V600E",
       x = "Regulator (Non-Shadowed)",
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")
VIPER_nonV600E_vs_V600E.WT_shadow
# ggsave("VIPER_nonV600E_vs_V600E.WT_shadow.pdf", VIPER_nonV600E_vs_V600E.WT_shadow)

#WT VS NONV600E/V600E
dge_WT_vs_V600E.nonV600E <- DGEList(counts = MTAB_counts_wt_vs_v600e_nonv600e)

dge_WT_vs_V600E.nonV600E <- calcNormFactors(dge_WT_vs_V600E.nonV600E)

cutoff <- 18
keep_genes_WT_vs_V600E.nonV600E <- apply(cpm(dge_WT_vs_V600E.nonV600E), 1, max) >= cutoff
dge_WT_vs_V600E.nonV600E_filtered <- dge_WT_vs_V600E.nonV600E[keep_genes_WT_vs_V600E.nonV600E, ]


metadata_MTAB_wt_vs_v600e_nonv600e$mutation <- relevel(factor(metadata_MTAB_wt_vs_v600e_nonv600e$mutation), ref = "nonV600E_V600E")

design_WT_vs_V600E.nonV600E <- model.matrix(~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation, data = metadata_MTAB_wt_vs_v600e_nonv600e)

v_WT_vs_V600E.nonV600E <- voom(dge_WT_vs_V600E.nonV600E_filtered, design_WT_vs_V600E.nonV600E, plot = TRUE)
fit_WT_vs_V600E.nonV600E <- lmFit(v_WT_vs_V600E.nonV600E, design_WT_vs_V600E.nonV600E)
fit_WT_vs_V600E.nonV600E <- eBayes(fit_WT_vs_V600E.nonV600E)

topTable_WT_vs_V600E.nonV600E <- topTable(fit_WT_vs_V600E.nonV600E, coef = "mutationwt", adjust.method = "BH", number = Inf)
topTable_WT_vs_V600E.nonV600E$ENSEMBL <- rownames(topTable_WT_vs_V600E.nonV600E)


gencode_v36 <- rtracklayer::import("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.annotation.gtf")
entrez_annot_v36 <- data.table::fread("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.metadata.EntrezGene.gz",
                                      col.names = c("transcript_id", "entrez_id"))


temp_v36 <- gencode_v36 %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  dplyr::select(gene_id, transcript_id) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  merge(entrez_annot_v36, by = "transcript_id", all.x = TRUE)

gencode_v36_genes <- gencode_v36[gencode_v36$type == "gene"]
gencode_v36_genes$entrez_id <- temp_v36$entrez_id[match(gencode_v36_genes$gene_id, temp_v36$gene_id, nomatch = NA)]


annot_df <- as.data.frame(gencode_v36_genes)[, c("gene_id", "entrez_id", "gene_name")]
annot_df$gene_id <- sub("\\..*", "", annot_df$gene_id)
topTable_WT_vs_V600E.nonV600E$gene_name <- annot_df$gene_name[match(topTable_WT_vs_V600E.nonV600E$ENSEMBL, annot_df$gene_id)]

#filter valid genes
topTable_WT_vs_V600E.nonV600E <- topTable_WT_vs_V600E.nonV600E[!is.na(topTable_WT_vs_V600E.nonV600E$gene_name), ]
topTable_WT_vs_V600E.nonV600E <- topTable_WT_vs_V600E.nonV600E[!duplicated(topTable_WT_vs_V600E.nonV600E$gene_name), ]

topTable_WT_vs_V600E.nonV600E$zscore <- qnorm(topTable_WT_vs_V600E.nonV600E$P.Value / 2, lower.tail = FALSE) * sign(topTable_WT_vs_V600E.nonV600E$t)
signature_WT_vs_V600E.nonV600E <- setNames(topTable_WT_vs_V600E.nonV600E$zscore, topTable_WT_vs_V600E.nonV600E$gene_name)

# set.seed(123)
# n_perm <- 1
# genes <- rownames(v_WT_vs_V600E.nonV600E$E)
# null_matrix_WT_vs_V600E.nonV600E <- matrix(NA, nrow = length(genes), ncol = n_perm)
# rownames(null_matrix_WT_vs_V600E.nonV600E) <- genes
# 
# for (i in 1:n_perm) {
#   metadata_MTAB_wt_vs_v600e_nonv600e$mutation_perm <- sample(metadata_MTAB_wt_vs_v600e_nonv600e$mutation)
#   design_perm <- model.matrix(~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation_perm,
#                               data = metadata_MTAB_wt_vs_v600e_nonv600e)
# 
#   fit_perm <- lmFit(v_WT_vs_V600E.nonV600E, design_perm)
#   fit_perm <- eBayes(fit_perm)
# 
#   null_matrix_WT_vs_V600E.nonV600E[, i] <- fit_perm$t[, "mutation_permwt"]
# }
# 
# save(v_WT_vs_V600E.nonV600E, metadata_MTAB_wt_vs_v600e_nonv600e, file = "VIPER_WT_vs_V600E.nonV600E.RData")

load("Null_Matrix_WT_vs_V600E.nonV600E.RData")


gencode_v36 <- rtracklayer::import("/CCBdata/users/arnau/MOBER/extdata/gencode.v36.annotation.gtf")

gencode_v36_genes <- gencode_v36[gencode_v36$type == "gene"]
gencode_df <- as.data.frame(gencode_v36_genes)[, c("gene_id", "gene_name")]
gencode_df <- gencode_df[!duplicated(gencode_df$gene_id), ]


ensembl_ids <- rownames(null_matrix_WT_vs_V600E.nonV600E)
gencode_df$gene_id <- sub("\\..*", "", gencode_df$gene_id)
gene_names <- gencode_df$gene_name[match(ensembl_ids, gencode_df$gene_id)]


null_matrix_mapped <- null_matrix_WT_vs_V600E.nonV600E[!is.na(gene_names), ]
rownames(null_matrix_mapped) <- gene_names[!is.na(gene_names)]

# remove duplicates
null_matrix_mapped_WT_vs_V600E.nonV600E <- null_matrix_mapped[!duplicated(rownames(null_matrix_mapped)), ]

df_residual <- min(fit_WT_vs_V600E.nonV600E$df.residual)

#convert entire matrix of t-stats to two-sided p-values
p_null <- 2 * pt(-abs(null_matrix_mapped_WT_vs_V600E.nonV600E), df = df_residual)

#convert to z-scores
null_matrix_mapped_WT_vs_V600E.nonV600E_ZSCORE <- qnorm(p_null / 2, lower.tail = FALSE) * sign(null_matrix_mapped_WT_vs_V600E.nonV600E)


network_MTAB <- read.delim("/CCBdata/projects/BRAF_regulon/Aracne_COAD_output_MTAB_protein_coding/network.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
network_MTAB <- network_MTAB[, 1:3]
write.table(network_MTAB, "network_MTAB.adj", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
Expression_Data_MTAB <- read_tsv("MTAB_ARACNE_samples.tsv")
Expression_Data_MTAB <- as.data.frame(Expression_Data_MTAB)
rownames(Expression_Data_MTAB) <- Expression_Data_MTAB[[1]]
Expression_Data_MTAB[[1]] <- NULL
Expression_Data_MTAB <- as.matrix(Expression_Data_MTAB)

# regulon_object <- aracne2regulon("network_MTAB.adj", Expression_Data_MTAB, verbose = FALSE)

# Run msVIPER
mrs_WT_vs_V600E.nonV600E <- msviper(signature_WT_vs_V600E.nonV600E, regulon_object, nullmodel = null_matrix_mapped_WT_vs_V600E.nonV600E_ZSCORE, verbose = TRUE)


mrs_WT_vs_V600E.nonV600E_df <- summary(mrs_WT_vs_V600E.nonV600E, 50)
mrs_WT_vs_V600E.nonV600E_df$TF <- rownames(mrs_WT_vs_V600E.nonV600E_df)
mrs_WT_vs_V600E.nonV600E_df <- mrs_WT_vs_V600E.nonV600E_df[order(mrs_WT_vs_V600E.nonV600E_df$NES), ]
mrs_WT_vs_V600E.nonV600E_df$TF <- factor(mrs_WT_vs_V600E.nonV600E_df$TF, levels = mrs_WT_vs_V600E.nonV600E_df$TF)


VIPER_wt_vs_nonv600e_v600e <- ggplot(mrs_WT_vs_V600E.nonV600E_df, aes(x = TF, y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "NES") +
  theme_minimal(base_size = 12) +
  labs(title = "VIPER: BRAF WT vs V600E/nonV600E",
       x = "Regulator",
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")
VIPER_wt_vs_nonv600e_v600e
# ggsave("VIPER_wt_vs_nonv600e_v600e.pdf", VIPER_wt_vs_nonv600e_v600e)

#shadow analysis
mrshadow_WT_vs_V600E.nonV600E <- shadow(mrs_WT_vs_V600E.nonV600E, regulators = 25, verbose = FALSE)
summary(mrshadow_WT_vs_V600E.nonV600E)

shadowed_TF_WT_vs_V600E.nonV600E <- sapply(strsplit(summary(mrshadow_WT_vs_V600E.nonV600E)$Shadow.pairs, " -> "), `[`, 2)

all_TF_WT_vs_V600E.nonV600E <- summary(mrshadow_WT_vs_V600E.nonV600E)$msviper.results$Regulon

non_shadowed_TF_WT_vs_V600E.nonV600E <- setdiff(all_TF_WT_vs_V600E.nonV600E, shadowed_TF_WT_vs_V600E.nonV600E)

mrs_WT_vs_V600E.nonV600E_df <- summary(mrs_WT_vs_V600E.nonV600E, length(mrs_WT_vs_V600E.nonV600E$es$nes))
mrs_WT_vs_V600E.nonV600E_df$TF <- rownames(mrs_WT_vs_V600E.nonV600E_df)

#filter only non-shadowed TFs
filtered_WT_vs_V600E.nonV600E_df <- mrs_WT_vs_V600E.nonV600E_df[mrs_WT_vs_V600E.nonV600E_df$TF %in% non_shadowed_TF_WT_vs_V600E.nonV600E, ]


filtered_WT_vs_V600E.nonV600E_df <- filtered_WT_vs_V600E.nonV600E_df[order(filtered_WT_vs_V600E.nonV600E_df$NES), ]
filtered_WT_vs_V600E.nonV600E_df$TF <- factor(filtered_WT_vs_V600E.nonV600E_df$TF, levels = filtered_WT_vs_V600E.nonV600E_df$TF)

VIPER_WT_vs_nonV600E_V600E_shadow <- ggplot(filtered_WT_vs_V600E.nonV600E_df, aes(x = TF, y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "NES") +
  theme_minimal(base_size = 12) +
  labs(title = "VIPER (Non-Shadowed TFs): BRAF WT vs V600E/nonV600E",
       x = "Regulator (Non-Shadowed)",
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")
VIPER_WT_vs_nonV600E_V600E_shadow
# ggsave("VIPER_WT_vs_V600E.nonV600E_shadow.pdf", VIPER_WT_vs_nonV600E_V600E_shadow)
