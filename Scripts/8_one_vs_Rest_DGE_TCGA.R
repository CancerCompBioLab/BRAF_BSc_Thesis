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


#1 vs "rest" comparison in the TCGA cohort

#V600E vs WT+nonV600E

diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e <- final_diff_exp_TCGA_metadata

diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e$condition[!(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e$condition == "BRAFV600E")] <- "nonV600E_wt"

diff_exp_TCGA_counts_V600E_VS_wt_nonv600e <- final_diff_exp_TCGA_counts[, colnames(final_diff_exp_TCGA_counts) %in% rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e)]



diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e <- diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e[
  rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e) %in% colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e), 
]

#subset count matrix 
diff_exp_TCGA_counts_V600E_VS_wt_nonv600e <- diff_exp_TCGA_counts_V600E_VS_wt_nonv600e[, 
                                                                                       colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e) %in% rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e)]

diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e <- diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e[
  match(colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e), rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e)), 
]

#check
all(colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e) == rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e))

diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e <- diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e %>%
  filter(!is.na(gender) & !is.na(MsiStatus) & !is.na(TumorPurity) & 
           !is.na(age_at_initial_pathologic_diagnosis) & !is.na(Site))

#subset 
diff_exp_TCGA_counts_V600E_VS_wt_nonv600e <- diff_exp_TCGA_counts_V600E_VS_wt_nonv600e[, 
                                                                                       colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e) %in% rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e)]

#order
diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e <- diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e[
  match(colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e), rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e)), ]


all(colnames(diff_exp_TCGA_counts_V600E_VS_wt_nonv600e) == rownames(diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e))



dds_TCGA_v600e_vs_nonv600e_wt <- DESeqDataSetFromMatrix(countData = diff_exp_TCGA_counts_V600E_VS_wt_nonv600e,
                                                        colData = diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e,
                                                        design = ~ gender + MsiStatus + TumorPurity + age_at_initial_pathologic_diagnosis + Site + condition)

# save(dds_TCGA_v600e_vs_nonv600e_wt, diff_exp_TCGA_counts_V600E_VS_wt_nonv600e, diff_exp_TCGA_metadata_V600E_VS_wt_nonv600e, file = "DESeq2_TCGA_v600e_vs_nonv600e_wt.Rdata")

# load("DESeq2_TCGA_v600e_vs_nonv600e_wt.Rdata")
# 
# keep <- rowSums(counts(dds_TCGA_v600e_vs_nonv600e_wt) >= 10) >= 45
# dds_TCGA_v600e_vs_nonv600e_wt <- dds_TCGA_v600e_vs_nonv600e_wt[keep,]
#  
# dds_TCGA_v600e_vs_nonv600e_wt <- DESeq(dds_TCGA_v600e_vs_nonv600e_wt)
# save(dds_TCGA_v600e_vs_nonv600e_wt, file = "DESeq2_TCGA_v600e_vs_nonv600e_wt_RESULTS.Rdata")

## ALL THE RESULTS (data frames) FROM THE DIFFERENTIAL EXPRESSION ANALYSIS ARE IN THE RESULTS FOLDER. FOR TCGA, E-MTAB-12862, AND ONE-VS-REST (TCGA). These codes were run in the cluster. You can use the loaded dataframes to check the GSEA.
load("DESeq2_TCGA_v600e_vs_nonv600e_wt_RESULTS.Rdata")

dds_TCGA_v600e_vs_nonv600e_wt
res_v600e_vs_nonv600e_wt <- results(dds_TCGA_v600e_vs_nonv600e_wt, alpha = 0.05, contrast = c('condition', 'BRAFV600E','nonV600E_wt'))
res_v600e_vs_nonv600e_wt <- res_v600e_vs_nonv600e_wt[order(res_v600e_vs_nonv600e_wt$padj), ]
res_v600e_vs_nonv600e_wt_df <- as.data.frame(res_v600e_vs_nonv600e_wt)
res_v600e_vs_nonv600e_wt_df$Gene <- rownames(res_v600e_vs_nonv600e_wt_df)

res_v600e_vs_nonv600e_wt_df <- res_v600e_vs_nonv600e_wt_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_v600e_vs_nonv600e_wt_df, "DESeq2_TCGA_v600e_vs_nonv600e_wt_results.csv", row.names = FALSE)


DEG_V600E_vs_nonV600E_wt <- read.csv("DESeq2_TCGA_v600e_vs_nonv600e_wt_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)


DEG_V600E_vs_nonV600E_wt <- DEG_V600E_vs_nonV600E_wt %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))


DEG_V600E_vs_nonV600E_wt <- DEG_V600E_vs_nonV600E_wt %>% filter(!is.na(log2FoldChange))
DEG_V600E_vs_nonV600E_wt <- DEG_V600E_vs_nonV600E_wt %>% filter(!is.na(gene_symbol))
# DEG_wt_vs_V600E <- DEG_wt_vs_V600E %>% filter(DEG_wt_vs_V600E$pvalue < 0.05)#get rid of the
ranked_genes_V600E_vs_nonV600E_wt <- DEG_V600E_vs_nonV600E_wt$log2FoldChange
names(ranked_genes_V600E_vs_nonV600E_wt) <- DEG_V600E_vs_nonV600E_wt$gene_symbol

ranked_genes_V600E_vs_nonV600E_wt <- sort(ranked_genes_V600E_vs_nonV600E_wt, decreasing=TRUE)
ranked_genes_V600E_vs_nonV600E_wt <- ranked_genes_V600E_vs_nonV600E_wt[!duplicated(names(ranked_genes_V600E_vs_nonV600E_wt))]

library(EnhancedVolcano)

#load TCGA_DGE_BRAFV600E_vs_BRAFnonV600E_WT.csv from github (Results/TCGA/One_vs_Rest). use: DEG_V600E_vs_nonV600E_wt <- read_csv("Results/TCGA/One_vs_Rest/TCGA_DGE_BRAFV600E_vs_BRAFnonV600E_WT.csv")
v600e_vs_nonv600e_wt <- EnhancedVolcano(DEG_V600E_vs_nonV600E_wt,
                                        lab = DEG_V600E_vs_nonV600E_wt$gene_symbol,
                                        x = "log2FoldChange",
                                        y = "padj",
                                        pCutoff = 0.05,
                                        FCcutoff = 1,
                                        title = "BRAF V600E vs BRAF WT+ nonV600E TCGA",
                                        legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                                        pointSize = 2.0,
                                        labSize = 3.0)
v600e_vs_nonv600e_wt
# ggsave("REPORT_VOLCANO_V600E_WT_NONV600E_tcga.pdf", v600e_vs_nonv600e_wt, height = 10, width = 12)


#hallmarks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_tcga_V600E_vs_nonV600E_wt <- fgsea(pathways = pathway_list, 
                                                 stats = ranked_genes_V600E_vs_nonV600E_wt) #not specify permutations

fgsea_results_tcga_V600E_vs_nonV600E_wt <- fgsea_results_tcga_V600E_vs_nonV600E_wt[order(fgsea_results_tcga_V600E_vs_nonV600E_wt$padj), ]

head(fgsea_results_tcga_V600E_vs_nonV600E_wt %>% arrange(padj), n=15)

plot_TCGA_V600E_vs_nonV600E_wt <- fgsea_results_tcga_V600E_vs_nonV600E_wt %>% filter(padj < 0.05)
hallmarks_TCGA_V600E_vs_nonV600E_wt <- ggplot(plot_TCGA_V600E_vs_nonV600E_wt, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAF V600E vs BRAFWT+nonV600E TCGA", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_TCGA_V600E_vs_nonV600E_wt
# ggsave("REPORT_HALLMARKS_V600VSWT_nonV600E.pdf", hallmarks_TCGA_V600E_vs_nonV600E_wt)

#kegg
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_MTAB_v600e_vs_nonv600_wt_e_kegg <- fgsea(pathways = kegg_pathway_list, 
                                                       stats = ranked_genes_V600E_vs_nonV600E_wt)

fgsea_results_MTAB_v600e_vs_nonv600_wt_e_kegg <- fgsea_results_MTAB_v600e_vs_nonv600_wt_e_kegg[order(fgsea_results_MTAB_v600e_vs_nonv600_wt_e_kegg$padj), ]

fgsea_results_MTAB_v600e_vs_nonv600e_wt_kegg_PLOT <- fgsea_results_MTAB_v600e_vs_nonv600_wt_e_kegg %>% filter(padj < 0.05)
tcga_combo_v600e_vs_non_wt <- ggplot(fgsea_results_MTAB_v600e_vs_nonv600e_wt_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="MTAB BRAFV600E vs nonV600E KEGG", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
tcga_combo_v600e_vs_non_wt



#nonV600E vs WT+V600E

diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e <- final_diff_exp_TCGA_metadata

diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e$condition[!(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e$condition == "BRAFnonV600E")] <- "V600E_wt"

diff_exp_TCGA_counts_nonV600E_VS_wt_v600e <- final_diff_exp_TCGA_counts[, colnames(final_diff_exp_TCGA_counts) %in% rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e)]

diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e <- diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e[
  rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e) %in% colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e), 
]

#subset count matrix
diff_exp_TCGA_counts_nonV600E_VS_wt_v600e <- diff_exp_TCGA_counts_nonV600E_VS_wt_v600e[, 
                                                                                       colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e) %in% rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e)]

diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e <- diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e[
  match(colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e), rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e)), 
]

#check
all(colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e) == rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e))

diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e <- diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e %>%
  filter(!is.na(gender) & !is.na(MsiStatus) & !is.na(TumorPurity) & 
           !is.na(age_at_initial_pathologic_diagnosis) & !is.na(Site))

#subset
diff_exp_TCGA_counts_nonV600E_VS_wt_v600e <- diff_exp_TCGA_counts_nonV600E_VS_wt_v600e[, 
                                                                                       colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e) %in% rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e)]

#order
diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e <- diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e[
  match(colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e), rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e)), ]


all(colnames(diff_exp_TCGA_counts_nonV600E_VS_wt_v600e) == rownames(diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e))


dds_TCGA_nonv600e_vs_v600e_wt <- DESeqDataSetFromMatrix(countData = diff_exp_TCGA_counts_nonV600E_VS_wt_v600e,
                                                        colData = diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e,
                                                        design = ~ gender + MsiStatus + TumorPurity + age_at_initial_pathologic_diagnosis + Site + condition)

#save(dds_TCGA_nonv600e_vs_v600e_wt, diff_exp_TCGA_counts_nonV600E_VS_wt_v600e, diff_exp_TCGA_metadata_nonV600E_VS_wt_v600e, file = "DESeq2_TCGA_nonv600e_vs_v600e_wt.Rdata")

load("DESeq2_TCGA_nonv600e_vs_v600e_wt_RESULTS.Rdata")

dds_TCGA_nonv600e_vs_v600e_wt
res_nonv600e_vs_v600e_wt <- results(dds_TCGA_nonv600e_vs_v600e_wt, alpha = 0.05, contrast = c('condition', 'BRAFnonV600E','V600E_wt'))
res_nonv600e_vs_v600e_wt <- res_nonv600e_vs_v600e_wt[order(res_nonv600e_vs_v600e_wt$padj), ]
res_nonv600e_vs_v600e_wt_df <- as.data.frame(res_nonv600e_vs_v600e_wt)
res_nonv600e_vs_v600e_wt_df$Gene <- rownames(res_nonv600e_vs_v600e_wt_df)

res_nonv600e_vs_v600e_wt_df <- res_nonv600e_vs_v600e_wt_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_nonv600e_vs_v600e_wt_df, "DESeq2_TCGA_nonv600e_vs_v600e_wt_results.csv", row.names = FALSE)

DEG_nonV600E_vs_V600E_wt <- read.csv("DESeq2_TCGA_nonv600e_vs_v600e_wt_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_nonV600E_vs_V600E_wt <- DEG_nonV600E_vs_V600E_wt %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))

DEG_nonV600E_vs_V600E_wt <- DEG_nonV600E_vs_V600E_wt %>% filter(!is.na(log2FoldChange))
DEG_nonV600E_vs_V600E_wt <- DEG_nonV600E_vs_V600E_wt %>% filter(!is.na(gene_symbol))
# DEG_wt_vs_V600E <- DEG_wt_vs_V600E %>% filter(DEG_wt_vs_V600E$pvalue < 0.05)#get rid of the
ranked_genes_nonV600E_vs_V600E_wt <- DEG_nonV600E_vs_V600E_wt$log2FoldChange
names(ranked_genes_nonV600E_vs_V600E_wt) <- DEG_nonV600E_vs_V600E_wt$gene_symbol

ranked_genes_nonV600E_vs_V600E_wt <- sort(ranked_genes_nonV600E_vs_V600E_wt, decreasing=TRUE)
ranked_genes_nonV600E_vs_V600E_wt <- ranked_genes_nonV600E_vs_V600E_wt[!duplicated(names(ranked_genes_nonV600E_vs_V600E_wt))]

library(EnhancedVolcano)
#load TCGA_DGE_BRAFnonV600E_vs_BRAFV600E_WT.csv from github (Results/TCGA/One_vs_Rest). use: DEG_nonV600E_vs_V600E_wt <- read_csv("Results/TCGA/One_vs_Rest/TCGA_DGE_BRAFnonV600E_vs_BRAFV600E_WT.csv")
nonv600e_vs_v600e_wt <- EnhancedVolcano(DEG_nonV600E_vs_V600E_wt,
                                        lab = DEG_nonV600E_vs_V600E_wt$gene_symbol,
                                        x = "log2FoldChange",
                                        y = "padj",
                                        pCutoff = 0.05,
                                        FCcutoff = 1,
                                        title = "BRAF nonV600E vs BRAFWT + V600E TCGA",
                                        legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                                        pointSize = 2.0,
                                        labSize = 3.0)
nonv600e_vs_v600e_wt
# ggsave("REPORT_VOLCANO_nonv600e_vs_wt_v600e.pdf", nonv600e_vs_v600e_wt, width = 10, height = 10)

#hallmarks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_tcga_nonV600E_vs_V600E_wt <- fgsea(pathways = pathway_list, 
                                                 stats = ranked_genes_nonV600E_vs_V600E_wt) #not specify permutations


fgsea_results_tcga_nonV600E_vs_V600E_wt <- fgsea_results_tcga_nonV600E_vs_V600E_wt[order(fgsea_results_tcga_nonV600E_vs_V600E_wt$padj), ]

head(fgsea_results_tcga_nonV600E_vs_V600E_wt %>% arrange(padj), n=15)

plot_TCGA_nonV600E_vs_V600E_wt <- fgsea_results_tcga_nonV600E_vs_V600E_wt %>% filter(padj < 0.05)
hallmarks_TCGA_nonV600E_vs_V600E_wt <- ggplot(plot_TCGA_nonV600E_vs_V600E_wt, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAF nonV600E vs BRAFWT+V600E TCGA", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_TCGA_nonV600E_vs_V600E_wt
# ggsave("REPORT_HALLAMRKS_nonv600e_vs_wt_v600e.pdf", hallmarks_TCGA_nonV600E_vs_V600E_wt)
#kegg
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg <- fgsea(pathways = kegg_pathway_list, 
                                                       stats = ranked_genes_nonV600E_vs_V600E_wt)

fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg <- fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg[order(fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg$padj), ]

fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg_PLOT <- fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg %>% filter(pval < 0.05)
tcga_combo_nonv600e_vs_v600e_wt <- ggplot(fgsea_results_MTAB_nonv600e_vs_v600_wt_e_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="MTAB BRAFnonV600E vs V600E_wt KEGG", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
tcga_combo_nonv600e_vs_v600e_wt


#WT VS NONV600E + V600E

diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e <- final_diff_exp_TCGA_metadata

diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e$condition[!(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e$condition == "wt")] <- "V600E_nonV600E"

diff_exp_TCGA_counts_wt_VS_v600e_nonv600e <- final_diff_exp_TCGA_counts[, colnames(final_diff_exp_TCGA_counts) %in% rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e)]

diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e <- diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e[
  rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e) %in% colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e), 
]

#subset
diff_exp_TCGA_counts_wt_VS_v600e_nonv600e <- diff_exp_TCGA_counts_wt_VS_v600e_nonv600e[, 
                                                                                       colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e) %in% rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e)]

diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e <- diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e[
  match(colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e), rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e)), 
]


all(colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e) == rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e))


diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e <- diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e %>%
  filter(!is.na(gender) & !is.na(MsiStatus) & !is.na(TumorPurity) & 
           !is.na(age_at_initial_pathologic_diagnosis) & !is.na(Site))

# Subset the counts matrix to match the filtered metadata
diff_exp_TCGA_counts_wt_VS_v600e_nonv600e <- diff_exp_TCGA_counts_wt_VS_v600e_nonv600e[, 
                                                                                       colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e) %in% rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e)]

# Ensure correct order
diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e <- diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e[
  match(colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e), rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e)), ]


all(colnames(diff_exp_TCGA_counts_wt_VS_v600e_nonv600e) == rownames(diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e))

dds_TCGA_wt_vs_v600e_nonv600e <- DESeqDataSetFromMatrix(countData = diff_exp_TCGA_counts_wt_VS_v600e_nonv600e,
                                                        colData = diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e,
                                                        design = ~ gender + MsiStatus + TumorPurity + age_at_initial_pathologic_diagnosis + Site + condition)

# save(dds_TCGA_wt_vs_v600e_nonv600e, diff_exp_TCGA_counts_wt_VS_v600e_nonv600e, diff_exp_TCGA_metadata_wt_vs_v600e_nonv600e, file = "DESeq2_TCGA_wt_vs_v600e_nonv600e.Rdata")

load("DESeq2_TCGA_wt_vs_v600e_nonv600e_RESULTS.Rdata")

dds_TCGA_wt_vs_v600e_nonv600e
res_wt_vs_v600e_nonv600e <- results(dds_TCGA_wt_vs_v600e_nonv600e, alpha = 0.05, contrast = c('condition', 'wt','V600E_nonV600E'))
res_wt_vs_v600e_nonv600e <- res_wt_vs_v600e_nonv600e[order(res_wt_vs_v600e_nonv600e$padj), ]
res_wt_vs_v600e_nonv600e_df <- as.data.frame(res_wt_vs_v600e_nonv600e)
res_wt_vs_v600e_nonv600e_df$Gene <- rownames(res_wt_vs_v600e_nonv600e_df)

res_wt_vs_v600e_nonv600e_df <- res_wt_vs_v600e_nonv600e_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_wt_vs_v600e_nonv600e_df, "DESeq2_TCGA_wt_vs_v600e_nonv600e_results.csv", row.names = FALSE)

DEG_wt_vs_V600E_nonV600E <- read.csv("DESeq2_TCGA_wt_vs_v600e_nonv600e_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_wt_vs_V600E_nonV600E <- DEG_wt_vs_V600E_nonV600E %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))

DEG_wt_vs_V600E_nonV600E <- DEG_wt_vs_V600E_nonV600E %>% filter(!is.na(log2FoldChange))
DEG_wt_vs_V600E_nonV600E <- DEG_wt_vs_V600E_nonV600E %>% filter(!is.na(gene_symbol))
# DEG_wt_vs_V600E <- DEG_wt_vs_V600E %>% filter(DEG_wt_vs_V600E$pvalue < 0.05)#get rid of the
ranked_genes_wt_vs_V600E_nonV600E <- DEG_wt_vs_V600E_nonV600E$log2FoldChange
names(ranked_genes_wt_vs_V600E_nonV600E) <- DEG_wt_vs_V600E_nonV600E$gene_symbol

ranked_genes_wt_vs_V600E_nonV600E <- sort(ranked_genes_wt_vs_V600E_nonV600E, decreasing=TRUE)
ranked_genes_wt_vs_V600E_nonV600E <- ranked_genes_wt_vs_V600E_nonV600E[!duplicated(names(ranked_genes_wt_vs_V600E_nonV600E))]

library(EnhancedVolcano)

#load TCGA_DGE_BRAFWT_vs_BRAFV600E_nonV600E.csv from github (Results/TCGA/One_vs_Rest). use: DEG_wt_vs_V600E_nonV600E <- read_csv("Results/TCGA/One_vs_Rest/TCGA_DGE_BRAFWT_vs_BRAFV600E_nonV600E.csv")
nonv600e_vs_v600e_wt <- EnhancedVolcano(DEG_wt_vs_V600E_nonV600E,
                                        lab = DEG_wt_vs_V600E_nonV600E$gene_symbol,
                                        x = "log2FoldChange",
                                        y = "padj",
                                        pCutoff = 0.05,
                                        FCcutoff = 1,
                                        title = "BRAF wt vs V600E_nonV600E",
                                        legendLabels = c("No sig", "Log2 FC", "p-value", "Both"),
                                        pointSize = 2.0,
                                        labSize = 3.0)
nonv600e_vs_v600e_wt

#hallmarks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_tcga_wt_vs_nonV600E_v600e <- fgsea(pathways = pathway_list, 
                                                 stats = ranked_genes_wt_vs_V600E_nonV600E) #not specify permutations

fgsea_results_tcga_wt_vs_nonV600E_v600e <- fgsea_results_tcga_wt_vs_nonV600E_v600e[order(fgsea_results_tcga_wt_vs_nonV600E_v600e$padj), ]

head(fgsea_results_tcga_wt_vs_nonV600E_v600e %>% arrange(padj), n=15)

plot_TCGA_wt_vs_V600E_nonv600e <- fgsea_results_tcga_wt_vs_nonV600E_v600e %>% filter(pval < 0.05)
hallmarks_TCGA_wt_vs_nonV600E_v600e <- ggplot(plot_TCGA_wt_vs_V600E_nonv600e, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="TCGA BRAF wt vs nonV600E+V600E", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_TCGA_wt_vs_nonV600E_v600e

msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg <- fgsea(pathways = kegg_pathway_list, 
                                                       stats = ranked_genes_wt_vs_V600E_nonV600E)

fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg <- fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg[order(fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg$padj), ]

fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg_PLOT <- fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg %>% filter(pval < 0.05)
tcga_combo_wt_vs_v600e_nonv600e <- ggplot(fgsea_results_MTAB_wt_vs_v600_nonv600e_e_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="MTAB wy vs V600E_nonv600e KEGG", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
tcga_combo_wt_vs_v600e_nonv600e


