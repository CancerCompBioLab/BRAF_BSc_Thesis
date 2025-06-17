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

#differential expression analysis using TCGA cohorts. GSEA are also included.

#please execute previous scripts in order to have all the files/objects created.

diff_exp_TCGA_metadata_wt_vs_V600E <- final_diff_exp_TCGA_metadata %>% filter(final_diff_exp_TCGA_metadata$condition == "wt" | final_diff_exp_TCGA_metadata$condition == "BRAFV600E")

diff_exp_TCGA_counts_wt_vs_V600E <- final_diff_exp_TCGA_counts[, colnames(final_diff_exp_TCGA_counts) %in% rownames(diff_exp_TCGA_metadata_wt_vs_V600E)]

diff_exp_TCGA_metadata_wt_vs_V600E <- diff_exp_TCGA_metadata_wt_vs_V600E[
  rownames(diff_exp_TCGA_metadata_wt_vs_V600E) %in% colnames(diff_exp_TCGA_counts_wt_vs_V600E), 
]

#subset count matrix to include only samples present in the metadata
diff_exp_TCGA_counts_wt_vs_V600E <- diff_exp_TCGA_counts_wt_vs_V600E[, 
                                                                     colnames(diff_exp_TCGA_counts_wt_vs_V600E) %in% rownames(diff_exp_TCGA_metadata_wt_vs_V600E)]

diff_exp_TCGA_metadata_wt_vs_V600E <- diff_exp_TCGA_metadata_wt_vs_V600E[
  match(colnames(diff_exp_TCGA_counts_wt_vs_V600E), rownames(diff_exp_TCGA_metadata_wt_vs_V600E)), 
]

#check
all(colnames(diff_exp_TCGA_counts_wt_vs_V600E) == rownames(diff_exp_TCGA_metadata_wt_vs_V600E))

diff_exp_TCGA_metadata_wt_vs_V600E <- diff_exp_TCGA_metadata_wt_vs_V600E %>%
  filter(!is.na(gender) & !is.na(MsiStatus) & !is.na(TumorPurity) & 
           !is.na(age_at_initial_pathologic_diagnosis) & !is.na(Site))

#subset
diff_exp_TCGA_counts_wt_vs_V600E <- diff_exp_TCGA_counts_wt_vs_V600E[, 
                                                                     colnames(diff_exp_TCGA_counts_wt_vs_V600E) %in% rownames(diff_exp_TCGA_metadata_wt_vs_V600E)]

#order
diff_exp_TCGA_metadata_wt_vs_V600E <- diff_exp_TCGA_metadata_wt_vs_V600E[
  match(colnames(diff_exp_TCGA_counts_wt_vs_V600E), rownames(diff_exp_TCGA_metadata_wt_vs_V600E)), ]

dds_wt_vs_v600e <- DESeqDataSetFromMatrix(
  countData = diff_exp_TCGA_counts_wt_vs_V600E,
  colData = diff_exp_TCGA_metadata_wt_vs_V600E,
  design = ~ gender + MsiStatus + TumorPurity + age_at_initial_pathologic_diagnosis + Site + condition
)
#remove stage

# keep <- rowSums(counts(dds_wt_vs_v600e) >= 10) >= 45
# dds_wt_vs_v600e <- dds_wt_vs_v600e[keep,]
# 
# dds_wt_vs_v600e <- DESeq(dds_wt_vs_v600e)
levels(dds_wt_vs_v600e$condition)

#ran this in the cluster.
#V600E VS WT
load("~/OnDemand/OnDemand_Results_TCGA/DESeq2_wt_vs_v600e_v2_results.RData")

dds_wt_vs_v600e
res_wt_vs_v600e <- results(dds_wt_vs_v600e, alpha = 0.05, contrast = c('condition', 'BRAFV600E','wt'))
res_wt_vs_v600e <- res_wt_vs_v600e[order(res_wt_vs_v600e$padj), ]
res_wt_vs_v600e_df <- as.data.frame(res_wt_vs_v600e)
res_wt_vs_v600e_df$Gene <- rownames(res_wt_vs_v600e_df)

res_wt_vs_v600e_df <- res_wt_vs_v600e_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_wt_vs_v600e_df, "DESeq2_wt_vs_v600e_v2_results.csv", row.names = FALSE)



DEG_wt_vs_V600E_v2 <- read.csv("DESeq2_wt_vs_v600e_v2_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

#ensmbl to gene name
DEG_wt_vs_V600E_v2 <- DEG_wt_vs_V600E_v2 %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))


DEG_wt_vs_V600E_v2 <- DEG_wt_vs_V600E_v2 %>% filter(!is.na(log2FoldChange))
DEG_wt_vs_V600E_v2 <- DEG_wt_vs_V600E_v2 %>% filter(!is.na(gene_symbol))
# DEG_wt_vs_V600E <- DEG_wt_vs_V600E %>% filter(DEG_wt_vs_V600E$pvalue < 0.05)#get rid of the
ranked_genes_wt_vs_v600e_v2 <- DEG_wt_vs_V600E_v2$log2FoldChange
names(ranked_genes_wt_vs_v600e_v2) <- DEG_wt_vs_V600E_v2$gene_symbol

ranked_genes_wt_vs_v600e_v2 <- sort(ranked_genes_wt_vs_v600e_v2, decreasing=TRUE)
ranked_genes_wt_vs_v600e_v2 <- ranked_genes_wt_vs_v600e_v2[!duplicated(names(ranked_genes_wt_vs_v600e_v2))]


library(EnhancedVolcano)

v600e_vs_wt <- EnhancedVolcano(DEG_wt_vs_V600E_v2,
                               lab = DEG_wt_vs_V600E_v2$gene_symbol,
                               x = "log2FoldChange",
                               y = "padj",
                               pCutoff = 0.05,
                               FCcutoff = 1,
                               title = "BRAF V600E vs BRAF WT - TCGA",
                               legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                               pointSize = 2.0,
                               labSize = 3.0)
v600e_vs_wt

# ggsave("REPORT_VOLCANO_BRAFV600E_VS_WT_TCGA.pdf", v600e_vs_wt, height = 10, width = 12)

#GSEA
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_wt_vs_v600e_v2 <- fgsea(pathways = pathway_list, 
                                      stats = ranked_genes_wt_vs_v600e_v2) #not specify permutations

fgsea_results_wt_vs_v600e_v2 <- fgsea_results_wt_vs_v600e_v2[order(fgsea_results_wt_vs_v600e_v2$padj), ]

# head(fgsea_results_wt_vs_v600e_v2 %>% arrange(padj), n=15)

#hallmarks
plot_res <- fgsea_results_wt_vs_v600e_v2 %>% filter(padj < 0.05)
hallmarks_results <- ggplot(plot_res, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E vs BRAF WT TCGA pval<0.05", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
hallmarks_results
ggsave("REPORT_HALLMARKS_BRAFV600E_vs_WT_TCGA_pval.pdf", hallmarks_results)

#KEGG
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)


fgsea_results_wt_vs_v600e_kegg <- fgsea(pathways = kegg_pathway_list, 
                                        stats = ranked_genes_wt_vs_v600e_v2)

fgsea_results_wt_vs_v600e_kegg <- fgsea_results_wt_vs_v600e_kegg[order(fgsea_results_wt_vs_v600e_kegg$padj), ]

plot_res_kegg <- fgsea_results_wt_vs_v600e_kegg %>% filter(padj < 0.05)
tcga_kegg_v600e_vs_wt <- ggplot(plot_res_kegg, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E vs BRAF WT TCGA", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
# ggsave("REPORT_KEGG_BRAFV600E_vs_WT.pdf", tcga_kegg_v600e_vs_wt)

#REACTOME
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_pathway_list <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

fgsea_results_wt_vs_v600e_reactome <- fgsea(pathways = reactome_pathway_list, 
                                            stats = ranked_genes_wt_vs_v600e_v2)

fgsea_results_wt_vs_v600e_reactome <- fgsea_results_wt_vs_v600e_reactome[order(fgsea_results_wt_vs_v600e_reactome$padj), ]

fgsea_results_wt_vs_v600e_reactome_PLOT <- fgsea_results_wt_vs_v600e_reactome %>% filter(padj < 0.05)
tcga_reactome_v600e_vs_wt<- ggplot(fgsea_results_wt_vs_v600e_reactome_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E VS BRAF WT TCGA", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
# ggsave("REPORT_REACTOME_BRAFV600E_WT.pdf", tcga_reactome_v600e_vs_wt, width = 10)




# nonV600E vs WT
load("~/OnDemand/OnDemand_Results_TCGA/DESeq2_wt_vs_nonv600e_v2_results.RData")
dds_wt_vs_nonv600e

res_wt_vs_nonv600e <- results(dds_wt_vs_nonv600e, alpha = 0.05, contrast = c('condition', 'BRAFnonV600E','wt'))
res_wt_vs_nonv600e <- res_wt_vs_nonv600e[order(res_wt_vs_nonv600e$padj), ]
res_wt_vs_nonv600e_df <- as.data.frame(res_wt_vs_nonv600e)
res_wt_vs_nonv600e_df$Gene <- rownames(res_wt_vs_nonv600e_df)

res_wt_vs_nonv600e_df <- res_wt_vs_nonv600e_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_wt_vs_nonv600e_df, "DESeq2_wt_vs_nonv600e_v2_results.csv", row.names = FALSE)

DEG_wt_vs_nonV600E_v2 <- read.csv("DESeq2_wt_vs_nonv600e_v2_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_wt_vs_nonV600E_v2 <- DEG_wt_vs_nonV600E_v2 %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))

DEG_wt_vs_nonV600E_v2 <- DEG_wt_vs_nonV600E_v2 %>% filter(!is.na(log2FoldChange))
DEG_wt_vs_nonV600E_v2 <- DEG_wt_vs_nonV600E_v2 %>% filter(!is.na(gene_symbol))
# DEG_wt_vs_V600E <- DEG_wt_vs_V600E %>% filter(DEG_wt_vs_V600E$pvalue < 0.05)#get rid of the
ranked_genes_wt_vs_nonv600e_v2 <- DEG_wt_vs_nonV600E_v2$log2FoldChange
names(ranked_genes_wt_vs_nonv600e_v2) <- DEG_wt_vs_nonV600E_v2$gene_symbol

ranked_genes_wt_vs_nonv600e_v2 <- sort(ranked_genes_wt_vs_nonv600e_v2, decreasing=TRUE)
ranked_genes_wt_vs_nonv600e_v2 <- ranked_genes_wt_vs_nonv600e_v2[!duplicated(names(ranked_genes_wt_vs_nonv600e_v2))]

nonv600e_vs_wt <- EnhancedVolcano(DEG_wt_vs_nonV600E_v2,
                                  lab = DEG_wt_vs_nonV600E_v2$gene_symbol,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  title = "BRAF nonV600E vs BRAF WT TCGA",
                                  legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                                  pointSize = 2.0,
                                  labSize = 3.0)
nonv600e_vs_wt
# ggsave("REPORT_VOLCANO_nonV600E_WT.pdf", nonv600e_vs_wt, width = 10, height = 10)
#hallamrks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_wt_vs_nonv600e_v2 <- fgsea(pathways = pathway_list, 
                                         stats = ranked_genes_wt_vs_nonv600e_v2) #not specify permutations

fgsea_results_wt_vs_nonv600e_v2 <- fgsea_results_wt_vs_nonv600e_v2[order(fgsea_results_wt_vs_nonv600e_v2$padj), ]

head(fgsea_results_wt_vs_nonv600e_v2)

fgsea_results_wt_vs_nonv600e_v2_plot <- fgsea_results_wt_vs_nonv600e_v2 %>% filter(pval< 0.05)
hallmars_nonv600e <- ggplot(fgsea_results_wt_vs_nonv600e_v2_plot, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAF nonV600E vs BRAFWT TCGA pval<0.05", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
hallmars_nonv600e
# ggsave("REPORT_HALLMARKS_nonv600e_vs_wt_tcga.pdf", hallmars_nonv600e)
#kegg
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_wt_vs_nonv600e_kegg <- fgsea(pathways = kegg_pathway_list, 
                                           stats = ranked_genes_wt_vs_nonv600e_v2)

fgsea_results_wt_vs_nonv600e_kegg <- fgsea_results_wt_vs_nonv600e_kegg[order(fgsea_results_wt_vs_nonv600e_kegg$padj), ]

fgsea_results_wt_vs_nonv600e_kegg_PLOT <- fgsea_results_wt_vs_nonv600e_kegg %>% filter(padj < 0.05)
tcga_kegg_nonv600e_vvs_wt <- ggplot(fgsea_results_wt_vs_nonv600e_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAF nonV600E vs BRAF WT TCGA", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
# ggsave("REPORT_KEGG_nonv600e_vs_wt.pdf", tcga_kegg_nonv600e_vvs_wt)
#reactome
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_pathway_list <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

fgsea_results_wt_vs_nonv600e_reactome <- fgsea(pathways = reactome_pathway_list, 
                                               stats = ranked_genes_wt_vs_nonv600e_v2)

fgsea_results_wt_vs_nonv600e_reactome <- fgsea_results_wt_vs_nonv600e_reactome[order(fgsea_results_wt_vs_nonv600e_reactome$padj), ]

fgsea_results_wt_vs_nonv600e_reactome_PLOT <- fgsea_results_wt_vs_nonv600e_reactome %>% filter(padj < 0.005)
reactome_nonv600e <- ggplot(fgsea_results_wt_vs_nonv600e_reactome_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRA FnonV600E vs BRAFWT TCGA", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

# ggsave("REPORT_REACTOME_nonV600E_wt.pdf", reactome_nonv600e, width = 10)
#V600E vs nonv600E
load("~/OnDemand/OnDemand_Results_TCGA/DESeq2_v600e_vs_nonv600e_v2_results.RData")
dds_wt_v600e_nonv600e

res_v600e_vs_nonv600e <- results(dds_wt_v600e_nonv600e, alpha = 0.05, contrast = c('condition', 'BRAFV600E','BRAFnonV600E'))
res_v600e_vs_nonv600e <- res_v600e_vs_nonv600e[order(res_v600e_vs_nonv600e$padj), ]
res_v600e_vs_nonv600e_df <- as.data.frame(res_v600e_vs_nonv600e)
res_v600e_vs_nonv600e_df$Gene <- rownames(res_v600e_vs_nonv600e_df)

res_v600e_vs_nonv600e_df <- res_v600e_vs_nonv600e_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_v600e_vs_nonv600e_df, "DESeq2_v600e_vs_nonv600e_v2_results.csv", row.names = FALSE)

DEG_V600E_vs_nonV600E_v2 <- read.csv("DESeq2_v600e_vs_nonv600e_v2_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_V600E_vs_nonV600E_v2 <- DEG_V600E_vs_nonV600E_v2 %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))
DEG_V600E_vs_nonV600E_v2 <- DEG_V600E_vs_nonV600E_v2 %>% filter(!is.na(log2FoldChange))
DEG_V600E_vs_nonV600E_v2 <- DEG_V600E_vs_nonV600E_v2 %>% filter(!is.na(gene_symbol))
# DEG_wt_vs_V600E <- DEG_wt_vs_V600E %>% filter(DEG_wt_vs_V600E$pvalue < 0.05)#get rid of the
ranked_genes_v600e_vs_nonv600e_v2 <- DEG_V600E_vs_nonV600E_v2$log2FoldChange
names(ranked_genes_v600e_vs_nonv600e_v2) <- DEG_V600E_vs_nonV600E_v2$gene_symbol

ranked_genes_v600e_vs_nonv600e_v2 <- sort(ranked_genes_v600e_vs_nonv600e_v2, decreasing=TRUE)
ranked_genes_v600e_vs_nonv600e_v2 <- ranked_genes_v600e_vs_nonv600e_v2[!duplicated(names(ranked_genes_v600e_vs_nonv600e_v2))]

v600e_nonv600e <- EnhancedVolcano(DEG_V600E_vs_nonV600E_v2,
                                  lab = DEG_V600E_vs_nonV600E_v2$gene_symbol,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  title = "BRAF V600E vs BRAF nonV600E TCGA",
                                  legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                                  pointSize = 2.0,
                                  labSize = 3.0)
v600e_nonv600e
# ggsave("REPORT_VOCLANO_V600E_VS_NONV600E_TCGA.pdf", v600e_nonv600e, width = 10, height = 10)
#hallmarks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)



fgsea_results_v600e_vs_nonv600e_v2 <- fgsea(pathways = pathway_list, 
                                            stats = ranked_genes_v600e_vs_nonv600e_v2) #not specify permutations

fgsea_results_v600e_vs_nonv600e_v2 <- fgsea_results_v600e_vs_nonv600e_v2[order(fgsea_results_v600e_vs_nonv600e_v2$padj), ]

head(fgsea_results_v600e_vs_nonv600e_v2)
#hallamrks
fgsea_results_v600e_vs_nonv600e_v2_plot <- fgsea_results_v600e_vs_nonv600e_v2 %>% filter(pval < 0.05)
hallmarks_combo <- ggplot(fgsea_results_v600e_vs_nonv600e_v2_plot, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E vs BRAFnonV600E", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_combo
#kegg
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_v600e_vs_nonv600e_v2_kegg <- fgsea(pathways = kegg_pathway_list, 
                                                 stats = ranked_genes_v600e_vs_nonv600e_v2)

fgsea_results_v600e_vs_nonv600e_v2_kegg <- fgsea_results_v600e_vs_nonv600e_v2_kegg[order(fgsea_results_v600e_vs_nonv600e_v2_kegg$padj), ]

fgsea_results_v600e_vs_nonv600e_v2_kegg_plot <- fgsea_results_v600e_vs_nonv600e_v2_kegg %>% filter(padj < 0.05)
tcga_kegg_v600e_nonv600e <- ggplot(fgsea_results_v600e_vs_nonv600e_v2_kegg_plot, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E VS BRAFnonV600E", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
#reactome
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_pathway_list <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

fgsea_results_v600e_vs_nonv600e_reactome <- fgsea(pathways = reactome_pathway_list, 
                                                  stats = ranked_genes_v600e_vs_nonv600e_v2)

fgsea_results_v600e_vs_nonv600e_reactome <- fgsea_results_v600e_vs_nonv600e_reactome[order(fgsea_results_v600e_vs_nonv600e_reactome$padj), ]

fgsea_results_v600e_vs_nonv600e_reactome_PLOT <- fgsea_results_v600e_vs_nonv600e_reactome %>% filter(padj < 0.05)
ggplot(fgsea_results_v600e_vs_nonv600e_reactome_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E vs BRAFnonV600E REACTOME", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()


