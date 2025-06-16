#loading the mtab mutation data information from the vcf files. And performing DEA and GSEA

library(stringr)
sample_ids_MTAB <- metadata_MTAB$`DNA Tumor Sample Barcode`
grep_lines <- readLines("/CCBdata/users/lea/BRAF_Postdoc_Project/Datasets/RNA-seq/E-MTAB-12862/Mutation/results_CRC_Braf/annotation/project/braf_grep_results.txt")


parse_braf_line <- function(line) {
  sample <- str_match(line, "^([^/]+)/")[,2]
  protein_change <- str_match(line, "p\\.\\w+\\d+\\w+|p\\.\\w+\\*|p\\.\\w+fs")[,1]
  data.frame(sample = sample, mutation = protein_change, stringsAsFactors = FALSE)
}

parsed <- do.call(rbind, lapply(grep_lines, parse_braf_line))

parsed_unique <- parsed[!duplicated(parsed$sample), ]

final_df <- data.frame(sample = sample_ids_MTAB, stringsAsFactors = FALSE)
final_df <- merge(final_df, parsed_unique, by.x = "sample", by.y = "sample", all.x = TRUE)
final_df$mutation[is.na(final_df$mutation)] <- "wt"

colnames(metadata_MTAB)[3] <- "DNA_Tumor_Sample_Barcode"
colnames(final_df)[1] <- "DNA_Tumor_Sample_Barcode"
metadata_MTAB_mutation <- metadata_MTAB %>%
  inner_join(final_df, by = "DNA_Tumor_Sample_Barcode")

table(metadata_MTAB_mutation$mutation)
#39 nonV600E
#208 V600E
#816 wt
#including rectal

#differential gene expression analysis MTAB

metadata_MTAB_mutation_filtered <- metadata_MTAB_mutation %>% filter(!metadata_MTAB_mutation$`Anatomic Organ Subdivision` == "Rectum")
metadata_MTAB_mutation_filtered <- metadata_MTAB_mutation_filtered[, -1]

MTAB_V600E <- metadata_MTAB_mutation_filtered %>% filter(metadata_MTAB_mutation_filtered$mutation == "p.Val640Glu")

MTAB_nonV600E <- metadata_MTAB_mutation_filtered %>% filter(!(metadata_MTAB_mutation_filtered$mutation == "p.Val640Glu") & !(metadata_MTAB_mutation_filtered$mutation == "wt"))

MTAB_BRAF_WT <- metadata_MTAB_mutation_filtered %>% filter(metadata_MTAB_mutation_filtered$mutation == "wt")
#V600E 196 samples
#nonV600E 28 samples
#wt 558 samples

#V600E vs WT

MTAB_V600E_metadata_vs_wt <- metadata_MTAB_mutation_filtered %>% filter((metadata_MTAB_mutation_filtered$mutation == "p.Val640Glu") | (metadata_MTAB_mutation_filtered$mutation == "wt"))

MTAB_counts_diffexp <- MTAB_counts
colnames(MTAB_counts_diffexp) <- gsub("\\.", "-", colnames(MTAB_counts_diffexp))
rownames(MTAB_counts_diffexp) <- MTAB_counts_diffexp$`ENSG-ID`
MTAB_counts_diffexp <- MTAB_counts_diffexp[, -1]

MTAB_counts_filtered_v600e_vs_wt <- MTAB_counts_diffexp[, colnames(MTAB_counts_diffexp) %in% MTAB_V600E_metadata_vs_wt$DNA_Tumor_Sample_Barcode]

colnames(MTAB_V600E_metadata_vs_wt)[40] <- "MSI_status"
colnames(MTAB_V600E_metadata_vs_wt)[33] <- "Tumour_Cell_Content_Pathology"
colnames(MTAB_V600E_metadata_vs_wt)[17] <- "Age_at_diagnosis"
colnames(MTAB_V600E_metadata_vs_wt)[21] <- "Primary_Site_Disease"
colnames(MTAB_V600E_metadata_vs_wt)[23] <- "Tumour_Site"

rownames(MTAB_V600E_metadata_vs_wt) <- MTAB_V600E_metadata_vs_wt$DNA_Tumor_Sample_Barcode

remove_col_mtab <- setdiff(rownames(MTAB_V600E_metadata_vs_wt), colnames(MTAB_counts_filtered_v600e_vs_wt))
MTAB_V600E_metadata_vs_wt <- MTAB_V600E_metadata_vs_wt %>% filter(!(MTAB_V600E_metadata_vs_wt$DNA_Tumor_Sample_Barcode %in% remove_col_mtab))


MTAB_V600E_metadata_vs_wt$mutation <- as.factor(MTAB_V600E_metadata_vs_wt$mutation)
MTAB_V600E_metadata_vs_wt$MSI_status <- as.factor(MTAB_V600E_metadata_vs_wt$MSI_status)
MTAB_V600E_metadata_vs_wt$Sex <- as.factor(MTAB_V600E_metadata_vs_wt$Sex)
MTAB_V600E_metadata_vs_wt$Tumour_Cell_Content_Pathology <- as.factor(MTAB_V600E_metadata_vs_wt$Tumour_Cell_Content_Pathology)
MTAB_V600E_metadata_vs_wt$Age_at_diagnosis <- as.factor(MTAB_V600E_metadata_vs_wt$Age_at_diagnosis)
MTAB_V600E_metadata_vs_wt$Primary_Site_Disease <- as.factor(MTAB_V600E_metadata_vs_wt$Primary_Site_Disease)

rownames(MTAB_V600E_metadata_vs_wt) <- MTAB_V600E_metadata_vs_wt$DNA_Tumor_Sample_Barcode

dds_MTAB_diffexp <- DESeqDataSetFromMatrix(countData = MTAB_counts_filtered_v600e_vs_wt,
                                           colData = MTAB_V600E_metadata_vs_wt,
                                           design = ~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation)

# save(dds_MTAB_diffexp, MTAB_counts_filtered_v600e_vs_wt, MTAB_V600E_metadata_vs_wt, file = "DESeq_MTAB_V600E_vs_WT.Rdata")

# load("DESeq_MTAB_V600E_vs_WT.Rdata")
# 
# keep <- rowSums(counts(dds_MTAB_diffexp) >= 10) >= 195
# dds_MTAB_diffexp <- dds_MTAB_diffexp[keep,]
#  
# dds_MTAB_diffexp <- DESeq(dds_MTAB_diffexp)
# save(dds_MTAB_diffexp, file = "DESeq_MTAB_V600E_vs_WT.Rdata")

load("DESeq_MTAB_V600E_vs_WT_results.Rdata")

dds_MTAB_diffexp
res_dds_MTAB_diffexp <- results(dds_MTAB_diffexp, alpha = 0.05, contrast = c('mutation', 'p.Val640Glu','wt'))
res_dds_MTAB_diffexp <- res_dds_MTAB_diffexp[order(res_dds_MTAB_diffexp$padj), ]
res_dds_MTAB_diffexp_df <- as.data.frame(res_dds_MTAB_diffexp)
res_dds_MTAB_diffexp_df$Gene <- rownames(res_dds_MTAB_diffexp_df)

res_dds_MTAB_diffexp_df <- res_dds_MTAB_diffexp_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_dds_MTAB_diffexp_df, "DESeq2_MTAB_V600E_vs_WT_results.csv", row.names = FALSE)

DEG_MTAB_V600E_VS_WT <- read.csv("DESeq2_MTAB_V600E_vs_WT_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_MTAB_V600E_VS_WT <- DEG_MTAB_V600E_VS_WT %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))

DEG_MTAB_V600E_VS_WT <- DEG_MTAB_V600E_VS_WT %>% filter(!is.na(log2FoldChange))
DEG_MTAB_V600E_VS_WT <- DEG_MTAB_V600E_VS_WT %>% filter(!is.na(gene_symbol))


ranked_genes_mtab_v600e_vs_wt <- DEG_MTAB_V600E_VS_WT$log2FoldChange
names(ranked_genes_mtab_v600e_vs_wt) <- DEG_MTAB_V600E_VS_WT$gene_symbol

ranked_genes_mtab_v600e_vs_wt <- sort(ranked_genes_mtab_v600e_vs_wt, decreasing=TRUE)
ranked_genes_mtab_v600e_vs_wt <- ranked_genes_mtab_v600e_vs_wt[!duplicated(names(ranked_genes_mtab_v600e_vs_wt))]

library(EnhancedVolcano)

MTAB_v600e_vs_wt_volcano <- EnhancedVolcano(DEG_MTAB_V600E_VS_WT,
                                            lab = DEG_MTAB_V600E_VS_WT$gene_symbol,
                                            x = "log2FoldChange",
                                            y = "padj",
                                            pCutoff = 0.05,
                                            FCcutoff = 1,
                                            title = "BRAF V600E vs BRAF WT MTAB",
                                            legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                                            pointSize = 2.0,
                                            labSize = 3.0)
MTAB_v600e_vs_wt_volcano
# ggsave("REPORT_VOLCANO_BRAFV600E_vs_WT_MTAB.pdf", MTAB_v600e_vs_wt_volcano, height = 12, width = 12)

#HALLMARKS
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_mtab_v600e_vs_wt <- fgsea(pathways = pathway_list, 
                                        stats = ranked_genes_mtab_v600e_vs_wt) #not specify permutations

fgsea_results_mtab_v600e_vs_wt <- fgsea_results_mtab_v600e_vs_wt[order(fgsea_results_mtab_v600e_vs_wt$padj), ]

head(fgsea_results_mtab_v600e_vs_wt %>% arrange(padj), n=15)

plot_mtab_V600e_wt <- fgsea_results_mtab_v600e_vs_wt %>% filter(padj < 0.05)
hallmarks_results_mtab_v600e_wt <- ggplot(plot_mtab_V600e_wt, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E vs BRAF WT MTAB", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_results_mtab_v600e_wt
# ggsave("REPORT_HALLMARKS_V600E_VS_WT_MTAB.pdf", hallmarks_results_mtab_v600e_wt)
#KEGG
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_MTAB_v600e_vs_wt_kegg <- fgsea(pathways = kegg_pathway_list, 
                                             stats = ranked_genes_mtab_v600e_vs_wt)

fgsea_results_MTAB_v600e_vs_wt_kegg <- fgsea_results_MTAB_v600e_vs_wt_kegg[order(fgsea_results_MTAB_v600e_vs_wt_kegg$padj), ]

fgsea_results_MTAB_v600e_vs_wt_kegg_PLOT <- fgsea_results_MTAB_v600e_vs_wt_kegg %>% filter(padj < 0.05)
kegg_mtab_V600e_wt <- ggplot(fgsea_results_MTAB_v600e_vs_wt_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E vs BRAF WT MTAB", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

# ggsave("REPORT_KEGG_BRAFV600E_WT_MTAB.pdf", kegg_mtab_V600e_wt)
#reactome
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_pathway_list <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

fgsea_results_mtab_v600e_vs_wt_reactome <- fgsea(pathways = reactome_pathway_list, 
                                                 stats = ranked_genes_mtab_v600e_vs_wt)

fgsea_results_mtab_v600e_vs_wt_reactome <- fgsea_results_mtab_v600e_vs_wt_reactome[order(fgsea_results_mtab_v600e_vs_wt_reactome$padj), ]

fgsea_results_mtab_v600e_vs_wt_reactome_PLOT <- fgsea_results_mtab_v600e_vs_wt_reactome %>% filter(padj < 0.00005)
reactome_mtab <- ggplot(fgsea_results_mtab_v600e_vs_wt_reactome_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAFV600E VS BRAF WT MTAB", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
# ggsave("REPORT_REACTOME_BAFV600E_WT_MTAB.pdf", reactome_mtab, width = 12)


#nonV600E vs WT
MTAB_metadata_nonv600e_vs_wt <- metadata_MTAB_mutation_filtered %>% filter((metadata_MTAB_mutation_filtered$mutation == "wt") | !(metadata_MTAB_mutation_filtered$mutation == "p.Val640Glu"))
MTAB_metadata_nonv600e_vs_wt$mutation[MTAB_metadata_nonv600e_vs_wt$mutation != "wt"] <- "nonV600E"

MTAB_counts_diffexp <- MTAB_counts
colnames(MTAB_counts_diffexp) <- gsub("\\.", "-", colnames(MTAB_counts_diffexp))
rownames(MTAB_counts_diffexp) <- MTAB_counts_diffexp$`ENSG-ID`
MTAB_counts_diffexp <- MTAB_counts_diffexp[, -1]

MTAB_counts_filtered_nonv600e_vs_wt <- MTAB_counts_diffexp[, colnames(MTAB_counts_diffexp) %in% MTAB_metadata_nonv600e_vs_wt$DNA_Tumor_Sample_Barcode]


colnames(MTAB_metadata_nonv600e_vs_wt)[40] <- "MSI_status"
colnames(MTAB_metadata_nonv600e_vs_wt)[33] <- "Tumour_Cell_Content_Pathology"
colnames(MTAB_metadata_nonv600e_vs_wt)[17] <- "Age_at_diagnosis"
colnames(MTAB_metadata_nonv600e_vs_wt)[21] <- "Primary_Site_Disease"
colnames(MTAB_metadata_nonv600e_vs_wt)[23] <- "Tumour_Site"

rownames(MTAB_metadata_nonv600e_vs_wt) <- MTAB_metadata_nonv600e_vs_wt$DNA_Tumor_Sample_Barcode
remove_col_mtab_nonv600e_vs_wt <- setdiff(rownames(MTAB_metadata_nonv600e_vs_wt), colnames(MTAB_counts_filtered_nonv600e_vs_wt))

MTAB_metadata_nonv600e_vs_wt <- MTAB_metadata_nonv600e_vs_wt %>% filter(!(MTAB_metadata_nonv600e_vs_wt$DNA_Tumor_Sample_Barcode %in% remove_col_mtab_nonv600e_vs_wt))

dds_MTAB_diffexp_nonv600e_vs_wt <- DESeqDataSetFromMatrix(countData = MTAB_counts_filtered_nonv600e_vs_wt,
                                                          colData = MTAB_metadata_nonv600e_vs_wt,
                                                          design = ~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation)

# save(dds_MTAB_diffexp_nonv600e_vs_wt, MTAB_counts_filtered_nonv600e_vs_wt, MTAB_metadata_nonv600e_vs_wt, file = "DESeq2_MTAB_nonV600E_vs_wt.RData")

# load("DESeq2_MTAB_nonV600E_vs_wt.RData")
# 
# keep <- rowSums(counts(dds_MTAB_diffexp_nonv600e_vs_wt) >= 10) >= 28
# dds_MTAB_diffexp_nonv600e_vs_wt <- dds_MTAB_diffexp_nonv600e_vs_wt[keep,]
#  
# dds_MTAB_diffexp_nonv600e_vs_wt <- DESeq(dds_MTAB_diffexp_nonv600e_vs_wt)
# save(dds_MTAB_diffexp_nonv600e_vs_wt, file = "DESeq2_MTAB_nonV600E_vs_wt_results.Rdata")


load("DESeq2_MTAB_nonV600E_vs_wt_results.Rdata")

dds_MTAB_diffexp_nonv600e_vs_wt
res_dds_MTAB_diffexp_nonv600e_vs_wt <- results(dds_MTAB_diffexp_nonv600e_vs_wt, alpha = 0.05, contrast = c('mutation', 'nonV600E','wt'))
res_dds_MTAB_diffexp_nonv600e_vs_wt <- res_dds_MTAB_diffexp_nonv600e_vs_wt[order(res_dds_MTAB_diffexp_nonv600e_vs_wt$padj), ]
res_dds_MTAB_diffexp_nonv600e_vs_wt_df <- as.data.frame(res_dds_MTAB_diffexp_nonv600e_vs_wt)
res_dds_MTAB_diffexp_nonv600e_vs_wt_df$Gene <- rownames(res_dds_MTAB_diffexp_nonv600e_vs_wt_df)

res_dds_MTAB_diffexp_nonv600e_vs_wt_df <- res_dds_MTAB_diffexp_nonv600e_vs_wt_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_dds_MTAB_diffexp_nonv600e_vs_wt_df, "DESeq2_MTAB_nonV600E_vs_WT_results.csv", row.names = FALSE)

DEG_MTAB_nonV600E_VS_WT <- read.csv("DESeq2_MTAB_nonV600E_vs_WT_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_MTAB_nonV600E_VS_WT <- DEG_MTAB_nonV600E_VS_WT %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))

DEG_MTAB_nonV600E_VS_WT <- DEG_MTAB_nonV600E_VS_WT %>% filter(!is.na(log2FoldChange))
DEG_MTAB_nonV600E_VS_WT <- DEG_MTAB_nonV600E_VS_WT %>% filter(!is.na(gene_symbol))


ranked_genes_mtab_nonv600e_vs_wt <- DEG_MTAB_nonV600E_VS_WT$log2FoldChange
names(ranked_genes_mtab_nonv600e_vs_wt) <- DEG_MTAB_nonV600E_VS_WT$gene_symbol

ranked_genes_mtab_nonv600e_vs_wt <- sort(ranked_genes_mtab_nonv600e_vs_wt, decreasing=TRUE)
ranked_genes_mtab_nonv600e_vs_wt <- ranked_genes_mtab_nonv600e_vs_wt[!duplicated(names(ranked_genes_mtab_nonv600e_vs_wt))]

library(EnhancedVolcano)

MTAB_nonv600e_vs_wt_volcano <- EnhancedVolcano(DEG_MTAB_nonV600E_VS_WT,
                                               lab = DEG_MTAB_nonV600E_VS_WT$gene_symbol,
                                               x = "log2FoldChange",
                                               y = "padj",
                                               pCutoff = 0.05,
                                               FCcutoff = 1,
                                               title = "BRAF nonV600E vs BRAFWT MTAB",
                                               legendLabels = c("No sig", "Log2 FC", "p-adj", "Both"),
                                               pointSize = 2.0,
                                               labSize = 3.0)
MTAB_nonv600e_vs_wt_volcano
# ggsave("REPORT_VOLCANO_nonv600e_wt_mtab.pdf", MTAB_nonv600e_vs_wt_volcano, width = 10, height = 10)

#hallmarks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_mtab_nonv600e_vs_wt <- fgsea(pathways = pathway_list, 
                                           stats = ranked_genes_mtab_nonv600e_vs_wt) #not specify permutations

fgsea_results_mtab_nonv600e_vs_wt <- fgsea_results_mtab_nonv600e_vs_wt[order(fgsea_results_mtab_nonv600e_vs_wt$padj), ]

head(fgsea_results_mtab_nonv600e_vs_wt %>% arrange(padj), n=15)

plot_mtab_nonV600e_wt <- fgsea_results_mtab_nonv600e_vs_wt %>% filter(padj < 0.05)
hallmarks_results_mtab_nonv600e_wt <- ggplot(plot_mtab_nonV600e_wt, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="BRAF nonV600E vs BRAF WT MTAB", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_results_mtab_nonv600e_wt
# ggsave("REPORT_HALLMARKS_nonv600e_wt_mtab.pdf", hallmarks_results_mtab_nonv600e_wt)

#KEGG
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_MTAB_nonv600e_vs_wt_kegg <- fgsea(pathways = kegg_pathway_list, 
                                                stats = ranked_genes_mtab_nonv600e_vs_wt)

fgsea_results_MTAB_nonv600e_vs_wt_kegg <- fgsea_results_MTAB_nonv600e_vs_wt_kegg[order(fgsea_results_MTAB_nonv600e_vs_wt_kegg$padj), ]

fgsea_results_MTAB_nonv600e_vs_wt_kegg_PLOT <- fgsea_results_MTAB_nonv600e_vs_wt_kegg %>% filter(pval < 0.05)
ggplot(fgsea_results_MTAB_nonv600e_vs_wt_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="MTAB BRAFnonV600E vs WT KEGG", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

#rectome
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_pathway_list <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

fgsea_results_mtab_nonv600e_vs_wt_reactome <- fgsea(pathways = reactome_pathway_list, 
                                                    stats = ranked_genes_mtab_nonv600e_vs_wt)

fgsea_results_mtab_nonv600e_vs_wt_reactome <- fgsea_results_mtab_nonv600e_vs_wt_reactome[order(fgsea_results_mtab_nonv600e_vs_wt_reactome$padj), ]

fgsea_results_mtab_nonv600e_vs_wt_reactome_PLOT <- fgsea_results_mtab_nonv600e_vs_wt_reactome %>% filter(padj < 0.05)
ggplot(fgsea_results_mtab_nonv600e_vs_wt_reactome_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title=" MTAB BRAF nonV600E VS WT REACTOME", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

#V600E vs nonV600E

MTAB_metadata_v600e_vs_nonv600e <- metadata_MTAB_mutation_filtered %>% filter(!(metadata_MTAB_mutation_filtered$mutation == "wt") | (metadata_MTAB_mutation_filtered$mutation == "p.Val640Glu"))

MTAB_metadata_v600e_vs_nonv600e$mutation[MTAB_metadata_v600e_vs_nonv600e$mutation != "p.Val640Glu"] <- "nonV600E"

MTAB_counts_diffexp <- MTAB_counts
colnames(MTAB_counts_diffexp) <- gsub("\\.", "-", colnames(MTAB_counts_diffexp))
rownames(MTAB_counts_diffexp) <- MTAB_counts_diffexp$`ENSG-ID`
MTAB_counts_diffexp <- MTAB_counts_diffexp[, -1]

MTAB_counts_filtered_v600e_vs_nonv600e <- MTAB_counts_diffexp[, colnames(MTAB_counts_diffexp) %in% MTAB_metadata_v600e_vs_nonv600e$DNA_Tumor_Sample_Barcode]

colnames(MTAB_metadata_v600e_vs_nonv600e)[40] <- "MSI_status"
colnames(MTAB_metadata_v600e_vs_nonv600e)[33] <- "Tumour_Cell_Content_Pathology"
colnames(MTAB_metadata_v600e_vs_nonv600e)[17] <- "Age_at_diagnosis"
colnames(MTAB_metadata_v600e_vs_nonv600e)[21] <- "Primary_Site_Disease"
colnames(MTAB_metadata_v600e_vs_nonv600e)[23] <- "Tumour_Site"

rownames(MTAB_metadata_v600e_vs_nonv600e) <- MTAB_metadata_v600e_vs_nonv600e$DNA_Tumor_Sample_Barcode
remove_col_mtab_v600e_vs_nonv600e <- setdiff(rownames(MTAB_metadata_v600e_vs_nonv600e), colnames(MTAB_counts_filtered_v600e_vs_nonv600e))

MTAB_metadata_v600e_vs_nonv600e <- MTAB_metadata_v600e_vs_nonv600e %>% filter(!(MTAB_metadata_v600e_vs_nonv600e$DNA_Tumor_Sample_Barcode %in% remove_col_mtab_v600e_vs_nonv600e))

MTAB_metadata_v600e_vs_nonv600e$Sex <- factor(MTAB_metadata_v600e_vs_nonv600e$Sex)
MTAB_metadata_v600e_vs_nonv600e$MSI_status <- factor(MTAB_metadata_v600e_vs_nonv600e$MSI_status)
MTAB_metadata_v600e_vs_nonv600e$Tumour_Cell_Content_Pathology <- as.numeric(MTAB_metadata_v600e_vs_nonv600e$Tumour_Cell_Content_Pathology)
MTAB_metadata_v600e_vs_nonv600e$Tumour_Site <- factor(MTAB_metadata_v600e_vs_nonv600e$Tumour_Site)
MTAB_metadata_v600e_vs_nonv600e$mutation <- factor(MTAB_metadata_v600e_vs_nonv600e$mutation)


dds_MTAB_diffexp_v600e_vs_nonv600e <- DESeqDataSetFromMatrix(countData = MTAB_counts_filtered_v600e_vs_nonv600e,
                                                             colData = MTAB_metadata_v600e_vs_nonv600e,
                                                             design = ~ Sex + MSI_status + Tumour_Cell_Content_Pathology + Age_at_diagnosis + Tumour_Site + mutation)

save(dds_MTAB_diffexp_v600e_vs_nonv600e, MTAB_counts_filtered_v600e_vs_nonv600e, MTAB_metadata_v600e_vs_nonv600e, file = "DESeq2_MTAB_V600E_vs_nonV600E.RData")

# load("DESeq2_MTAB_V600E_vs_nonV600E.RData")
# 
# keep <- rowSums(counts(dds_MTAB_diffexp_v600e_vs_nonv600e) >= 10) >= 28
# dds_MTAB_diffexp_v600e_vs_nonv600e <- dds_MTAB_diffexp_v600e_vs_nonv600e[keep,]
#  
# dds_MTAB_diffexp_v600e_vs_nonv600e <- DESeq(dds_MTAB_diffexp_v600e_vs_nonv600e)
# save(dds_MTAB_diffexp_v600e_vs_nonv600e, file = "DESeq2_MTAB_V600E_vs_nonV600E_results.RData")


load("DESeq2_MTAB_V600E_vs_nonV600E_results.RData")

dds_MTAB_diffexp_v600e_vs_nonv600e
res_dds_MTAB_diffexp_v600e_vs_nonv600e <- results(dds_MTAB_diffexp_v600e_vs_nonv600e, alpha = 0.05, contrast = c('mutation', 'p.Val640Glu','nonV600E'))
res_dds_MTAB_diffexp_v600e_vs_nonv600e <- res_dds_MTAB_diffexp_v600e_vs_nonv600e[order(res_dds_MTAB_diffexp_v600e_vs_nonv600e$padj), ]
res_dds_MTAB_diffexp_v600e_vs_nonv600e_df <- as.data.frame(res_dds_MTAB_diffexp_v600e_vs_nonv600e)
res_dds_MTAB_diffexp_v600e_vs_nonv600e_df$Gene <- rownames(res_dds_MTAB_diffexp_v600e_vs_nonv600e_df)

res_dds_MTAB_diffexp_v600e_vs_nonv600e_df <- res_dds_MTAB_diffexp_v600e_vs_nonv600e_df %>%
  mutate(Gene = sub("\\..*", "", Gene))

write.csv(res_dds_MTAB_diffexp_v600e_vs_nonv600e_df, "DESeq2_MTAB_V600E_vs_nonV600E_results.csv", row.names = FALSE)

DEG_MTAB_V600E_VS_nonV600E <- read.csv("DESeq2_MTAB_V600E_vs_nonV600E_results.csv")
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

DEG_MTAB_V600E_VS_nonV600E <- DEG_MTAB_V600E_VS_nonV600E %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, 
                              keys = Gene, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first"))

DEG_MTAB_V600E_VS_nonV600E <- DEG_MTAB_V600E_VS_nonV600E %>% filter(!is.na(log2FoldChange))
DEG_MTAB_V600E_VS_nonV600E <- DEG_MTAB_V600E_VS_nonV600E %>% filter(!is.na(gene_symbol))


ranked_genes_mtab_v600e_vs_nonv600e <- DEG_MTAB_V600E_VS_nonV600E$log2FoldChange
names(ranked_genes_mtab_v600e_vs_nonv600e) <- DEG_MTAB_V600E_VS_nonV600E$gene_symbol

ranked_genes_mtab_v600e_vs_nonv600e <- sort(ranked_genes_mtab_v600e_vs_nonv600e, decreasing=TRUE)
ranked_genes_mtab_v600e_vs_nonv600e <- ranked_genes_mtab_v600e_vs_nonv600e[!duplicated(names(ranked_genes_mtab_v600e_vs_nonv600e))]

library(EnhancedVolcano)

MTAB_v600e_vs_nonv600e_volcano <- EnhancedVolcano(DEG_MTAB_V600E_VS_nonV600E,
                                                  lab = DEG_MTAB_V600E_VS_nonV600E$gene_symbol,
                                                  x = "log2FoldChange",
                                                  y = "padj",
                                                  pCutoff = 0.05,
                                                  FCcutoff = 1,
                                                  title = "MTAB BRAF V600E vs nonV600E",
                                                  legendLabels = c("No sig", "Log2 FC", "p-value", "Both"),
                                                  pointSize = 2.0,
                                                  labSize = 3.0)
MTAB_v600e_vs_nonv600e_volcano

#hallmarks
msigdb_pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathway_list <- split(msigdb_pathways$gene_symbol, msigdb_pathways$gs_name)

fgsea_results_mtab_v600e_vs_nonv600e <- fgsea(pathways = pathway_list, 
                                              stats = ranked_genes_mtab_v600e_vs_nonv600e) #not specify permutations

fgsea_results_mtab_v600e_vs_nonv600e <- fgsea_results_mtab_v600e_vs_nonv600e[order(fgsea_results_mtab_v600e_vs_nonv600e$padj), ]

head(fgsea_results_mtab_v600e_vs_nonv600e %>% arrange(padj), n=15)

plot_mtab_v600e_nonv600e <- fgsea_results_mtab_v600e_vs_nonv600e %>% filter(padj < 0.05)
hallmarks_results_mtab_v600e_nonv600e <- ggplot(plot_mtab_v600e_nonv600e, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="MTAB BRAF V600E vs nonV600E", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()

hallmarks_results_mtab_v600e_nonv600e
#kegg
msigdb_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_pathway_list <- split(msigdb_kegg$gene_symbol, msigdb_kegg$gs_name)

fgsea_results_MTAB_v600e_vs_nonv600e_kegg <- fgsea(pathways = kegg_pathway_list, 
                                                   stats = ranked_genes_mtab_v600e_vs_nonv600e)

fgsea_results_MTAB_v600e_vs_nonv600e_kegg <- fgsea_results_MTAB_v600e_vs_nonv600e_kegg[order(fgsea_results_MTAB_v600e_vs_nonv600e_kegg$padj), ]

fgsea_results_MTAB_v600e_vs_nonv600e_kegg_PLOT <- fgsea_results_MTAB_v600e_vs_nonv600e_kegg %>% filter(padj < 0.05)
ggplot(fgsea_results_MTAB_v600e_vs_nonv600e_kegg_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title="MTAB BRAFV600E vs nonV600E KEGG", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()


#reactome
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_pathway_list <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

fgsea_results_mtab_v600e_vs_nonv600e_reactome <- fgsea(pathways = reactome_pathway_list, 
                                                       stats = ranked_genes_mtab_v600e_vs_nonv600e)

fgsea_results_mtab_v600e_vs_nonv600e_reactome <- fgsea_results_mtab_v600e_vs_nonv600e_reactome[order(fgsea_results_mtab_v600e_vs_nonv600e_reactome$padj), ]

fgsea_results_mtab_v600e_vs_nonv600e_reactome_PLOT <- fgsea_results_mtab_v600e_vs_nonv600e_reactome %>% filter(padj < 0.05)
ggplot(fgsea_results_mtab_v600e_vs_nonv600e_reactome_PLOT, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill=NES > 0)) + 
  coord_flip() + 
  labs(title=" MTAB BRAF V600E VS nonV600E REACTOME", x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()
















