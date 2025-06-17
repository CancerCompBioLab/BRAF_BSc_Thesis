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

#loading all the data TCGA. Please change paths to your local directories
gene_expression_data <- read_rds("/path_to_the_file")
mutation_data <- read_rds("/path_to_the_file")
metadata <-readr::read_tsv("Datasets/215982clinic_sup_grassso_cancerDisc18.tsv")
metadata_COAD <- read.csv("Datasets/COAD_clinic.csv")
metadata_PANCAN<-readr::read_tsv("Datasets/clinical_PANCAN_patient_with_followup.tsv")

#reading
mutation_data <- read.maf(mutation_data)
#raw counts
expression_unstranded_data <- assays(gene_expression_data)[["unstranded"]]
#mutation data for TCGA
mutation_data <- mutation_data@data
V600E_BRAF_mutations <- mutation_data %>% filter(mutation_data$Hugo_Symbol == "BRAF" & mutation_data$HGVSp_Short == "p.V640E") %>% dplyr::select(Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)

# double check because there might be classes.
non_V600E_BRAF_mutations <- mutation_data %>% filter(mutation_data$Hugo_Symbol == "BRAF" & mutation_data$HGVSp_Short != "p.V640E") %>% dplyr::select(Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)
wt_BRAF_mutations <- mutation_data %>% filter(mutation_data$Hugo_Symbol != "BRAF") %>% dplyr::select(Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)

#V600E BRAF
#I create a shorter barcode to match them with the sample names of the expression (they are unique if we take 16 characters)
V600E_BRAF_mutations$short_barcode <- substr(V600E_BRAF_mutations$Tumor_Sample_Barcode, 1, 16)
expression_colnames_short <- substr(colnames(expression_unstranded_data), 1, 16)

#And I selected the samples of the V600E BRAF
selected_columns_V600E <- expression_colnames_short %in% V600E_BRAF_mutations$short_barcode
filtered_V600E_expression <- expression_unstranded_data[, selected_columns_V600E]

#non-V600E BRAF
non_V600E_BRAF_mutations$short_barcode <- substr(non_V600E_BRAF_mutations$Tumor_Sample_Barcode, 1, 16)
#expression_colnames_short <- substr(colnames(expression_unstranded_data), 1, 16)

selected_columns_nonV600E <- expression_colnames_short %in% non_V600E_BRAF_mutations$short_barcode
filtered_nonV600E_expression <- expression_unstranded_data[, selected_columns_nonV600E]



#Wild type BRAF
wt_BRAF_mutations$short_barcode <- substr(wt_BRAF_mutations$Tumor_Sample_Barcode, 1, 16)
#expression_colnames_short <- substr(colnames(expression_unstranded_data), 1, 16)

selected_columns_wtBRAF <- expression_colnames_short %in% wt_BRAF_mutations$short_barcode
filtered_wtBRAF_expression <- expression_unstranded_data[, selected_columns_wtBRAF]


duplicated_samples_V600E_nonV600E <- intersect(colnames(filtered_nonV600E_expression), colnames(filtered_V600E_expression))
duplicated_samples_V600E_Wild <- intersect(colnames(filtered_V600E_expression), colnames(filtered_wtBRAF_expression))
duplicated_samples_nonV600E_Wild <- intersect(colnames(filtered_nonV600E_expression), colnames(filtered_wtBRAF_expression))

columns_to_remove <- duplicated_samples_V600E_nonV600E

filtered_V600E_expression <- filtered_V600E_expression[, 
                                                       !(colnames(filtered_V600E_expression) %in% columns_to_remove)]

filtered_nonV600E_expression <- filtered_nonV600E_expression[, 
                                                             !(colnames(filtered_nonV600E_expression) %in% columns_to_remove)]

columns_to_remove_V600E_Wildtype <- duplicated_samples_V600E_Wild

filtered_wtBRAF_expression <- filtered_wtBRAF_expression[, 
                                                         !(colnames(filtered_wtBRAF_expression) %in% columns_to_remove_V600E_Wildtype)]

columns_to_remove_nonV600E_Wildtype <- duplicated_samples_nonV600E_Wild


filtered_wtBRAF_expression <- filtered_wtBRAF_expression[, 
                                                         !(colnames(filtered_wtBRAF_expression) %in% columns_to_remove_nonV600E_Wildtype)]


#building metadata tcga
sample_info_wt_vs_v600e <- data.frame(
  sample = c(colnames(filtered_wtBRAF_expression),
             colnames(filtered_V600E_expression),
             colnames(filtered_nonV600E_expression)),
  condition = c(rep("wt", ncol(filtered_wtBRAF_expression)),
                rep("BRAFV600E", ncol(filtered_V600E_expression)),
                rep("BRAFnonV600E", ncol(filtered_nonV600E_expression))))

sample_info_wt_vs_v600e$short_sample <- substr(sample_info_wt_vs_v600e$sample, 1, 12)
metadata_sex <- metadata_COAD %>% filter(metadata_COAD$bcr_patient_barcode %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(bcr_patient_barcode, gender)
metadata_msistatus <- metadata %>% filter(metadata$Sample %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(Sample, MsiStatus)
metadata_tumorpurity <- metadata %>% filter(metadata$Sample %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(Sample, TumorPurity)
metadata_age <- metadata_COAD %>% filter(metadata_COAD$bcr_patient_barcode %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis)
metadata_stage <- metadata_PANCAN %>% filter(metadata_PANCAN$bcr_patient_barcode %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(bcr_patient_barcode, pathologic_stage)
metadata_location <- metadata_PANCAN %>% filter(metadata_PANCAN$bcr_patient_barcode %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(bcr_patient_barcode, anatomic_neoplasm_subdivision)
# metadata_location_2 <- metadata_PANCAN %>% filter(metadata_PANCAN$bcr_patient_barcode %in% sample_info_wt_vs_v600e$short_sample) %>% dplyr::select(bcr_patient_barcode, anatomic_neoplasm_subdivision)

#Left
metadata_location$Site <- ""
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision == "Sigmoid Colon"] = "left"
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision=="Descending Colon"] = "left"
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision=="Splenic Flexure"] = "left"

#Right
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision=="Cecum"] = "right"
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision=="Ascending Colon"] = "right"
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision=="Hepatic Flexure"] = "right"
metadata_location$Site[metadata_location$anatomic_neoplasm_subdivision=="Transverse Colon"] = "right"

#Stages
metadata_stage$STAGEn <- 0
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage I"] = 1
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IA"] = 1
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IB"] = 1
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage II"] = 2
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IIA"] = 2
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IIB"] = 2
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IIC"] = 2
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage III"] = 3
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IIIA"] = 3
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IIIB"] = 3
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IIIC"] = 3
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IV"] = 4
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IVA"] = 4
metadata_stage$STAGEn[metadata_stage$pathologic_stage=="Stage IVB"] = 4

#final metadata TCGA
final_metadata_wt_vs_v600e <- sample_info_wt_vs_v600e %>%
  left_join(metadata_sex, by = c("short_sample" = "bcr_patient_barcode")) %>%
  left_join(metadata_msistatus, by = c("short_sample" = "Sample")) %>%
  left_join(metadata_tumorpurity, by = c("short_sample" = "Sample")) %>%
  left_join(metadata_age, by = c("short_sample" = "bcr_patient_barcode")) %>%
  left_join(metadata_stage, by = c("short_sample" = "bcr_patient_barcode")) %>%
  left_join(metadata_location, by = c("short_sample" = "bcr_patient_barcode"))

final_metadata_wt_vs_v600e <- final_metadata_wt_vs_v600e %>% dplyr::select(-short_sample, -pathologic_stage)
