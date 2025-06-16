# Transcriptional and Regulatory Divergence Between BRAF-V600E and non-V600E Colon Cancers

# Folders

## Scripts

### 1_Data_Filtering_TCGA.R

Data loading, sample stratification and building metadata of TCGA (discovery cohort).

### 2_PreMSI.R

Using PreMSI from the WangX-lab (https://github.com/WangX-Lab/PreMSIm) to predicit missing metadata in TCGA.

### 3_Exploratory_Analysis_TCGA.R

Performing exploratory analysis in the TCGA cohort.

### 4_E_MTAB_12862_loading_and_Exploratory_analyisis.R

Loading E-MTAB-12862 dataset (validation cohort) and perfroming exploratory analysis.

### 5_PAMr_classifier.R

Training PAMr classifier and externally validating in annotated (left & right) samples of TCGA.

### 6_DEA_GSEA_TCGA.R

Performing Differential Expression Analysis and Gene Set Enrichment Anaylsis between BRAF samples in TCGA.

### 7_DEA_GSEA_MTAB.R

Performing Differential Expression Analysis and Gene Set Enrichment Anaylsis between BRAF samples in E-MTAB-12862.

### 8_one_vs_Rest_DGE_TCGA.R

Performing Differential Expression Analysis and Gene Set Enrichment Anaylsis using "one-vs-rest" method in TCGA.

### 9_VIPERxLIMMA_analysis.R

Loading regulons from ARACNe's output and performing VIPER analysis (limma was used to get the differentially expressed genes signatures)



