# Transcriptional and Regulatory Divergence Between BRAF-V600E and non-V600E Colon Cancers
In this project I aimed to look into the transcriptional and regulatory landscapes of BRAF-V600E and BRAF non-V600E colon cancers. By integrating gene expression data from multiple cohorts, we analyzed differential expression profiles and pathway enrichments between mutated (BRAF-V600E and BRAF nonV600E) and wild-type (BRAF WT) tumor samples. Primary analysis were conducted using TCGA-COAD dataset as discovery cohort, with replication in the independent E-MTAB-12862 cohort. In addition, we applied the VIPER algorithm to infer the activity of transcription factors, signaling proteins and chromatin regulatory genes.

# Folders

## Scripts
All the scripts are intended to be run sequentially. However the outputs of computationally intensive steps (e.g., DGE and VIPER analyses) are provided separately (Results folder) in case you wish to explore specific sections of the code without rerunning the full pipeline. (Please read the comments provided in the code to load specific data frames from the Results folder). R version 4.1.3.


### 0_TCGAbiolinks.R

Download expression and mutation data from TCGAbiolinks. In case you can not download the data using the code (you should be able), do it manually (https://portal.gdc.cancer.gov/projects/TCGA-COAD)

### 1_Data_Filtering_TCGA.R

Data loading, sample stratification and building metadata of TCGA (discovery cohort).

### 2_PreMSI.R

Using PreMSI from the WangX-lab (https://github.com/WangX-Lab/PreMSIm) to predicit missing metadata in TCGA.

### 3_Exploratory_Analysis_TCGA.R

Performing exploratory analysis in the TCGA cohort.

### 4_E_MTAB_12862_loading_and_Exploratory_analyisis.R

Loading E-MTAB-12862 dataset (validation cohort) and performing exploratory analysis.

### 5_PAMr_classifier.R

Training PAMr classifier and validating in annotated (left & right) samples of TCGA.

### 6_DEA_GSEA_TCGA.R

Performing Differential Expression Analysis and Gene Set Enrichment Anaylsis between BRAF samples in TCGA.

### 7_DEA_GSEA_MTAB.R

Performing Differential Expression Analysis and Gene Set Enrichment Anaylsis between BRAF samples in E-MTAB-12862.

### 8_one_vs_Rest_DGE_TCGA.R

Performing Differential Expression Analysis and Gene Set Enrichment Anaylsis using "one-vs-rest" method in TCGA.

### 9_VIPERxLIMMA_analysis.R

Loading regulons from ARACNe's output and performing VIPER analysis (limma was used to get the differentially expressed genes signatures).

### ARACNe

In this folder I included all the scripts I used to generate the regulons for TCGA and E-MTAB-12862 following the ARACNe pipeline (https://github.com/califano-lab/ARACNe-AP). Please follow the ARACNe pipeline in order to download all the software/scripts (e.g., ARACNE-AP) needed to build the regulons. You can also see the TF_COAD.txt file which contains the transcription factors, signaling proteins and chromatin regulatory genes used to build the regulons.

## Datasets

In this folder you can find the metadata for TCGA and E-MTAB-12862 fully annotated cohorts (Annotated_Metadatas folder). In addition, other metadata files, necessary to build the TCGA metadata, are included. I also uploaded some intermediate files generated from the scripts (please read the comments of my code). Raw counts and additional datasets are not in the folder due to their large size (you have the instructions in the scripts to download them).

## Results

This folder contains the results of the multiple Differential Expression Analysis performed (TCGA, E-MTAB-12862, 1vsRest TCGA), as well as the VIPER analysis results (non-Shadowed results included). Regulons from E-MTAB-12862 and TCGA are not included due to its large size.

## Manuscript

Here I uploaded my Supplementary Information and Main Manuscript of the project (PDF).

# Contact

For any questions, feedback, or other opportunities, please don't be afraid to reach out to me! You can email me at taras.yuziv@alum.esci.upf.edu or tarasfriend2004@gmail.com



