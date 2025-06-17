library(TCGAbiolinks)
library(dplyr)
library(maftools)
library(SummarizedExperiment)
#Query TCGA-COAD gene expression data
query_exp <- GDCquery(project = "TCGA-COAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query_exp)
#Create R object
TCGA.COAD.exp <- GDCprepare(query_exp)

#Raw counts are in unstranded
exp_COAD_raw<-assays(TCGA.COAD.exp)$unstranded
#Let's select the TPM counts
exp_COAD<-assays(TCGA.COAD.exp)$tpm_unstrand
#Convert to dataframe
exp_COAD<-as.data.frame(exp_COAD)


# Query available files for TCGA-COAD
query_mut <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Simple Nucleotide Variation",
  access='open'
)
GDCdownload(query_mut)
output_query_mutation<-getResults(query_mut)
#Create a SummarizedExperiment
mut_COAD<-GDCprepare(query_mut,summarizedExperiment = TRUE)
#Shorten the sample name so that we can combine it will clinical_data file
mut_COAD$Tumor_Sample_Barcode<-substr(mut_READ$Tumor_Sample_Barcode, 1, 12)
#Create a MAF file
# A MAF file is a tab-delimited text format designed specifically to store somatic mutation data for cancer genomics.
#It is a standardized format used by The Cancer Genome Atlas (TCGA).
maf <- read.maf(maf = mut_READ)
