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
#TPM counts
exp_COAD<-assays(TCGA.COAD.exp)$tpm_unstrand



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
#Create a MAF file
maf <- read.maf(maf = mut_READ)
