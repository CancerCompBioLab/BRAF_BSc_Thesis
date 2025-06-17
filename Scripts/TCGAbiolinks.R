library(TCGAbiolinks)
library(dplyr)
library(maftools)
library(SummarizedExperiment)
#Query TCGA-READ gene expression data
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

#We need to shorten the sample name to 12 characters to merge it after with the clinical data file,
#which only contains 12 characters.
colnames(exp_COAD)<-substr(colnames(exp_COAD), 1, 16)
#Now we see that the sample names are shorter.
colnames(exp_COAD)[1:5]
#We want to keep only the tumour samples and remove any other samples. The tumour samples end with 01A.
#We therefore filter for these samples
exp_COAD<- dplyr::select(exp_COAD,contains("01A"))
#Now we have 163 cases, which is more similar to 167 cases stated above (difference possible as we select only 01A)
#Now we can further shorten the sample name so that it is identical with the sample name in clinical data
colnames(exp_COAD)<-substr(colnames(exp_COAD), 1, 12)
#All cases are unique
duplicated(colnames(exp_COAD))
#Now we see that the sample names are shorter.

colnames(exp_COAD)[1:5]


# Query available files for TCGA-READ
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