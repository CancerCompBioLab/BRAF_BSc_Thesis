#!/bin/bash

Aracne=./ARACNe-AP/Aracne.jar

ls -l 
[ -f "$Aracne" ] || { echo "ARACNe-AP Path $Aracne doesn't exist"; exit 1; } 
echo "ARACNe-AP Path:$Aracne"
echo "Calculating Threshold"

java -Xmx20G -jar $Aracne -e ./ARACNe_log2transformed_TPM_MTAB_782samples.tsv -o Aracne_COAD_output_MTAB_protein_coding -t ./TF_COAD.txt --pvalue 1E-8 \
--seed 1 --calculateThreshold
