#!/bin/bash

Aracne=ARACNe-AP/Aracne.jar

#threads_to_use=62
ls -l

[ -f "$Aracne" ] || { echo "ARACNe-AP Path $Aracne doesn't exist"; exit 1; } 
echo "ARACNe-AP Path: $Aracne"
echo "Using # threads: $1"

for i in {1..500}
do
  echo "Bootstrap $i started"
  date
  java -Xmx60G -jar $Aracne -e ./ARACNe_log2transformed_TPM_MTAB_782samples.tsv -o Aracne_COAD_output_MTAB_protein_coding -t ./TF_COAD.txt --pvalue 1E-8 \
  --seed $i --threads $1
  echo "Bootstrap $i completed"
  date

done
