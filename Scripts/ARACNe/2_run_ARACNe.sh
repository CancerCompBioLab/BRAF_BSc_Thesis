#!/bin/bash

Aracne=path/to/your/ARACNE-AP/file

#threads_to_use=62
ls -l

[ -f "$Aracne" ] || { echo "ARACNe-AP Path $Aracne doesn't exist"; exit 1; } 
echo "ARACNe-AP Path: $Aracne"
echo "Using # threads: $1"

for i in {1..500}
do
  echo "Bootstrap $i started"
  date
  java -Xmx60G -jar $Aracne -e ./samples -o output_folder -t .Scripts/ARACNe/TF_COAD.txt --pvalue 1E-8 \
  --seed $i --threads $1
  echo "Bootstrap $i completed"
  date

done
