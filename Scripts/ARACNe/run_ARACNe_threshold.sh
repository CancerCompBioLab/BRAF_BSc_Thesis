#!/bin/bash

Aracne=./path/ARACNe-AP/file

ls -l 
[ -f "$Aracne" ] || { echo "ARACNe-AP Path $Aracne doesn't exist"; exit 1; } 
echo "ARACNe-AP Path:$Aracne"
echo "Calculating Threshold"

java -Xmx20G -jar $Aracne -e ./samples -o output_folder -t Scripts/ARACNe/TF_COAD.txt --pvalue 1E-8 \
--seed 1 --calculateThreshold
