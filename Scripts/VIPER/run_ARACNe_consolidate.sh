#!/bin/bash

Aracne=./ARACNe-AP/Aracne.jar

ls -l 
[ -f "$Aracne" ] || { echo "ARACNe-AP Path $Aracne doesn't exist"; exit 1; } 
echo "ARACNe-AP Path:$Aracne"
echo "Calculating Threshold"

java -Xmx16G -jar $Aracne -o Aracne_COAD_output_TCGA_protein_coding --consolidate
