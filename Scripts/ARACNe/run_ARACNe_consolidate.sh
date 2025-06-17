#!/bin/bash

Aracne=./PATH/ARACNE-AP/file

ls -l 
[ -f "$Aracne" ] || { echo "ARACNe-AP Path $Aracne doesn't exist"; exit 1; } 
echo "ARACNe-AP Path:$Aracne"
echo "Calculating Threshold"

java -Xmx16G -jar $Aracne -o output_folder --consolidate
