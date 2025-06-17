#!/usr/bin/env bash
#SBATCH --job-name=ARACNe
#SBATCH --output=ARACNe_COAD_TCGA_run_%j.log
#SBATCH --error=ARACNe_COAD_TCGA_run_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=65gb
#SBATCH --partition=long
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tarasyuziv@vhio.net

set -euxo pipefail

#Aracne=/mnt/bioinfnas/computational/allinas/Programs/ARACNe-AP/Aracne.jar
ARACNe=/mnt/CCBdata/projects/BRAF_regulon/ARACNe-AP
CCBdata=/mnt/CCBdata/projects/BRAF_regulon
[ -d "$CCBdata" ] || { echo "$CCBdata doesn't exist"; exit 1; } && { echo "$CCBdata exists"; ls $CCBdata; }

singularity run --bind ${PWD},"$CCBdata","$ARACNe" --pwd "$CCBdata" /mnt/CCBdata/projects/BRAF_regulon/openjdk_8u151-jdk.sif ./run_ARACNe_consolidate.sh  

