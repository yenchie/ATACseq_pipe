#!/bin/bash
# Program: ATACseq_peakAnnt
#SBATCH -J ATACseq_peakAnnt
#SBATCH -o ATACseq_peakAnnt.Gm724_%A.out
#SBATCH -n 4
#SBATCH --mem=8GB

set -euo pipefail

echo "$(scontrol show job $SLURM_JOBID)" | grep "JobId="
echo "$(scontrol show job $SLURM_JOBID)" | grep "Command="
echo ""
echo ""
echo ""
# Name: ATACseq.peakAnno.v.1.sh
# Description: This script is for peak files annotation
# Author: YCW
# Date: 2023/12/28
# Update: 
# Updata log:
# 1. 
# 2.

#
# Dependence:
module add R/4.1.2
#
# Configure Setting (if needed):
# # - all option
# #configure loading
# source $configure
# THREADS=$core
#
################################################################################
### USAGE: 
################################################################################
#		sbatch ATACseq.peakAnno.v.1.sh $1 $2 $3
#
# Input:
# - $1: the path of peak files.
# - $2: the path of output folder 
# - $3: the format of peak files
#
# Output:
# - .ACR.PeakAnnotation.txt named by sampleID will be created in the output folder ($3)
# - ChIPseeker.plot.pdf named by sampleID will be created in the output folder ($3)
# - ChIPseeker.plot.pdf named by a folder which peak files in will be created in the output folder ($3)
################################################################################


################################################################################
echo "peak files in ${1} ATACseqQC start"
echo "peak files format: ${3}"
echo "output path: ${2}"
echo ""
echo ""
echo ""

[ -d ${2} ] || mkdir -p ${2}

time Rscript --vanilla /bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/peakAnnt.Gm724.v.1.0.R \
            ${1} ${2} ${3} 
            
#USAGE:Rscript --vanilla /bcst/JYL/YCWang/Rscript/peakAnnt.Gm724.v.1.0.R [arg1] [arg2] [arg3]\n
          #Arg1: input folder path, ex. /bcst/JYL/YCWang/testing/ATAC_test/peakannt/bed \n
          #Arg2: output file path ex. /bcst/JYL/YCWang/testing/ATAC_test/peakannt \n
          #Arg3: type of input files ex. bed \n
echo ""   
echo ""
echo "peak files in ${1} ATACseqQC done"
