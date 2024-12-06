#!/bin/bash
### tool:
# samtools/1.13
# R/4.1.2
# spp/1.16 (R package)
module clear -f
module add samtools/1.13
module add R/4.2.1

sampleID=${1}
input_bam=${2}


echo "input: $input_bam"

##########
### step1: Converting BAM to TAGALIGN files
# To converting BAM to TAGALIGN files for SPP
# input: *.bam (non-clonal and uniq mapped reads)
# output: *.tagAlign.gz

# use pre-filtered BAM files (after removing unmapped, low quality reads, multimapping reads and duplicates)
# TagAlign format is text-based and contains six tab delimited columns:
## 1. Name of the chromosome
## 2. The starting position of the feature in the chromosome. zero-based(start from zero).
## 3. The ending position of the feature in the chromosome or scaffold.
## 4. read sequence (You can just put the letter "N" in here to saves space)
## 5. score (indicates uniqueness or quality; preferably 1000)
## 6. strand
# Reference: https://github.com/kundajelab/phantompeakqualtools

echo -e "$sampleID BAM_to_TAGALIGN"
if [ -f ./$sampleID.tagAlign.gz ]; then
  echo "$sampleID.tagAlign.gz exist"
else
  echo "run $sampleID.tagAlign.gz"
  samtools view ./$input_bam | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > ./$sampleID.tagAlign.gz
  echo "$sampleID.tagAlign.gz done"
fi

# # samtools view: used to view BAM file
# # awk: used to process or extract data
# # BEGIN: before any input lines are read
# # OFS: Output Field Separator
# # \t: tab
# # and($2,16): when $2==16(equal), print 16; when $2!=16(not equal), print 0
# ## using {and (A,B)}, A and B should be number
# # $2 in BAM file (SAM format) is SAM flag
# # SAM flag 16 means "read reverse strand"
# # $3 in BAM is chromosome name
# # $4 in BAM is one-based (start from one) mapped position
# # $10 in BAM is read sequence
# # $4-1+length($10): ending position
# # use "N" for read sequence to saves space
# # use "1000" for score
# # meaning of this awk command line: if read reverse strand, then col6 of tagAlign is "-"; if not, then col6 is "+"
# # gzip -c: compress file to *.gz
# echo -e "$sampleID BAM_to_TAGALIGN_done"
# echo ""
# echo ""
# echo ""



##########
### step2: run SPP
# NOTICE: SPP is only for NARROW peak (https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)
# run SPP to predict the fragment length based on strand cross-correlation analysis
# tools: spp/1.16, R/4.1.2
# input: *.tagAlign.gz
# output: *.qc (peakshift/phantomPeak results, col3 is predicted fragment length) and *.pdf (cross-correlation plot)
# predicted fragment length is used for next step, MACS2 peak calling "--extsize"

# R version should be 3.1 or higher for SPP
# run_spp.R (SPP we used) is in spp/1.14 downloaded from "https://github.com/hms-dbmi/spp" and different from Phantompeakqualtools
# It is EXTREMELY important to filter out multi-mapping reads from the BAM/tagAlign files. Large number of multimapping reads can severly affect the phantom peak coefficient and peak calling results.
# Input file can be BAM or tagAlign. But You MUST have samtools installed to use run_spp.R with BAM files.
# *.qc contains 11 tab delimited columns:
##  1. tagAlign/BAM filename
##  2. total read number in input file
##  3. estFragLen(Estimated fragment length); take the top (first) value
##  4. corr_estFragLen
##  5. phantomPeak (Read length/phantom peak strand shift)
##  6. corr_phantomPeak (Correlation value at phantom peak)
##  7. strand shift at which cross-correlation is lowest
##  8. minimum value of cross-correlation
##  9. NSC: Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8
## 10. RSC: Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)
## 11. Quality tag based on thresholded RSC (codes= -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)
## Reference 1: https://github.com/kundajelab/phantompeakqualtools
## Reference 2: ENCODE3 pipeline v1 specifications

if [ -d ./SPP ]; then
  echo ""
else
  mkdir ./SPP 1>/dev/null 2>/dev/null
  chmod 770 ./SPP
fi
cd ./SPP

echo -e "$sampleID SPP"
if [ -f $sampleID.spp ]; then
  echo "$sampleID SPP already finished"
else
  time Rscript /bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/run_spp.r -c="../$sampleID.tagAlign.gz" -savp="$sampleID.cc.plot.pdf" -out="tmp.spp" -x=-500:100

  echo -e "Filename\tnumReads\testFragLen\tcorr_estFragLen\tPhantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tphantomPeakCoef\trelPhantomPeakCoef\tQualityTag" >> ./tmp.spp
  perl /bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/transpose.pl ./tmp.spp >> ./$sampleID.spp   

  # Extract the fragment size estimate from the SPP output file
  #   Normalized strand cross-correlation coefficient (NSC) = col9 in outFile
  #   Relative strand cross-correlation coefficient (RSC) = col10 in outFile
  #   Estimated fragment length = col3 in outFile, take the top value
  #   Important columns highlighted, but all/whole file can be stored for display

  extsize=$(cat "./$sampleID.spp" | awk 'FNR == 3 {print}' | awk '{print $1}')
  echo "fragment size estimate from the SPP: $extsize"
  
  NSC=$(cat "./$sampleID.spp" | awk 'FNR == 9 {print}' | awk '{print $1}')
  echo "Normalized strand cross-correlation coefficient (NSC) from the SPP: $NSC"
  
  RSC=$(cat "./$sampleID.spp" | awk 'FNR == 10 {print}' | awk '{print $1}')
  echo "Relative strand cross-correlation coefficient (RSC) from the SPP: $RSC"
  
fi
# Rscript: used to run R script
# -c: ChIP_alignFile, tagAlign/BAM file (can be gzipped)
# -savp: save cross-correlation plot
# -out: append peakshift/phantomPeak results to a file
# -x=<min>:<max>: strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10), the number decided according to ENCODE3 pipeline
# echo -e "$sampleID SPP_done"
# echo ""
# echo ""
# echo ""
# output file *.cc.qc: col3 is predicted fragment length
# predicted fragment length is used for next step, MACS2 peak calling "--extsize"

echo "$sampleID SPP finished"