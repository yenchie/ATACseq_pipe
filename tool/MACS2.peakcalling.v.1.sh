#!/bin/bash
module clear -f
module add MACS2/2.2.9.1
module add samtools/1.13
module add UCSC/v369
module add subread/2.0.3

sampleID=${1}
input=${2}
e_genome_size=1.01e+9

### peakcaling ### 


if [ -f ./$sampleID.genome.info ]; then
  echo ""
else
  samtools view -H $input \
  | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
  > $sampleID.genome.info
fi

if [ -d ./MACS2 ]; then
  echo ""
else
  mkdir ./MACS2
fi

prefix=$(basename "$input" .bam)
echo "input: $prefix"

macs2 randsample -i $input -f BAMPE -p 100 -o ./MACS2/$prefix.bed

# -p PERCENTAGE, --percentage PERCENTAGE
#                         Percentage of tags you want to keep. Input 80.0 for
#                         80%. This option can't be used at the same time with
#                         -n/--num. REQUIRED


#Call peaks
macs2 callpeak -g $e_genome_size -f BEDPE --nomodel --shift -37 --extsize 73 \
  -B --broad --broad-cutoff 0.01 --keep-dup all --cutoff-analysis -n $sampleID -t ./MACS2/$prefix.bed \
  --outdir ./MACS2 2> macs2.log


## from https://github.com/CebolaLab/ATAC-seq
# If using paired-end reads, MACS2 will be used with the -f BEDPE option. Depending on the analysis aims, there are several different options that can be used. The ENCODE3 pipeline uses the --nomodel --shift -37 --extsize 73 options for analysing ATAC-seq data, to account for the size of nucleosomes. Nucleosomes cover ~145 bp and the ATAC-seq reads need to be shifted towards the 5' end by half this distance.


# The following commands require an input file detailing the chromosome sizes. 
# Use the UCSC tool fetchChromSizes (install via conda): fetchChromSizes hg38 > hg38.chrom.sizes. 
# Generate the fold-change bedGraph
# bash ~/fetchChromSizes.sh hg38 > hg38.chrom.sizes

bedGraphToBigWig ./MACS2/$sampleID"_treat_pileup.bdg" $sampleID.genome.info ./MACS2/$sampleID"_treat_pileup.bw"

# FE bw
macs2 bdgcmp -t ./MACS2/$sampleID"_treat_pileup.bdg" -c ./MACS2/$sampleID"_control_lambda.bdg" -m FE -o ./MACS2/$sampleID"_FE.bdg" 
sort -k1,1 -k2,2n ./MACS2/$sampleID"_FE.bdg" > ./MACS2/$sampleID"_FE.sorted.bdg"
bedGraphToBigWig ./MACS2/$sampleID"_FE.sorted.bdg" $sampleID.genome.info ./MACS2/$sampleID"_macs2_FE.bw"

rm ./MACS2/$sampleID"_FE.bdg"
rm ./MACS2/$sampleID"_FE.sorted.bdg"

# logLR bw
macs2 bdgcmp -t ./MACS2/$sampleID"_treat_pileup.bdg" -c ./MACS2/$sampleID"_control_lambda.bdg" -m logLR -p 0.00001 -o ./MACS2/$sampleID"_logLR.bdg" 
sort -k1,1 -k2,2n ./MACS2/$sampleID"_logLR.bdg" > ./MACS2/$sampleID"_logLR.sorted.bdg"
bedGraphToBigWig ./MACS2/$sampleID"_logLR.sorted.bdg" $sampleID.genome.info ./MACS2/$sampleID"_macs2_logLR.bw"

rm ./MACS2/$sampleID"_logLR.bdg"
rm ./MACS2/$sampleID"_logLR.sorted.bdg"

# usage: macs2 bdgcmp [-h] -t TFILE -c CFILE [-S SFACTOR] [-p PSEUDOCOUNT]
#                    [-m {ppois,qpois,subtract,logFE,FE,logLR,slogLR,max} [{ppois,qpois,subtract,logFE,FE,logLR,slogLR,max} ...]] [--outdir OUTDIR]
## from https://github.com/macs3-project/MACS/wiki/Build-Signal-Track
# -m FE means to calculate fold enrichment. Other options can be logLR for log likelihood, subtract for subtracting noise from treatment sample.
# -p sets pseudocount. This number will be added to 'pileup per million reads' value. You don't need it while generating fold enrichment track because control lambda will always >0. But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.



#############################################
# QC - fraction of reads in peak (FRiP score)

# One quality metric for peak calling is to calculate the fraction of reads in peak (FRiP) score. 
# For ATAC-seq, the FRiP score is recommended to be >0.2, with >0.3 as optimal. 
# We will use featureCounts from the SourceForge Subread package. 
# Install Subread using conda install -c bioconda subread (or see this link to install from the source).

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ./MACS2/$sampleID"_peaks.broadPeak" \
  > ./MACS2/$sampleID"_peaks.saf"

### count
echo "$input"
featureCounts -p --countReadPairs -a ./MACS2/$sampleID"_peaks.saf" -F SAF -o ./MACS2/$sampleID-readCountInPeaks.txt $input

# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 
# -p                  Specify that input data contain paired-end reads. To
#                      perform fragment counting (ie. counting read pairs), the
#                      '--countReadPairs' parameter should also be specified in
#                      addition to this parameter.
# --countReadPairs    Count read pairs (fragments) instead of reads. This option
#                      is only applicable for paired-end reads.


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr ./MACS2/$sampleID"_peaks.broadPeak" > ./MACS2/$sampleID"_peaks_sorted.broadPeak"

echo -e "MACS2 peakcalling done, $sampleID"


