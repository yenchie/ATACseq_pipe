#!/bin/bash
module clear -f
module add MACS2/2.2.9.1
module add samtools/1.13
module add UCSC/v369
module add subread/2.0.3

sampleID=${1}
input=${2}
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
macs2 callpeak -q 0.01 -g 7.88e8 -f BEDPE --nomodel --shift -37 --extsize 73 \
  -B --broad --keep-dup all --cutoff-analysis -n $sampleID -t ./MACS2/$prefix.bed \
  --outdir ./MACS2


# The following commands require an input file detailing the chromosome sizes. 
# Use the UCSC tool fetchChromSizes (install via conda): fetchChromSizes hg38 > hg38.chrom.sizes. 
# Generate the fold-change bedGraph
# bash ~/fetchChromSizes.sh hg38 > hg38.chrom.sizes
macs2 bdgcmp -t ./MACS2/$sampleID"_treat_pileup.bdg" -c ./MACS2/$sampleID"_control_lambda.bdg" -m qpois -o ./MACS2/$sampleID"_qpois.bdg" 

#Sort the bedGraph file and convert to bigWig
sort -k1,1 -k2,2n ./MACS2/$sampleID"_qpois.bdg" > ./MACS2/$sampleID"_qpois.sorted.bdg"

bedGraphToBigWig ./MACS2/$sampleID"_qpois.sorted.bdg" $sampleID.genome.info ./MACS2/$sampleID"_macs2_qpois.bw"

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
featureCounts -p -a ./MACS2/$sampleID"_peaks.saf" -F SAF -o ./MACS2/$sampleID-readCountInPeaks.txt $input


echo -e "MACS2 peakcalling done, $sampleID"


