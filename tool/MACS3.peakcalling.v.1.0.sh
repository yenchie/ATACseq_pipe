#!/bin/bash
module clear -f
module load Python/3.10.8
module add samtools/1.13
module add UCSC/v369
module add subread/2.0.3

python3 -m venv ~/envs/MACS3env/
# Activate your environment
source ~/envs/MACS3env/bin/activate

sampleID=${1}
input=${2}
#e_genome_size=1.01e+9

### peakcaling ### 


if [ -f ./$sampleID.genome.info ]; then
  echo ""
else
  samtools view -H $input \
  | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
  > $sampleID.genome.info
fi


if [ -d ./MACS3 ]; then
  echo ""
else
  mkdir ./MACS3
fi

prefix=$(basename "$input" .bam)
echo "input: $prefix"

## run cutoff-analysis ------
macs3 hmmratac \
  --cutoff-analysis-only \
  -i $input \
  -f BAMPE \
  --outdir ./MACS3 \
  -n $sampleID.poisson.ATAC \
  --hmm-type poisson 2> macs3.log



# Run MACS3 HMMRATAC
macs3 hmmratac \
  -i $input \
  -f BAMPE \
  --outdir ./MACS3 \
  -n $sampleID.poisson.ATAC \
  --save-digested \
  --save-states \
  --save-training-data \
  --hmm-type poisson \
  -c 6 2> macs3.log


## generate bigwig -----------
# Loop through all digested bdg files
for bdg in $(find ./MACS3/ -name "*digested*.bdg"); do
    bdg_prefix=$(basename "$bdg" .bdg)  # remove path and .bdg extension
    bedGraphToBigWig "$bdg" "$sampleID.genome.info" "./MACS3/${bdg_prefix}.bw"
    echo -e "removing $bdg"
    rm $bdg
    echo -e "removed $bdg"
done



#############################################
# QC - fraction of reads in peak (FRiP score)

# One quality metric for peak calling is to calculate the fraction of reads in peak (FRiP) score. 
# For ATAC-seq, the FRiP score is recommended to be >0.2, with >0.3 as optimal. 
# We will use featureCounts from the SourceForge Subread package. 
# Install Subread using conda install -c bioconda subread (or see this link to install from the source).

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ./MACS3/$sampleID".poisson.ATAC_accessible_regions.narrowPeak" \
  > ./MACS3/$sampleID"_peaks.saf"

### count
echo "$input"
featureCounts -p --countReadPairs -a ./MACS3/$sampleID"_peaks.saf" -F SAF -o ./MACS3/$sampleID-readCountInPeaks.txt $input

# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 
# -p                  Specify that input data contain paired-end reads. To
#                      perform fragment counting (ie. counting read pairs), the
#                      '--countReadPairs' parameter should also be specified in
#                      addition to this parameter.
# --countReadPairs    Count read pairs (fragments) instead of reads. This option
#                      is only applicable for paired-end reads.

RiP=$(cat ./MACS3/$sampleID-readCountInPeaks.txt.summary | awk '$1 == "Assigned" { sum += $2 } END { print sum }')
          total=$(cat ./MACS3/$sampleID-readCountInPeaks.txt.summary | awk 'BEGIN { sum = 0 } NR > 1 { sum += $2 } END { print sum }')
  
          FRiP=$(echo "scale=4; $RiP / $total" | bc)

          
echo "peakcalling RiP by MACS3: $RiP" >> "./macs3.log"
echo "peakcalling total input alignments by MACS3: $total" >> "./macs3.log"
echo "peakcalling FRiP by MACS3: $FRiP" >> "./macs3.log"


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k5,5nr "./MACS3/"$sampleID".poisson.ATAC_accessible_regions.narrowPeak" > "./MACS3/"$sampleID".poisson.ATAC_accessible_regions.sorted.narrowPeak"

## generate split digested bw profile --------

time bash /bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/deeptools.heatmap.4.macs3.splited.bws.v0.1.sh $sampleID  ./ \
          /bcst/JYL/JYL_qnap_2/YCWang/testing/macs3_test/sampled.gene.ref.bed
          #/nas/cluster_data/homes/yenching/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/geneID_convert_table/Gm.a6.v1.ISU.v2.1.a4.v1.gene.bed


echo -e "MACS3 peakcalling done, $sampleID"



