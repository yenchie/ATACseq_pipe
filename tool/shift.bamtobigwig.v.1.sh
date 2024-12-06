#!/bin/bash
module clear -f
module add deepTools/3.5.4
sampleID=${1}
input=${2}

echo "input: $input"
samtools index ./$input

alignmentSieve --numberOfProcessors max --ATACshift --bam ./$input \
    -o ./$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.shifted.bam
    
samtools sort -@ 12 -T temp $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.shifted.bam -o $sampleID.trim.nochloro.nomt.nonclonal.unique.shifted.sort.bam
          
samtools index ./$sampleID.trim.nochloro.nomt.nonclonal.unique.shifted.sort.bam
        
        
binSize="10"
#smoothLength="20"
#--smoothLength "$smoothLength" \
#echo "binSize=$binSize \n smoothLength=$smoothLength"
echo "binSize=$binSize"
time bamCoverage -b ./$sampleID.trim.nochloro.nomt.nonclonal.unique.shifted.sort.bam \
    -o ./$sampleID.trim.nochloro.nomt.nonclonal.unique.shifted.sort.bw \
    --binSize "$binSize" \
    --normalizeUsing RPKM \
    -p max 2> "./bw.log"

prefix=$(basename "$input" .bam)
time bamCoverage -b ./$prefix.bam \
    -o ./$prefix.bw \
    --binSize "$binSize" \
    --normalizeUsing RPKM \
    -p max 2> "./bw.log"
    
echo "$input bam file shifted and bigwig file created"
