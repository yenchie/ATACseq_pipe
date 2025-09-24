#!/bin/bash
module clear -f
module add MACS2/2.2.9.1
module add samtools/1.13
module add UCSC/v369
module add subread/2.0.3

sampleID=${1}
input=${2}
geneome_info=${3}


macs2 callpeak -t /bcst/JYL/JYL_qnap_2/YCWang/Project/Gm/ATACseq/Peakcalling/HMMRATAC/bam/lm.pool.ATAC/lm.pool.ATAC.merge.sort.bam --keep-dup all -n TEMP_MACS --nolambda -B

macs2 callpeak -t /bcst/JYL/JYL_qnap_2/YCWang/Project/Gm/ATACseq/Peakcalling/HMMRATAC/bam/mm.pool.ATAC/mm.pool.ATAC.merge.sort.bam --keep-dup all -n TEMP_MACS --nolambda -B

bedGraphToBigWig TEMP_MACS_treat_pileup.bdg /bcst/JYL/JYL_qnap_2/YCWang/Project/Gm/ATACseq/Peakcalling/HMMRATAC/lm.pool.ATAC.merge/lm.pool.ATAC.merge.genome.info lm.pool.ATAC.merge_fromMACS2.bw

bedGraphToBigWig TEMP_MACS_treat_pileup.bdg /bcst/JYL/JYL_qnap_2/YCWang/Project/Gm/ATACseq/Peakcalling/HMMRATAC/mm.pool.ATAC.merge/mm.pool.ATAC.merge.genome.info mm.pool.ATAC.merge_fromMACS2.bw



rm TEMP_MACS*



