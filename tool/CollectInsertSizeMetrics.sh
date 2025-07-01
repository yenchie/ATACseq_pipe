#!/bin/bash
module add R/4.2.1

input_bam=${1}

    # if [ -f insert_size_metric_no_hist.tsv ]; then
    #   echo "$input_bam CollectInsertSizeMetrics already done"
    # else  
    #   java -jar /software/shared/apps/picard/3.1.1/build/libs/picard.jar CollectInsertSizeMetrics \
    #         INPUT=$input_bam \
    #         OUTPUT=insert_size_metric.tsv \
    #         H=insert_size_metric.pdf \
    #         2> CollectInsertSizeMetrics.log
      
    #   cat insert_size_metric.tsv \
    #           | grep -v "^#" | head -3 |  perl /bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/transpose.pl > insert_size_metric_no_hist.tsv      
    #   echo "$input_bam CollectInsertSizeMetrics done"     
    
    # fi

# if ls *.ATACseq.QC.pdf 1> /dev/null 2>&1; then
#   rm -r ./splited.bam
#   echo "Previous results of ATACseqQC removed."
# else
#   echo "Clean path ready for ATACseqQC."
# fi


time Rscript --vanilla /bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/ATACseqQC.Gm724.v.2.R $input_bam 2> ATACseqQC.log
echo "$input_bam ATACseqQC done"     
      