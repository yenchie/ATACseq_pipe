#!/bin/bash
# Program: ATACseq_DataProcessing
#SBATCH -J ATACseq_DataProcessing
#SBATCH -o ATACseq_DataProcessing_%A.out
#SBATCH -n 20
#SBATCH --mem=40GB

## not less than 12 THREADS
set -euo pipefail

echo "$(scontrol show job $SLURM_JOBID)" | grep "JobId="
echo "$(scontrol show job $SLURM_JOBID)" | grep "Command="

# Description: This script is for ATACseq data processing and peakcalling, including ATACseqQC.
# For Pair-End data
# For reference: GmaxISU.v2.1, Gm724
# ATAC-seq histone Data Processing - from raw data to declonal, and remove mt, choloplast. 
# also estimate insert size, fastQC.
# Peakcalling
# convert to biwig
# do data processing iterratively 
#########################################################################
# Author: YCW
# Date: 2023/10/30
# Update: 2024/07/22
# Updata log: 
# 2. 2024/07/22
## add config to store the parameters and dependency function path
#
# 1.
## 2023/01/02: revised file names.
#
# Dependence:
#module add R/4.2.1
module add FastQC/11.0.2
module add Trimmomatic/0.39
module add samtools/1.13
module add sambamba/0.8.2
module add bowtie2/2.5.1
#module add picard/2.18.11
module add bedtools2/2.30.0
# module add HMMRATAC/1.2.10
#module add deepTools/3.5.4
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
#		sbatch ATACseq_DataProcessing.v1.0.sh $1 $2 $3
#
# Input:
# - $1: the path of sample list: list sample names and raw data paths you'd like to process in this run.
# formate:
#  cotH3 /archive/JYL/cotH3.L1.R1.fastq.gz,/archive/JYL/cotH3.L2.R1.fastq.gz /archive/JYL/cotH3.L1.R2.fastq.gz,/archive/JYL/cotH3.L2.R1.fastq.gz 
#  (sample_name[space]/R1fastq_path[space]/R2fastq_path, if need to pool raw data split them by comma w/o space)
# 
# - $2: assign the path of output data.
#
# - $3: config file path
#
# Output:
# - Folders named by sampleID will be created in the output folder ($2)

################################################################################


#=======================================================================#
# define major variant
list=$(cat $1)
outfolder=$2
config=$3
totalcount=$(cat $1 | wc -l)
njob=6 #define num of sample max-load
source $config
#=======================================================================#

#########################################################################
### Method
#########################################################################
### step0: create job log file to trace work progress and test whether files exist
#########################################################################

## config:
echo -e "config path: $config"

echo -e "## check these config before analysis"
echo -e "#########################################################################"
echo -e "# ***trim : MINLEN:$MINLEN"
echo -e "#########################################################################"

# create a job record file
echo > $outfolder/jobr.log
echo > $outfolder/jobf.log
# mkdir $outfolder/summary 1>/dev/null 2>/dev/null
# chmod 770 $outfolder/summary

#check environment clean
jobr=$(($(cat $outfolder/jobr.log |wc -l)-1))
jobf=$(($(cat $outfolder/jobf.log |wc -l)-1))
echo -e "running job: $jobr , finished job: $jobf"
echo ""
# total sample counts
echo -e "How many sample to do: $totalcount"
echo ""
#test file exist function
 all_exist () {
    local filename
    for filename; do
        test -f "$filename" && continue
        echo "$0: $filename does not exist - aborting" >&2
        return 1
    done
    return 0
}

#########################################################################
### step1: create the main function to do data input
#########################################################################
#main function
function mappp(){

  data=$data_l
  sampleID=$(echo $data | cut -f1 -d$' ')
  echo $sampleID >> $outfolder/jobr.log
  echo "### order num: $cc,  ID=$sampleID  ###"
  dataL=$(echo $data | cut -f2 -d' ' | sed -e "s/,/ /g")
  dataR=$(echo $data | cut -f3 -d' ' | sed -e "s/,/ /g")
  check_sampleID=$(test $sampleID && echo "1" || echo "0")

  all_exist $dataL $dataR && echo "$0: all rawdata files exist" >&2
  #all_exist $dataR && echo "$0: all R2 files exist" >&2

  
  echo -e "===================================================================="
  echo -e "[Start Analysis] SAMPLE ID: $sampleID"
  echo -e "===================================================================="

  # checkpoint: distinguish whether sample id is empty.
  # if not empty ==> create a temp. folder for the sample named by the sample name under the output directory
  if [ $check_sampleID == "0" ]; then
          echo -e "[ERROR]       Empty Sample ID!\n"
          exit
    
  else
  
    if [ -f $outfolder/$sampleID/HMMRATAC/$sampleID"_summits.bed" ]; then
      echo -e "$sampleID already Done"
      echo $sampleID >> $outfolder/jobf.log
      exit
      
    else
      if [ -d $outfolder/$sampleID ]; then
          echo ""
      else
          mkdir $outfolder/$sampleID 1>/dev/null 2>/dev/null
          chmod 770 $outfolder/$sampleID
      fi
      cd $outfolder/$sampleID
          
      if [ -f $outfolder/$sampleID/$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam ]; then
        echo -e "$sampleID Dataprocessing Done"
        echo -e "" >> ./${sampleID}.summary.txt
        
      else

          if [ -f $sampleID.trim.nochloro.nomt.sortn.bam ];then 
              echo -e "$sampleID mapping done and index exist"
              echo -e "" >> ./${sampleID}.summary.txt

          else
          
              echo -e "* LOADING paired-end data: $(basename $dataL)"
              echo -e "* LOADING paired-end data: $(basename $dataR)" 
              
              echo > ./${sampleID}.summary.txt
              #########################################################################
              ### step2: pool FASTQ
              #########################################################################
              ### NOTICE: It is an additional step, used in situations below:
              ## (1) there are two runs for one sample. There is not enough reads so there are two runs, in this situation we need to pool.
              ## (2) replicates raw data pool
              # combine the files to make a pooled file
              # input: *.fastq (raw data)
              # oustageut: *.fastq (pooled raw data)
              #########################################################################
              
              echo -e "$sampleID pooling"
              echo ""
              echo ""

              if [[ $dataL == *.gz ]]; then
                  # for gzipped fastq file
                  echo "The file is a gzipped fastq file."
                  
                  time zcat $dataL > ./$sampleID.R1.pool.fq
                  echo -e "$sampleID R1 pool done"
                  echo ""
                  echo ""
                  time zcat $dataR > ./$sampleID.R2.pool.fq
                  echo -e "$sampleID R2 pool done"
                  
              else
                  # for non-gzipped fastq file or other types
                  echo "The file is not a gzipped fastq file."
                  
                  time cat $dataL > ./$sampleID.R1.pool.fq
                  echo -e "$sampleID R1 pool done"
                  echo ""
                  echo ""
                  time cat $dataR > ./$sampleID.R2.pool.fq
                  echo -e "$sampleID R2 pool done"
              fi

              
              fastqc ./$sampleID.R1.pool.fq -o ./
              fastqc ./$sampleID.R2.pool.fq -o ./
              echo ""
              echo ""
              echo ""
              # * represent any character
              # OR cat file1 file2 > pooled_file #cat: row bind
              
              
              #########################################################################
              ### step3: count raw reads number in FASTQ
              #########################################################################
              # input: *.pool.fq (raw data)
              # oustageut: read number
              # FASTQ format: each read in FASTQ consists of four lines:
              # (1) begins with a '@' and contains sequence information
              # (2) sequence
              # (3) begins with a '+' and contains sequence information
              # (4) Phred quality score for the sequence
              # Reference: htstages://en.wikipedia.org/wiki/FASTQ_format
              #########################################################################
              echo ""
              fr=$(($(cat ${sampleID}.R1.pool.fq | wc -l) /4))
              rr=$(($(cat ${sampleID}.R2.pool.fq | wc -l)/4))
              
              
              echo -e "${sampleID} forward raw_reads_number $fr" >> ./${sampleID}.summary.txt
              echo ""
              echo -e "${sampleID} reverse raw_reads_number $rr" >> ./${sampleID}.summary.txt
              # It must be "zcat" for compressed *.gz file and "cat" for uncompressed file
              # wc: print newline, word, and byte counts for each file
              # wc -l: count lines
              # /4: divide four, because one read consist four lines
              # bc: basic calculate
              echo ""
              echo ""
              echo ""
              
              #########################################################################
              ### step4: Trim FASTQ
              #########################################################################
              # To trim and crop FASTQ and remove adapters
              # tool: Trimmomatic/0.39
              # input: *.pool.fq (raw data)
              # oustageut: *.trim.fq (good quality data after filtering and trimming)
              # Input and oustageut file can be FASTQ or *.fastq.gz, *.fastq = *.fq
              # Next step, bowtie2 input could not be *.gz, we set Trimmomatic oustageut uncompress
              #########################################################################
              echo ""
              echo "$sampleID Trim"
              #adapter=/bcst/JYL/JYL_qnap/Reference/Trimmomatic/Gm/adapter/NexteraPE-PE.fa
              #adapter=/bcst/JYL/JYL_qnap/Reference/Trimmomatic/Gm/adapter/TruSeq3-PE-2.fa
              #adapter=/bcst/JYL/JYL_qnap/Reference/Trimmomatic/Gm/adapter/ActiveMotif.PE.fa
              #adapter=/bcst/JYL/JYL_qnap/Reference/Trimmomatic/Gm/adapter/Nextera.active.motif.PE.tn5.fa
              echo -e "adapter: $adapter" >> ./${sampleID}.summary.txt
              time java -jar /software/shared/apps/Trimmomatic/0.39/trimmomatic-0.39.jar \
                PE $sampleID.R1.pool.fq $sampleID.R2.pool.fq \
                $sampleID.R1.trim.fq $sampleID.R1.trimunpaired.fq $sampleID.R2.trim.fq $sampleID.R2.trimunpaired.fq \
                ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:$MINLEN
                  
                    # LOG cutadapt ##########
                    # cutadapt -a $adapter -j 4 -e 0.25 -O 10 -o out.trimmed.R1.fq -p out.trimmed.R2.fq \
                    #     --pair-filter=any --untrimmed-output out.untrimmed.R1.fq \
                    #     --untrimmed-paired-output out.untrimmed.R2.fq GB0003_1_6_NS_R1_001.fastq GB0003_1_6_NS_R2_001.fastq > report_pe.txt  
                    
                    #adapter_err_rate=0.1 
                    #(0.1 by default)
                    # cutadapt -m $MINLEN -e $adapter_err_rate -O 10 \
                    #     -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG \
                    #     -A CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG \
                    #     -o $sampleID.R1.trim.fq -p $sampleID.R2.trim.fq \
                    #     --pair-filter=any \
                    #     --untrimmed-output $sampleID.R1.trimunpaired.fq \
                    #     --untrimmed-paired-output $sampleID.R2.trimunpaired.fq \
                    #     $sampleID.R1.pool.fq $sampleID.R2.pool.fq > report_pe.txt
                    
                    # cutadapt -m $MINLEN -e $adapter_err_rate -O 10 \
                    #     -a CTGTCTCTTATACACATCT \
                    #     -A CTGTCTCTTATACACATCT \
                    #     -o $sampleID.R1.trim.fq -p $sampleID.R2.trim.fq \
                    #     --pair-filter=any \
                    #     --untrimmed-output $sampleID.R1.trimunpaired.fq \
                    #     --untrimmed-paired-output $sampleID.R2.trimunpaired.fq \
                    #     $sampleID.R1.pool.fq $sampleID.R2.pool.fq > report_pe.txt
                    # command: java -jar path_to_trimmomatic-0.38.jar input oustageut
                    # PE: input is PAIR-END data, two oustageut files
                    # LEADING: Cut bases off the start of a read, if below a threshold quality
                    # TRAILING: Cut bases off the end of a read, if below a threshold quality
                    # SLIDINGWINDOW:<windowSize>:<requiredQuality>, Performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold.
                    # MINLEN: Drop the read if it is below a specified length; our reads are 150bp *please adjust this for you data
                    # -threads: number of threads
                    # -trimlog: creates a log of all read trimmings, including: (1) read name, (2) surviving sequence length, (3) location of the first surviving base, (4) location of the last surviving base, (5) amount trimmed from the end. NOTICE: It is a large file.
                    ##########
                    
              echo "$sampleID Trim_done"
              echo ""
              echo ""
              echo ""
              # ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
              ## ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
              ## example of ILLUMINACLIP from Po-Xing Zheng: /software/shared/apps/Trimmomatic/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10

              
              
              fastqc $sampleID.R1.trim.fq -o ./
              fastqc $sampleID.R2.trim.fq -o ./
              adapter_contentR1=$(xmllint --html --xpath 'string(//a[contains(@href, "#M10") and .="Adapter Content"]/preceding-sibling::img/@alt)'  $sampleID.R1.trim_fastqc.html)
              adapter_contentR2=$(xmllint --html --xpath 'string(//a[contains(@href, "#M10") and .="Adapter Content"]/preceding-sibling::img/@alt)'  $sampleID.R2.trim_fastqc.html)

                  if [[ $adapter_contentR1 == "[FAIL]" || $adapter_contentR2 == "[FAIL]" ]]; then
                      echo -e "adapter trimming in $sampleID.R1.trim.fq : $adapter_contentR1"
                      echo -e "adapter trimming in $sampleID.R2.trim.fq : $adapter_contentR2"
                      echo -e "adapter trimming incompletely, please check and consider re-trim, job quit."
                      exit
                  else
                      echo -e "adapter trimming in $sampleID.R1.trim.fq : $adapter_contentR1"
                      echo -e "adapter trimming in $sampleID.R2.trim.fq : $adapter_contentR2"
                      echo "$sampleID Trimed completely"
                  fi
              
              #########################################################################
              ### step5: count good reads number in FASTQ
              #########################################################################
              # count reads in FASTQ after trimming
              # input: *.trim.fq (good quality data)
              # oustageut: read number
              #########################################################################
              echo ""
              fgr=$(($(cat $sampleID.R1.trim.fq | wc -l)/4))
              rgr=$(($(cat $sampleID.R2.trim.fq | wc -l)/4))
              
              echo -e "$sampleID forward good_reads_number $fgr" >> ./$sampleID.summary.txt
              echo ""
              echo -e "$sampleID reverse good_reads_number $rgr" >> ./$sampleID.summary.txt
              # use cat for uncompress file
              echo ""
              echo ""
              echo ""
              
              #########################################################################
              ### step6: index reference genome
              #########################################################################
              ### NOTICE: This step is NOT needed except we use new reference genome.
              # to build a Bowtie2 index from reference genome
              # tool: bowtie2/2.3.2
              # input: *.fa (FASTA format, reference genome)
              # oustageut: six index files including *.1.bt2, *.2.bt2, *.3.bt2, *.4.bt2, *.rev.1.bt2, *.rev.2.bt2
              #########################################################################
              
              #echo "Index"
              #module add bowtie2/2.5.1
              #time bowtie2-build /bcst/JYL/JYL_qnap/db_2023/Gm/GmaxWm82ISU_01/v2.1/assembly/GmaxWm82ISU_01_724_v2.0.fa /bcst/JYL/JYL_qnap/db_2023/Gm/GmaxWm82ISU_01/v2.1/Gmax_724_bowtie2index
              #Usage: bowtie2-build [options] <reference_in> <bt2_base>
              # --threads: number of threads, default is 1
              # Reference: htstage://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
              #echo "Index_done"
              #echo ""
              #echo ""
              #echo ""
              
              #########################################################################
              ### step7: mapping
              #########################################################################
              # To align sequencing reads to reference genome.
              # Run bowtie2 with maximum 1 reported alignment, using 12 threads
              # tool: bowtie2/2.5.1
              # input: *.trim.fq (good quality data)
              # output: *.bam (mapping result) and alignment summary record in $sampleID.summary.txt
              
              # SAM format is text-based, tab separated, usually have header starting with '@' and contain 12 columns:
              ##  1. Read Name
              ##  2. SAM flag
              ##  3. contig name (chr name) or * for unmapped
              ##  4. mapped position of base 1 of a read on the reference sequence
              ##  5. MAPQ mapping quality
              ##  6. CIGAR string describing insertions and deletions
              ##  7. Name of reference sequence where mate’s alignment occurs
              ##  8. Position of mate
              ##  9. Template length
              ## 10. Read Sequence
              ## 11. Read Quality
              ## 12. Additional information in TAG:TYPE:VALUE format
              
              # Alignment summary for pair-end data looks like:
              # ############################################
              # 98084089 reads; of these:
              #   98084089 (100.00%) were paired; of these:
              #     51881443 (52.89%) aligned concordantly 0 times
              #     46202646 (47.11%) aligned concordantly exactly 1 time
              #     0 (0.00%) aligned concordantly >1 times
              #     ----
              #     51881443 pairs aligned concordantly 0 times; of these:
              #       7728297 (14.90%) aligned discordantly 1 time
              #     ----
              #     44153146 pairs aligned 0 times concordantly or discordantly; of these:
              #       88306292 mates make up the pairs; of these:
              #         82596351 (93.53%) aligned 0 times
              #         1732777 (1.96%) aligned exactly 1 time
              #         3977164 (4.50%) aligned >1 times
              # 57.90% overall alignment rate
              # ############################################
              ## Reference: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
              echo ""
              echo -e "$sampleID Matching"
              
              ############# chloroplast removing ###############
                echo -e "$sampleID chloroplast removing"    
                
                time bowtie2 -p 12 -x $chloroplast_bowtie2index \
                  -1 $sampleID.R1.trim.fq \
                  -2 $sampleID.R2.trim.fq 2>> ./$sampleID.summary.txt \
                  | samtools sort -n -@ 12 -T temp - > $sampleID.trim.mark.chloro.sortn.bam
  
                samtools view -c --threads 12 $sampleID.trim.mark.chloro.sortn.bam
  
                #### chroloplast mapped reads ####
                echo -e "chroloplast mapped reads"
                samtools view -b -F 4 $sampleID.trim.mark.chloro.sortn.bam  \
                      >  $sampleID.trim.sortn.chloro.bam
                samtools flagstat --threads 12 $sampleID.trim.sortn.chloro.bam >> ./$sampleID.summary.txt  
                echo -e "$sampleID.trim.sortn.chloro.bam detail of flag\n" >> ./$sampleID.summary.txt  
                      
                #### non-chroloplast reads ####     
                echo -e "non-chroloplast reads"
                samtools view -b -f 4 $sampleID.trim.mark.chloro.sortn.bam \
                      >  $sampleID.trim.sortn.nochloro.bam 
                samtools flagstat --threads 12 $sampleID.trim.sortn.nochloro.bam >> ./$sampleID.summary.txt   
                echo -e "$sampleID.trim.sortn.nochloro.bam detail of flag\n" >> ./$sampleID.summary.txt 
                     
                 
                ### split to pair end fastq ### 
                time bamToFastq -i $sampleID.trim.sortn.nochloro.bam \
                  -fq $sampleID.trim.nochloro.R1.fq \
                  -fq2 $sampleID.trim.nochloro.R2.fq 2>/dev/null
                  
                dcfr=$(($(cat $sampleID.trim.nochloro.R1.fq | wc -l) /4))
                dcrr=$(($(cat $sampleID.trim.nochloro.R2.fq | wc -l) /4))
                
                
                chr=$(samtools view -c --threads 12 "$sampleID.trim.sortn.chloro.bam")
                echo -e "\n$sampleID chloroplast reads number $chr\n" >> ./$sampleID.summary.txt  
                echo -e "\n$sampleID dechloroplast reads number R1: $dcfr R2: $dcrr\n" >> ./$sampleID.summary.txt  
                echo -e "$sampleID chloroplast removed"    
              
              
              ############# mitochondria removing ###############
                echo -e "$sampleID mitochondria removing"   
                time bowtie2 -p 12 -x $mitochondria_bowtie2index \
                  -1 $sampleID.trim.nochloro.R1.fq \
                  -2 $sampleID.trim.nochloro.R2.fq 2>> ./$sampleID.summary.txt \
                  | samtools sort -n -@ 12 -T temp - > $sampleID.trim.nochloro.mark.mt.sortn.bam
                  
                samtools view -c --threads 12 $sampleID.trim.nochloro.mark.mt.sortn.bam
                
                #### mitochondria mapped reads ####
                echo -e "mitochondria mapped reads"
                samtools view -b -F 4 $sampleID.trim.nochloro.mark.mt.sortn.bam \
                      >  $sampleID.trim.nochloro.sortn.mt.bam
                samtools flagstat --threads 12 $sampleID.trim.nochloro.sortn.mt.bam >> ./$sampleID.summary.txt  
                echo -e "$sampleID.trim.nochloro.sortn.mt.bam detail of flag\n" >> ./$sampleID.summary.txt  
                      
                #### non-mitochondria reads ####    
                echo -e "non-mitochondria reads"
                samtools view -b -f 4 $sampleID.trim.nochloro.mark.mt.sortn.bam \
                      >  $sampleID.trim.nochloro.sortn.nomt.bam 
                samtools flagstat --threads 12 $sampleID.trim.nochloro.sortn.nomt.bam >> ./$sampleID.summary.txt     
                echo -e "$sampleID.trim.nochloro.sortn.nomt.bam detail of flag\n" >> ./$sampleID.summary.txt  
                
                 
                ### split to pair end fastq ### 
                time bamToFastq -i $sampleID.trim.nochloro.sortn.nomt.bam \
                  -fq $sampleID.trim.nochloro.nomt.R1.fq \
                  -fq2 $sampleID.trim.nochloro.nomt.R2.fq 2>/dev/null
                  
                dcmfr=$(($(cat $sampleID.trim.nochloro.nomt.R1.fq | wc -l) /4))
                dcmrr=$(($(cat $sampleID.trim.nochloro.nomt.R2.fq | wc -l) /4))
                
                
                mtr=$(samtools view -c --threads 12 "$sampleID.trim.nochloro.sortn.mt.bam")
                echo -e "\n$sampleID mitochondria reads number, after dechloroplast $mtr\n" >> ./$sampleID.summary.txt  
                echo -e "\n$sampleID dechloroplast, demitochondria reads number R1: $dcmfr R2: $dcmrr\n" >> ./$sampleID.summary.txt  
                echo -e "$sampleID mitochondria removed"    
              
              ############# mapping ###############
              echo -e "$sampleID Matching: $bowtie2index"
              time bowtie2 --no-unal -p 12 -x $bowtie2index \
                -1 $sampleID.trim.nochloro.nomt.R1.fq \
                -2 $sampleID.trim.nochloro.nomt.R2.fq  2>> ./$sampleID.summary.txt \
                | samtools view -@ 12 -h -F 0x08 | samtools sort -n -@ 12 -T temp  - > $sampleID.trim.nochloro.nomt.sortn.bam
        
              # -F 0x08: remove unmapped read (0x4)
              # -f 0x08: remove mapped read (0x4)
              # -F 0x08: remove singletons (mate unmapped (0x8))
              
              rm ./*mark*.bam
              
              # NOTICE: this command line is for Pair-END data and the resulte is sorted by name.
              # --no-unal: Suppress SAM records for reads that failed to align. (unmapped reads are not in oustageut SAM file)
              # -k: Maximum number of alignments to report per read. Default is 1, and if multiple equivalent alignments exist, it chooses one randomly.
              # -k 2: to identify read is truely uniq or not in next step (see next step)
              # -x: The basename of the index for the reference genome.
              # Gmax_724_index: index from Gmax_189.fa, reference genome for Glycine max
              # -1, -2: input file
              # -S: oustageut SAM file
              # -N: Sets the number of mismatches. Can be set to 0 or 1. Setting this higher makes alignment much slower but increases sensitivity. Default: 0.
              
              echo -e "$sampleID matching done"
              echo ""
              echo ""
              echo ""
              #########################################################################
              ### step8: count mapped reads
              #########################################################################
              # count mapped reads number in bam
              # tool: samtools/1.13
              # input: *.bam (mapping result)
              # output: read number
              # samtools view can be used to view and convert SAM (as well as *.sam.gz)/BAM/CRAM files
              # Reference: http://www.htslib.org/doc/samtools-view.html
              echo ""
              #samtools view -c --threads 12 $sampleID.trim.sam
              mr=$(samtools view -c --threads 12 "$sampleID.trim.nochloro.nomt.sortn.bam")
              echo -e "\n$sampleID mappable reads number $mr\n" >> ./$sampleID.summary.txt
              echo ""
              samtools flagstat --threads 12 $sampleID.trim.nochloro.nomt.sortn.bam >> ./$sampleID.summary.txt
              echo -e "\n$sampleID.trim.nochloro.nomt.sortn.bam detail of flag\n" >> ./$sampleID.summary.txt
              # samtools view: here is used to view bam file
              echo ""
              
            
          fi

          #########################################################################
          ### step9:  Compute library complexity
          #########################################################################
          # Sort by name
          # convert to bedPE and obtain fragment coordinates
          # sort by position and strand
          # Obtain unique count statistics
          if [ -f ./$sampleID.trim.nochloro.nomt.pbc.qc ]; then
                echo "" >> ./$sampleID.summary.txt
          else
                echo -e "TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair" > tmp.pbc.qc
                bedtools bamtobed -bedpe -i $sampleID.trim.nochloro.nomt.sortn.bam 2>/dev/null | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' \
                  | sort | uniq -c \
                  | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
                  >> tmp.pbc.qc
                
                echo -e "\n$sampleID.trim.nochloro.nomt.sortn.bam detail of flag\n" >> ./$sampleID.trim.nochloro.nomt.pbc.qc
                samtools flagstat $sampleID.trim.nochloro.nomt.sortn.bam >> ./$sampleID.trim.nochloro.nomt.pbc.qc
                echo "" >> $sampleID.trim.nochloro.nomt.pbc.qc
                perl $transpose tmp.pbc.qc >> $sampleID.trim.nochloro.nomt.pbc.qc
                echo "" >> $sampleID.trim.nochloro.nomt.pbc.qc
                echo "The non-redundant fraction (NRF) is the fraction of non-redundant mapped reads \n
                      in a dataset; it is the ratio between the number of positions in the genome \n
                      that uniquely mapped reads map to and the total number of uniquely mappable \n
                      reads. \n
                      The NRF should be > 0.8. \n
                      The PBC1 is the ratio of genomic locations with EXACTLY one read pair over the genomic locations with AT LEAST one read \n
                      pair. PBC1 is the primary measure, and the PBC1 should be close to 1. \n
                      Provisionally 0-0.5 is severe bottlenecking, 0.5-0.8 is moderate bottlenecking, \n
                      0.8-0.9 is mild bottlenecking, and 0.9-1.0 is no bottlenecking. \n
                      The PBC2 is the ratio of genomic locations with EXACTLY one read pair over the genomic \n
                      locations with EXACTLY two read pairs. The PBC2 should be significantly \n
                      greater than 1." >> $sampleID.trim.nochloro.nomt.pbc.qc
                      
                rm ./tmp.pbc.qc
          fi
          
          echo ""
          echo ""
        
          #########################################################################
          ### step10:  declonal and filter low quality (MAPQ<=30) reads
          #########################################################################
          # input=$sampleID.trim.nochloro.nomt.sortn.bam
          # output: 
          #       - $sampleID.trim.nochloro.nomt.sortn.nonclonal.bam
          #       - $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam **
          #       - $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam.bai **
          #  ** final output sort by position  
          
          if [ -f ./$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam.bai ]; then
              echo "$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam exist and index exist" >> ./$sampleID.summary.txt
          else
              #########################################################################
              ### step10-1:  declonal
              #########################################################################
              # -r, --remove-duplicates: remove duplicates instead of just marking them
              # -t, --nthreads=NTHREADS
              # --overflow-list-size=OVERFLOW_LIST_SIZE
              #            size of the overflow list where reads, thrown from the hash table,
              #            get a second chance to meet their pairs (default is 200000 reads);
              #            increasing the size reduces the number of temporary files created
              
              sambamba markdup --overflow-list-size 600000 -r -t 12 $sampleID.trim.nochloro.nomt.sortn.bam $sampleID.trim.nochloro.nomt.sortn.nonclonal.bam
              echo -e "$sampleID Remove duplicates (remove clonal reads)"
              
              # count non_clonal reads
              dec=$(samtools view -c "$sampleID.trim.nochloro.nomt.sortn.nonclonal.bam")
              echo -e "\n$sampleID #reads after rmdup (non-clonal reads) $dec\n" >> ./$sampleID.summary.txt
              samtools flagstat $sampleID.trim.nochloro.nomt.sortn.nonclonal.bam >> ./$sampleID.summary.txt
              echo -e "\n$sampleID.trim.nochloro.nomt.sortn.nonclonal.bam detail of flag\n" >> ./$sampleID.summary.txt
              echo ""
              echo ""
              
              #########################################################################
              ### step10-1:  filter low quality (MAPQ<=$mapq) reads
              #########################################################################
              # USAGE:
              # -h: Include the header in the output.
              # -t 12: Utilize 12 threads for processing.
              # -q 10: filter out the low quality reads lower than 30. 
          
              samtools view \
              -h \
              -t 12 \
              -q $mapq \
              -o $sampleID.trim.nochloro.nomt.nonclonal.unique.bam \
              $sampleID.trim.nochloro.nomt.sortn.nonclonal.bam
              
              echo -e "$sampleID Remove low quality (MAPQ<=$mapq) reads and multimatch (multiple locations) after remove duplicates"
              echo ""
              fir=$(samtools view -c "$sampleID.trim.nochloro.nomt.nonclonal.unique.bam")
              echo -e "\nsampleID #unique reads (unique location) $fir\n" >> ./$sampleID.summary.txt
              
              samtools sort -@ 12 -T temp $sampleID.trim.nochloro.nomt.nonclonal.unique.bam -o $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
          
              samtools index $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
              echo ""
              samtools flagstat --threads 12 $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam >> ./$sampleID.summary.txt
              echo -e "\n$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam detail of flag\n" >> ./$sampleID.summary.txt
              echo ""
              echo ""
          fi
        
          
          # Check if any .fq files exist
          echo "remove processing tmp data"
          if ls ./*.fq 1> /dev/null 2>&1; then
              # Remove all .fq files
              rm ./*.fq
              echo "remove processing tmp data, .fq files."
          else
              echo "No processing tmp data, .fq files found."
          fi

          if ls ./*trim.nochloro.nomt.nonclonal.unique.bam 1> /dev/null 2>&1; then
              # Remove all .fq files
              rm ./*trim.nochloro.nomt.nonclonal.unique.bam
              echo "remove processing tmp data, *trim.nochloro.nomt.nonclonal.unique.bam file."
          else
              echo "No processing tmp data, *trim.nochloro.nomt.nonclonal.unique.bam file found."
          fi

           
          echo -e "$sampleID data recorded in ./$sampleID.summary.txt"
          
          
          #########################################################################
          ### step11: estimate insert size
          #########################################################################
          # input=$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
          # output: 
          #       - ${input}.stats.txt
          #       - bam.stat.plot
          #       - $sampleID_multiqc_report_data
          #       - ${input}.coverage
          #       - insert_size_metric
          #       - $sampleID.ATACseq.QC.pdf
        
          samtools stats $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam \
            > $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.stats.txt
            
          mkdir ./bam.stat.plot 
          plot-bamstats -p ./bam.stat.plot/ $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.stats.txt 
          insert_size=$(grep "insert size " $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.stats.txt)
          echo -e "\n$insert_size\n" >> ./$sampleID.summary.txt
          
          /software/shared/apps/Python/3.7.15/bin/multiqc ./$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.stats.txt \
              -i $sampleID
      
          samtools coverage $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam -m \
              -o $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.coverage
            
          bash $CollectInsertSizeMetrics ./$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam         
      
          echo ""
          echo ""
          echo -e "$sampleID data processing Done"
                      
      fi

        #########################################################################
        ### step12: Shift read coordinates
        #########################################################################
        # input=$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
        # output: 
        #       - $sampleID.trim.nochloro.nomt.nonclonal.unique.shifted.sort.bw
        # Shift read coordinates
        # An optional step in analysing data generated using the Tn5 transposase (such as ATAC-seq, ChIPmentation etc.) is to account for a small DNA insertion, 
        # introducted as repair of the transposase-induced nick introduces a 9bp insertion. 
        # Reads aligning to the + strand should be offset by +4bp and reads aligned to the -ve strand should be offset by -5bp. 
        # For references, see the first ATAC-seq paper by Buenrostro et al., (2013) and the analysis by Adey et al., (2010) which showed this insertion bias. 
        # Shifting coordinates is only really important if single-base resolution is required, for example in the analysis of transcription factor motifs in ATAC-seq peak footprints. 
        # Be aware that some tools do this shifting themselves (so double check manuals!).
        # We can use the deeptools package.
        
        if [ -f ./$sampleID.trim.nochloro.nomt.nonclonal.unique.shifted.sort.bw ]; then
            echo "" >> ./$sampleID.summary.txt
            echo "$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.shifted.bam file and index exist"
        else
            bash $Shift_bam $sampleID $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
        fi
        
        #########################################################################
        ### step13: SPP
        #########################################################################
        # input=$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
        # output: 
        #       - lm30109A.tagAlign.gz
        #       - ./SPP
        
        if [ -f ./SPP/$sampleID.spp ]; then
            echo "" >> ./$sampleID.summary.txt
            echo "$sampleID spp already done"
        else
          echo "spp script path: $SPP"
          bash $SPP $sampleID $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
            # Extract the fragment size estimate from the SPP output file
            #   Normalized strand cross-correlation coefficient (NSC) = col9 in outFile
            #   Relative strand cross-correlation coefficient (RSC) = col10 in outFile
            #   Estimated fragment length = col3 in outFile, take the top value
            #   Important columns highlighted, but all/whole file can be stored for display
    
            extsize=$(cat "./SPP/$sampleID.spp" | awk 'FNR == 3 {print}' | awk '{print $1}')
            echo "fragment size estimate from the SPP: $extsize" >> ./$sampleID.summary.txt
            
            NSC=$(cat "./SPP/$sampleID.spp" | awk 'FNR == 9 {print}' | awk '{print $1}')
            echo "Normalized strand cross-correlation coefficient (NSC) from the SPP: $NSC" >> ./$sampleID.summary.txt
            
            RSC=$(cat "./SPP/$sampleID.spp" | awk 'FNR == 10 {print}' | awk '{print $1}')
            echo "Relative strand cross-correlation coefficient (RSC) from the SPP: $RSC" >> ./$sampleID.summary.txt
            echo "$sampleID spp done"
        fi
        
        if [ -f ./SPP/tmp.spp ]; then
          rm ./SPP/tmp.spp
          echo "spp tmp data removed"
        else
          echo ""
        fi
        
        
        #########################################################################
        ### step14: peakcaling
        #########################################################################
        # input=$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
        # output: 
        #       - $sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam, properly pairs
        #       - ./MACS2
        #       - ./HMMRATAC

        echo -e "get  properly pairs "
        samtools view -@ 12 -h -f 2 $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam | samtools sort -@ 12 -T temp -o $sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam -
        samtools index $sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam
        pcr=$(samtools view -c "$sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam")
        echo -e "\n$sampleID peak calling input reads: $pcr\n" >> ./$sampleID.summary.txt
        
        ## peak calling by MACS2 
          echo "peak calling tool: MACS2, path: $MACS2"
          bash $MACS2 $sampleID $sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam 2> MACS2.peakcalling.log
          
          RiP=$(cat ./MACS2/$sampleID-readCountInPeaks.txt.summary | awk '$1 == "Assigned" { sum += $2 } END { print sum }')
          total=$(cat ./MACS2/$sampleID-readCountInPeaks.txt.summary | awk 'BEGIN { sum = 0 } NR > 1 { sum += $2 } END { print sum }')
  
          FRiP=$(echo "scale=4; $RiP / $total" | bc)
          echo "peakcalling RiP by MACS2: $RiP" >> ./$sampleID.summary.txt
          echo "peakcalling total input alignments by MACS2: $total" >> ./$sampleID.summary.txt
          echo "peakcalling FRiP by MACS2: $FRiP" >> ./$sampleID.summary.txt

        # ## peak calling by MACS3
        # echo "peak calling tool: MACS3, path: $MACS3"
        # bash $MACS3 $sampleID $sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam 2> MACS3.peakcalling.log
        
        # RiP=$(cat ./MACS3/$sampleID-readCountInPeaks.txt.summary | awk '$1 == "Assigned" { sum += $2 } END { print sum }')
        # total=$(cat ./MACS3/$sampleID-readCountInPeaks.txt.summary | awk 'BEGIN { sum = 0 } NR > 1 { sum += $2 } END { print sum }')

        # FRiP=$(echo "scale=4; $RiP / $total" | bc)
        # echo "peakcalling RiP by MACS3: $RiP" >> ./$sampleID.summary.txt
        # echo "peakcalling total input alignments by MACS3: $total" >> ./$sampleID.summary.txt
        # echo "peakcalling FRiP by MACS3: $FRiP" >> ./$sampleID.summary.txt
          
          

                    ## retired --------------------------------  
                    # ## peak calling by HMMRATAC 
                    #   # The summit is 1bp, so subtract 50 from the start and add 50 to the stop for newStart and newStop and use that resulting bed file as the input for DiffBind
              
                    #   echo "peak calling tool: HMMRATAC, path: $HMMRATAC"
                    #   bash $HMMRATAC $sampleID $sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam 2> HMMRATAC.peakcalling.log

                    #   #bash /bcst/JYL/JYL_qnap_2/YCWang/0_Script/Project/Gm/ATACseq/tool/HMMRATAC.peakcalling.v.1.sh $sampleID.f2 $sampleID.trim.nochloro.nomt.nonclonal.unique.f2.sort.bam 2> HMMRATAC.peakcalling.2.log

                    #   # Filter HMMRATAC output by the score, if desired.
                    #   # Score threshold will depend on dataset, score type and user preference. A threshold of 10 would be:
                    #   #
                    #   # awk -v OFS="\t" '$13>=10 {print}' $sampleID"_peaks.gappedPeak" > $sampleID.filteredPeaks.gappedPeak
                    #   #
                    #   # # To filter the summit file by the same threshold:
                    #   #
                    #   # awk -v OFS="\t" '$5>=10 {print}' $sampleID"_summits.bed" > $sampleID.filteredSummits.bed
                    #   #
                    #   # NOTE: HMMRATAC will report all peaks that match the structure defined by the model, including weak peaks. 
                    #   #       Filtering by score can be used to retain stronger peaks. 
                    #   #       Lower score = higher sensitivity and lower precision, 
                    #   #       Higher score = lower sensitivity and higher precision.
                      
                    #   RiP=$(cat ./HMMRATAC/$sampleID-readCountInPeaks.txt.summary | awk '$1 == "Assigned" { sum += $2 } END { print sum }')
                    #   total=$(cat ./HMMRATAC/$sampleID-readCountInPeaks.txt.summary | awk 'BEGIN { sum = 0 } NR > 1 { sum += $2 } END { print sum }')
              
                    #   FRiP=$(echo "scale=4; $RiP / $total" | bc)
                    #   echo "peakcalling RiP by HMMRATAC: $RiP" >> ./$sampleID.summary.txt
                    #   echo "peakcalling total input reads by HMMRATAC: $total" >> ./$sampleID.summary.txt
                    #   echo "peakcalling FRiP by HMMRATAC: $FRiP" >> ./$sampleID.summary.txt
                    
        
        
      echo ""
      echo ""
      echo ""
      echo -e "$sampleID peakcalling Done"
      echo $sampleID >> $outfolder/jobf.log
      
    fi   
  fi
}





########### here will keep submit data to mapp till the floor of max($totalcount , $jobf) ########## 
echo "Data processing and peakcalling start"
declare -i cc=1

#check job record log
echo ""
#################################################################
#step9: Using until-loop loading sample and start Data processing
#################################################################
until  [ $jobf -eq $totalcount ]
do
  jobr=$(($(cat $outfolder/jobr.log |wc -l)-1))
  jobf=$(($(cat $outfolder/jobf.log |wc -l)-1))
  jobc=$(expr $jobr + 1 - $jobf)
 
  if [ $cc -le $totalcount ] && [ $jobc -le $njob ]; then
     #echo $cc
     h=$(($cc))
     data_l=$(cat $1 | head -n $h | tail -n +$h)
     echo ""
     #echo $data_l

     # excute main function #
     mappp $data_l &
      cc=$(($cc+1))
      #echo $cc
     sleep 10s
   else
   echo "please waiting, fully load"
     sleep 600s
  fi
  jobf=$(($(cat $outfolder/jobf.log |wc -l)-1))

done 
echo "ALL processing finish"