#!/bin/bash
# Program: Peakcall
#SBATCH -J Peakcall
#SBATCH -o ATAC_PeakCalling_idr_%A.out
#SBATCH -n 12
#SBATCH --mem=20GB

## not less than 12 THREADS
set -euo pipefail

echo "$(scontrol show job $SLURM_JOBID)" | grep "JobId="
echo "$(scontrol show job $SLURM_JOBID)" | grep "Command="
date 
echo ""
# Description: ATAC-seq peak calling analysis, 
# For Pair-End data
# convert to biwig
# estimate idr
# do peakcalling iterratively 
#########################################################################
# Author: YCW
# Date: 2025/02/10
# Update: 
# Updata log: 
# 

#
# Dependence:
module add samtools/1.13
module add deepTools/3.5.4 
module add idr/2.0.3


# Configure Setting (if needed):
# # - all option
# #configure loading
# source $configure
# THREADS=$core
#
################################################################################
### USAGE: 
################################################################################
#		sbatch ATACseq_PeakCalling.idr.v.0.1.sh $1 $2
#
# Input:
# - $1: the path of sample list: list sample names and raw data paths you'd like to process in this run.
# sample_list.csv formate:
# Stage(Group),rep1,rep1.bam.path,rep2,bam.path
##  Note:
# 
# - $2: assign the path of output data.
#
# Output:
# - Folders named by sampleID will be created in the output folder ($2)

################################################################################


#=======================================================================#
# define major variant
list=$(cat $1)
outfolder=$2
totalcount=$(cat $1 | wc -l)
njob=6 #define num of sample max-load
ref="Gm724"
peakcall_function="/bcst/JYL/JYL_qnap_2/YCWang/0_Script/000/Pipeline/ATAC_pipe/tool/HMMRATAC.peakcalling.v.1.sh"
#=======================================================================#
# config:

if [ "$ref" == "Gm724" ]; then
  e_genome_size=1010968777
elif [ "$ref" == "Gm189" ]; then
  e_genome_size=955054837
else
  echo "unknown reference efficiency size"
fi

echo -e "## check these config before analysis"
echo -e "#########################################################################"
echo -e "# ***peakcall : effective genome size: $e_genome_size, $ref"
echo -e "#########################################################################"

# create a job record file
echo >$outfolder/jobr.log
echo >$outfolder/jobf.log

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


# EXPT=""
# NAME=""
# INPUT=""
# peak_type=""
function rep_idr(){
    echo $EXPT
    echo $NAME
    echo $peak_type
    
    if [ -d ./$EXPT ]; then
    echo "./$EXPT exist"
    else
      mkdir ./$EXPT 1>/dev/null 2>/dev/null
      chmod 770 ./$EXPT
    fi
    
    cd ./${EXPT}
    samtools view -@ 12 -H ${NAME} > ./${EXPT}_header.sam
    nreads=$(samtools view -@ 12 -c "$NAME")
    nlines=$(( (nreads + 1) / 2 ))
    
    # This will shuffle the lines in the file and split itinto two SAM files
    samtools view -@ 12 ${NAME} | shuf - | split -d -l ${nlines} - "./${EXPT}" 
    cat ./${EXPT}_header.sam ./${EXPT}00 | samtools view -@ 12 -bS | samtools sort -@ 12 -T temp - > ./${EXPT}_00.bam
    cat ./${EXPT}_header.sam ./${EXPT}01 | samtools view -@ 12 -bS | samtools sort -@ 12 -T temp - > ./${EXPT}_01.bam

    samtools index ./${EXPT}_00.bam
    samtools index ./${EXPT}_01.bam

    rm ./${EXPT}_header.sam
    rm ./${EXPT}00
    rm ./${EXPT}01  
    
    #Peak calling on pseudoreplicates
 
    rep=$(printf "%s\n" "${EXPT}" | sed 's/N_//g')
    echo $rep

    echo "Calling peaks for pseudoreplicate0"
    sampleID="$GroupID.$rep.0"
    echo $sampleID
    peak1=./HMMRATAC/$sampleID"_peaks_sorted.$peak_type" 
    if [ -f $peak1 ]; then
      echo "$sampleID peakcalling already done"
    else
      ## peak calling by HMMRATAC 
        bash $peakcall_function $sampleID ./${EXPT}_00.bam 2> $sampleID.HMMRATAC.peakcalling.log
    fi
    
    echo "Calling peaks for pseudoreplicate1"
    sampleID="$GroupID.$rep.1"
    echo $sampleID
    peak2=./HMMRATAC/$sampleID"_peaks_sorted.$peak_type" 
    if [ -f $peak2 ]; then
      echo "$sampleID peakcalling already done"
    else
      ## peak calling by HMMRATAC 
        bash $peakcall_function $sampleID ./${EXPT}_01.bam 2> $sampleID.HMMRATAC.peakcalling.log
    fi


    #Independent replicate IDR (Nn) 
    sampleID=NULL
    echo $peak1 $peak2
    echo "Running IDR on pseudoreplicates..."
    time idr --samples $peak1 $peak2 \
      --input-file-type $peak_type --output-file ./${EXPT}_pseudorep-idr --rank score --allow-negative-scores --plot --log-output-file ./${EXPT}_idr.log
    echo $peak1 $peak2 >> ./${EXPT}_idr.log
    
    cd ../
    echo "$NAME idr finished"
    
  }
  
# Main function
function run_idr_peakcall(){
  data=$data_group
  GroupID=$(echo $data | cut -f1 -d$',')
  rep1=$(echo $data | cut -f2 -d$',')
  rep1_bam=$(echo $data | cut -f3 -d$',')
  rep2=$(echo $data | cut -f4 -d$',')
  rep2_bam=$(echo $data | cut -f5 -d$',')
  control=$(echo $data | cut -f6 -d',')
  echo "$GroupID" >> $outfolder/jobr.log
  echo "### order num: $cc,  ID=$GroupID, peak_type: ATAC ###"

  check_GroupID=$(test $GroupID && echo "1" || echo "0")
  all_exist $rep1_bam $rep2_bam $control && echo "$0: all bam files exist" >&2
  
  if [ $check_GroupID == "0" ]; then
    echo -e "[ERROR]       Empty Group ID!\n"
    exit
  else
   echo ""
  fi

  cd $outfolder
  
  peak_type="bed"
  
  if [ -d ./$GroupID ]; then
    echo "./$GroupID exist"
  else
    mkdir ./$GroupID 1>/dev/null 2>/dev/null
    chmod 770 ./$GroupID
  fi
  
  cd ./$GroupID
  
  ### Part1: rep peak calling ###
  # input:
  # 1.sampleID
  # 2.IP
  # 3.control
  # output:
  # sampleID="$GroupID.$repn"
  # ./sampleID/MACS2/$outname"_peaks.$peak_type"
  # ./sampleID/MACS2/$sampleID."_peaks_sorted.$peak_type"
  
  sampleID="$GroupID.$rep1"
  IP=$rep1_bam
  echo $sampleID
  echo $IP
  if [ -f ./HMMRATAC/$sampleID"_peaks_sorted.$peak_type" ]; then
    echo "$sampleID peakcalling already done"
  else
    ## peak calling by HMMRATAC 
          echo "peak calling tool: HMMRATAC, path: $peakcall_function"
          bash $peakcall_function $sampleID $IP 2> HMMRATAC.peakcalling.log
          echo -e "$sampleID Peakcalling Done"
  fi
  
  sampleID="$GroupID.$rep2"
  IP=$rep2_bam
  echo $sampleID
  echo $IP
  if [ -f ./HMMRATAC/$sampleID"_peaks_sorted.$peak_type" ]; then
    echo "$sampleID peakcalling already done"
  else
    ## peak calling by HMMRATAC 
        echo "peak calling tool: HMMRATAC, path: $peakcall_function"
        bash $peakcall_function $sampleID $IP 2> HMMRATAC.peakcalling.log
        echo -e "$sampleID Peakcalling Done"
  fi
  
  ########################
  ########## IDR #########
  ########################
  #Independent replicate IDR (Nt) 
  if [ -d ./IDR ]; then
    echo "./IDR exist"
  else
    mkdir ./IDR 1>/dev/null 2>/dev/null
    chmod 770 ./IDR
  fi
  
  cd ./IDR
  
  echo "Running IDR on relicates..."
  
  EXPT="Nt_$rep1.$rep2"
  if [ -d ./$EXPT ]; then
    echo "./$EXPT exist"
  else
    mkdir ./$EXPT 1>/dev/null 2>/dev/null
    chmod 770 ./$EXPT
  fi
  # [ -d "./$EXPT" ] && echo "./$EXPT exists" || { mkdir -p "./$EXPT" 1>/dev/null 2>/dev/null && chmod 770 "./$EXPT"; }
  
  if [ -f ./${EXPT}/${EXPT}-idr ]; then
    echo "./${EXPT}/${EXPT}-idr exist"
  else
    time idr --samples \
      $outfolder/$GroupID/HMMRATAC/$GroupID.$rep1"_peaks_sorted.$peak_type" \
      $outfolder/$GroupID/HMMRATAC/$GroupID.$rep2"_peaks_sorted.$peak_type" \
      --input-file-type $peak_type --output-file ./${EXPT}/${EXPT}-idr --rank score --allow-negative-scores --plot --log-output-file ./${EXPT}/${EXPT}_idr.log
      
      echo "Nt done"
  fi   
    
  ##########################  
  
  # N1, N2, self_pseudorep
  EXPT="N_$rep1"
  NAME="$rep1_bam"
  
  if [ -f ./${EXPT}/${EXPT}_pseudorep-idr ]; then
    echo "./${EXPT}/${EXPT}_pseudorep-idr exist"
  else
    rep_idr $EXPT $NAME
    echo "$EXPT done"
  fi
  
  EXPT="N_$rep2"
  NAME="$rep2_bam"
  
  if [ -f ./${EXPT}/${EXPT}_pseudorep-idr ]; then
    echo "./${EXPT}/${EXPT}_pseudorep-idr exist"
  else
    rep_idr $EXPT $NAME
    echo "$EXPT done"
  fi

  # Np, pooled_pseudorep
  #Merge treatment BAMS
  echo "Merging BAM files for pseudoreplicates..."
  EXPT="N_pooledPseudo_$rep1.$rep2"
  if [ -d ./$EXPT ]; then
    echo "./$EXPT exist"
  else
    mkdir ./$EXPT 1>/dev/null 2>/dev/null
    chmod 770 ./$EXPT
  fi
  
  
  NAME1="$rep1_bam"
  NAME2="$rep2_bam"
  
  #Merge treatment BAMS
  if [ -f ./$EXPT/$GroupID.${rep1}_${rep2}_merged.bam ]; then
    echo "Merge treatment BAM prepared"
  else
    echo "Merging BAM files for pseudoreplicates..."
    samtools merge -u ./$EXPT/$GroupID.${rep1}_${rep2}_merged.bam ${NAME1} ${NAME2}
    echo "Merge treatment BAM done"
  fi

  NAME="$outfolder/$GroupID/IDR/$EXPT/$GroupID.${rep1}_${rep2}_merged.bam"
   
  if [ -f ./${EXPT}/${EXPT}_pseudorep-idr ]; then
    echo "./${EXPT}/${EXPT}_pseudorep-idr exist"
  else
    rep_idr $EXPT $NAME
    echo "$EXPT done"
  fi
  
  echo "$GroupID $rep1, $rep2 idr done"
  
  ### summary ---
  idrs=$(find "./" -name "*-idr" -print)
  echo $idrs
 
  IDR_THRESH=0.05 
  summaryNAME="idr.$IDR_THRESH.summary"
  
      # # =============================
      # # Get peaks passing IDR threshold of 10%
      # # =============================
      # IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')
      # awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${IDR_OUTPUT} | sort | uniq | sort -k7n,7n | gzip -c > ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz
      # NPEAKS_IDR=$(zcat ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz | wc -l)

  # IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')
  # #IDR_THRESH_TRANSFORMED=$(printf "%.3f" $IDR_THRESH_TRANSFORMED)
  # echo $IDR_THRESH_TRANSFORMED
      
  #min(int(log2(-125IDR), 1000). peaks with an IDR of 0.05 have a score of int(-125log2(0.05)) = 540
  IDR_THRESH_TRANSFORMED=$(echo "-125 * l($IDR_THRESH)/l(2)" | bc -l)
  # Remove decimal part to get integer value
  #scaleIDR=${scaleIDR%.*}
  echo $IDR_THRESH_TRANSFORMED
  
    # Col 5: int(log2(-125IDR)
  if [ "$peak_type" == "broadPeak" ]; then
  
    awk 'BEGIN{OFS="\t"} $5>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' ./Nt_$rep1.$rep2/Nt_$rep1.$rep2-idr \
    | sort | uniq | sort -k7n,7n > ./$rep1.$rep2.conservative.IDR.$IDR_THRESH.$peak_type
    
    awk 'BEGIN{OFS="\t"} $5>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' ./N_pooledPseudo_$rep1.$rep2/N_pooledPseudo_${rep1}.${rep2}_pseudorep-idr \
    | sort | uniq | sort -k7n,7n > ./$rep1.$rep2.optimal.IDR.$IDR_THRESH.$peak_type

  elif [ "$peak_type" == "narrowPeak" ]; then
  
    awk 'BEGIN{OFS="\t"} $5>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ./Nt_$rep1.$rep2/Nt_$rep1.$rep2-idr \
    | sort | uniq | sort -k7n,7n > ./$rep1.$rep2.conservative.IDR.$IDR_THRESH.$peak_type
    
    awk 'BEGIN{OFS="\t"} $5>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ./N_pooledPseudo_$rep1.$rep2/N_pooledPseudo_${rep1}.${rep2}_pseudorep-idr \
      | sort | uniq | sort -k7n,7n > ./$rep1.$rep2.optimal.IDR.$IDR_THRESH.$peak_type

  elif [ "$peak_type" == "bed" ]; then
  
    awk 'BEGIN{OFS="\t"} $5>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' ./Nt_$rep1.$rep2/Nt_$rep1.$rep2-idr \
    | sort | uniq | sort -k7n,7n > ./$rep1.$rep2.conservative.IDR.$IDR_THRESH.$peak_type
    
    awk 'BEGIN{OFS="\t"} $5>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' ./N_pooledPseudo_$rep1.$rep2/N_pooledPseudo_${rep1}.${rep2}_pseudorep-idr \
      | sort | uniq | sort -k7n,7n > ./$rep1.$rep2.optimal.IDR.$IDR_THRESH.$peak_type

  else
    echo "Invalid peak type: $peak_type"
  fi
  
  
  IFS=$'\n'   # Set Internal Field Separator to newline to handle filenames with spaces
  echo "" >> ./$summaryNAME
  for idr in $idrs; do
      echo "$idr"
      
      peaks=$(awk -v IDR="$IDR_THRESH_TRANSFORMED" '{if($5 >= IDR) print $0}' "$idr" | wc -l)

      echo "$peaks"
      echo "$idr $peaks" >> ./$summaryNAME
  done
  echo "" >> ./$summaryNAME

  awk -v rep1="$rep1" -v rep2="$rep2" '
  $1 ~ "./N_" rep1 {N1=$2}
  $1 ~ "./N_" rep2 {N2=$2}
  END {
      if (N1 && N2) {
          maxN = (N1 > N2) ? N1 : N2;
          minN = (N1 > N2) ? N2 : N1;
          ratio = (N1 > N2) ? "N_" rep1 "/N_" rep2 : "N_" rep2 "/N_" rep1;
          print "self-consistency:", ratio, ":", maxN, "/", minN, "=", maxN/minN;
      }
  }' ./$summaryNAME >> ./$summaryNAME

  
  awk  -v rep1="$rep1" -v rep2="$rep2" '
      /\.\/N_pooledPseudo_'$rep1'.'$rep2'/ {Np=$2} /\.\/Nt_'$rep1'.'$rep2'/ {Nt=$2} END {
      if (Np && Nt) {
          maxN = (Np > Nt) ? Np : Nt;
          minN = (Np > Nt) ? Nt : Np;
          ratio = (Np > Nt) ? "Np/Nt" : "Nt/Np";
          print "peak consistency:", rep1, rep2, ratio, ":", maxN, "/", minN, "=", maxN/minN;
      }
  }' "$summaryNAME" >> "$summaryNAME"
  
  sort $summaryNAME | uniq | tee $summaryNAME.txt

  
  echo "$GroupID.$rep1.$rep2" >> $outfolder/jobf.log


}          

#Main run
###########here will keep submit data to mapp till the floor of max($totalcount , $jobf)########## 
echo "Peakcalling-idr start"

declare -i cc=1

#check job record log
echo ""
until  [ $jobf -eq $totalcount ]
do

  jobr=$(($(cat $outfolder/jobr.log |wc -l)-1))
  jobf=$(($(cat $outfolder/jobf.log |wc -l)-1))
  jobc=$(expr $jobr + 1 - $jobf)
 
  if [ $cc -le $totalcount ] && [ $jobc -le $njob ]; then
     #echo $cc
     h=$(($cc))
     data_group=$(cat $1 | head -n $h | tail -n +$h)
     echo ""
     #echo $data_l

     # excute main function #
     run_idr_peakcall $data_group &
      cc=$(($cc+1))
      #echo $cc
     sleep 10s
   else
   echo "please waiting, fully load"
     sleep 300s
  fi
  jobf=$(($(cat $outfolder/jobf.log |wc -l)-1))

done 

echo "Peakcalling-idr finish"