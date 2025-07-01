#!/bin/bash

module clear -f
module add deepTools/3.5.4



# date: 20250624
# author YCW

## following by MACS3.peakcalling.v.0.1.sh to generate the heatmap for short(NFR), mono-, di-, tri- around TSS.
sampleID=${1}
inputpath=${2}
ref=${3}

#ref="/nas/cluster_data/homes/yenching/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/geneID_convert_table/Gm.a6.v1.ISU.v2.1.a4.v1.gene.bed"
cd $inputpath

echo "Reference being used:"
echo -e $ref


#bws=$(find ./MACS3/ -name "*digested*.bw")
mapfile -t bws < <(find ./MACS3/ -name "$sampleID*digested*.bw")

#echo "$datapath"
#echo "$bws"

groups=()
new_order=("short." "mono." "di." "tri.")
reorder_file_list=()
for name in "${new_order[@]}"; do
  for file in "${bws[@]}"; do
    if echo "$file" | grep -q "digested_$name"; then
      groups+=("$name")
      reorder_file_list+=("$file")
    fi
  done
done

echo ""
echo "final ${groups[@]}"
echo "final ${reorder_file_list[@]}"



output_filename="$sampleID.deeptool.heatmap"

if [ -f $output_filename.output.gz ]; then
echo "$output_filename.output.gz already exists"
else
echo "$output_filename.output.gz generating"
echo "method: scale-regions"
time computeMatrix reference-point \
    --referencePoint TSS \
    -p max -b 1000 -a 1000 -bs 1 \
    --samplesLabel ${groups[@]} \
    --skipZeros \
    -S ${reorder_file_list[@]}  \
    -R $ref \
    -o $output_filename.output.gz 
    
fi



echo ""  
echo "deeptool profile plot start"
echo "" 


plotHeatmap -m $output_filename.output.gz \
-o $output_filename.heatmap.output.pdf \
--sortUsing mean \
--zMax 25 \
--zMin 0 \
--colorMap RdYlBu_r \
--missingDataColor 0.5 \
--legendLocation best \
--heatmapHeight 12 \
--heatmapWidth 4

plotHeatmap -m $output_filename.output.gz \
-o $output_filename.heatmap.output.perGroup.pdf \
--sortUsing mean \
--perGroup \
--zMax 25 \
--zMin 0 \
--colorMap RdYlBu_r \
--missingDataColor 0.5 \
--legendLocation best \
--heatmapHeight 12 \
--heatmapWidth 4

# --zMax 3 \
# --zMin 0 \
#--perGroup \
#--samplesLabel ${groups[@]} \
#--outFileSortedRegions $output_filename.heatmap.output.bed



echo "job finished"
