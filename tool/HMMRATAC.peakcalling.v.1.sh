#!/bin/bash
module clear -f
module add HMMRATAC/1.2.10
module add UCSC/v369
module add subread/2.0.3

sampleID=${1}
input=${2}
### peakcaling ### 
# The summit is 1bp, so subtract 50 from the start and add 50 to the stop for newStart and newStop and use that resulting bed file as the input for DiffBind
echo "HMMRATAC peakcalling ..."
#$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam
samtools view -H $input \
  | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
  > $sampleID.genome.info

if [ -d ./HMMRATAC ]; then
  echo ""
else
  mkdir ./HMMRATAC
fi

#

time java -jar /software/shared/apps/HMMRATAC/1.2.10/HMMRATAC_V1.2.10_exe.jar \
    -b $input \
    -i $input.bai \
    -g $sampleID.genome.info \
    -o ./HMMRATAC/$sampleID \
    --bedgraph True --peaks True --score max

grep E0 ./HMMRATAC/$sampleID.bedgraph > ./HMMRATAC/$sampleID.State0_regions.bed
grep E1 ./HMMRATAC/$sampleID.bedgraph > ./HMMRATAC/$sampleID.State1_regions.bed
grep E2 ./HMMRATAC/$sampleID.bedgraph > ./HMMRATAC/$sampleID.State2_regions.bed


# Filter HMMRATAC output by the score, if desired.
# Score threshold will depend on dataset, score type and user preference. A threshold of 10 would be:

#awk -v OFS="\t" '$13>=10 {print}' ./HMMRATAC/$sampleID"_peaks.gappedPeak" > ./HMMRATAC/$sampleID.filteredPeaks.gappedPeak

# To filter the summit file by the same threshold:

#awk -v OFS="\t" '$5>=10 {print}' ./HMMRATAC/$sampleID"_summits.bed" > ./HMMRATAC/$sampleID.filteredSummits.bed

# NOTE: HMMRATAC will report all peaks that match the structure defined by the model, including weak peaks. 
#       Filtering by score can be used to retain stronger peaks. 
#       Lower score = higher sensitivity and lower precision, 
#       Higher score = lower sensitivity and higher precision.

###################
# QC - fraction of reads in peak (FRiP score)

# One quality metric for peak calling is to calculate the fraction of reads in peak (FRiP) score. 
# For ATAC-seq, the FRiP score is recommended to be >0.2, with >0.3 as optimal. 
# We will use featureCounts from the SourceForge Subread package. 
# Install Subread using conda install -c bioconda subread (or see this link to install from the source).

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ./HMMRATAC/$sampleID"_peaks.gappedPeak" \
  > ./HMMRATAC/$sampleID"_peaks.saf"

### count
featureCounts -p -a ./HMMRATAC/$sampleID"_peaks.saf" -F SAF -o ./HMMRATAC/$sampleID-readCountInPeaks.txt $input

###################
## USAGE ####
#   java -jar HMMRATAC_V1.2_exe.jar -b <SortedBAM> -i <BAMIndex> -g <GenomeStatsFile> <options>
#
#   -m , --means <Double>
#                     Comma separated list of initial mean values for the fragment
#                     distribution, used to create the signal tracks.  The default values
#                     are 50,200,400,600.  It is recommended to change these values if
#                     you are working with non-human species, as other species have
#                     different nucleosome spacing.  These values will be updated with EM
#                     training, unless the -f option is set to false.  Also note that the
#                     first value, corresponding to the short read distribution, is NOT
#                     updated with EM training and will remain as it is set with this
#                     option.  It is recommended to set this first value to the read
#                     length, if the read length is under 100bp.
#   -s , --stddev <Double>
#                     Comma separated list of initial standard deviations values for the
#                     fragment distributions., used to create the signal tracks.  The
#                     default values are 20,20,20,20.  These values may also need to be
#                     updated for non-human species, although it may not be necessary, as
#                     the EM is robust in handling variances.  These values are also
#                     updated with EM training.  Note: the first value is not updated, as
#                     the short distribution uses an Exponential distribution that only
#                     uses the mean value, so the first value is meaningless.
#   -f , --fragem <True || False>
#                     Boolean to determine whether fragment EM training is to occur or
#                     not.  The default is true.  If this option is set to false, the
#                     initial values described with -m and -s are used as the
#                     parameters of the mixture model.  This is generally not
#                     recommended.  Setting this value to false can decrease the total
#                     runtime, but may result in a less accurate model.  If the data has
#                     been run already, and it is desired to re-run it using different
#                     reporting thresholds, you could set this option to false, provided
#                     that you reset the initial parameters to the updated ones created
#                     by the previous run.  These values are recorded in the .log file
#                     described later.
#   -q , --minmapq <int>
#                     This is the minimum mapping quality score for reads to be used in
#                     creating the signal tracks.  The default is 30.  Although this can
#                     be set to higher or lower values, it is generally not recommended.
#                     If few meet this threshold, that would indicate that the assay
#                     itself was compromised, and setting a lower value would likely
#                     still result in errors.
#   -u , --upper <int>
#                     Upper limit on fold change range for choosing training regions.
#                     HMMRATAC chooses training regions by finding genomic loci whose fold
#                     change above genomic background is within a certain range.  This
#                     option sets the upper limit of this range.  Default is 20.
#                     Generally speaking, the higher this range is, the more stringent
#                     the resulting model, and the lower the range, the less stringent
#                     the model.  If recall (sensitivity) is important for you, it is
#                     recommended to set lower ranges, while if precision or accuracy is
#                     more important, it is recommended to set higher ranges.  Once the
#                     regions within a certain range are found, HMMRATAC extends the regions
#                     by +/- 5KB and uses these extended regions to train the model.  The
#                     search ends once a maximum of 1000 regions are identified.
#   -l , --lower <int>
#                     Lower limit on fold change range for choosing
#                     training regions.  HMMRATAC chooses training regions by finding genomic
#                     loci whose fold change above genomic background is within a certain
#                     range.  This option sets the lower limit of this range.  Default is
#                     10.  Generally speaking, the higher this range is, the more
#                     stringent the resulting model, and the lower the range, the less
#                     stringent the model.  If recall (sensitivity) is important for you,
#                     it is recommended to set lower ranges, while if precision or
#                     accuracy is more important, it is recommended to set higher ranges.
#                     Once the regions within a certain range are found, HMMRATAC extends
#                     regions by +- 5KB and uses these extended regions to train the
#                     model.  The search ends once a maximum of 1000 regions are
#                     identified.
#   -z , --zscore <int>
#                     Zscored read depth to mask during Viterbi decoding.  Default is
#                     100.  In our experience, Viterbi has trouble decoding regions that
#                     Have very high read coverage.  If encountered, Viterbi will either
#                     call a Empty block, where no state is called, or will call one
#                     large open region.  To avoid this problem, HMMRATAC will skip over any
#                     regions whose centered Read coverage is equal to or greater than
#                     this value.  If the report peaks Option is set (-p True – described
#                     below), these regions are added back to The peak file and are
#                     labeled as “High_Coverage_Peak_#”.  These regions Can therefore be
#                     examined further, if desired.
#   -o , --output <Name>
#                     The prefix for all output files.  Default is “NA”.
#   -e , --blacklist <BED>
#                     BED file of regions to exclude from model creation and genome
#                     decoding.  These could be previously annotated blacklist regions.
#                     It is always recommended to include a blacklisted region list.
#   -p , --peaks <True || False>
#                     Boolean to determine whether or not to report a peak file in BED
#                     format.  Default is True.  The resulting file will be in gappedPeak
#                     format and will be described in detail below.
#   -k , --kmeans <int>
#                     Number of states in the model.  Default is 3.  It is generally not
#                     recommended to change this setting.  All of our tests and
#                     comparisons were done using the default setting.  It may be
#                     reasonable to set to a lower number of states when using 500 cell
#                     per-replicate data, but this is untested.  If this setting is
#                     changed, it is generally recommended to not report peaks, but only
#                     to report the genome-wide state annotation file.  This can then be
#                     interpreted manually.
#   -t , --training <BED>
#                     BED file of training regions to use instead of using fold-change
#                     ranges.  If this option is used, it is recommended to only use 1000
#                     regions and to extend them by +/- 5KB, as HMMRATAC would do for
#                     fold-change identified training regions.
#   --bedgraph <True || False>
#                     Boolean to determine whether a genome-wide state annotation
#                     bedgraph file should be reported.  Default is False.
#   --minlen <int>
#                     Minimum length of an open region state to call a peak.  Note: that
#                     the -p option must be set to true. Default is 200.
#   --score <max || ave || med || fc || zscore || all>
#                     What type of score system to use for scoring peaks.  “Max” refers
#                     to the maximum number of reads mapping to the open region.  “Med”
#                     refers to the median number of reads mapping to the open region.
#                     “Ave” refers to the average number of reads mapping to the open
#                     region.  “FC” refers to fold-change and is the average number of
#                     reads mapping to the open region divided by the genome average.
#                     “Zscore” refers to the average number of reads mapping to the open
#                     region minus the genome average divided by the genomic standard
#                     deviation.  “All” reports all these scores separated by a “_”
#                     (underscore). The order for "all" is MAX_MEAN_MEDIAN_ZSCORE_FC
#                     Default is “max”.
#   --bgscore <True || False>
#                     Boolean to determine whether to add a score to every state
#                     annotation in a bedgraph file.  Note that --bedgraph has to be set
#                     true.  This adds considerable time to the program and is generally
#                     not recommended. Default is false.
#   --trim <int>
#                     How many signals, or distributions, to trim from the signal tracks.
#                     It trims from the end of the matrix (IE 1 means trim the tri signal
#                     track, 2 means trim the tri  and di signal tracks…etc).  This could
#                     be useful if your data  doesn’t contain many longer fragments, such
#                     as  size selected  data.  Recommendations:  fragments <=  500bp set
#                     --trim 1; fragments <= 250bp set --trim 2; fragments <=
#                     150 set --trim 3.  Default is 0.
#   --window <int>
#                     Size of the bins to split the genome into for Viterbi decoding.  To
#                     save memory, HMMRATAC splits the genome into these sized bins and
#                     Performs Viterbi decoding on these bins separately.  Default is
#                     25000000.  It may be necessary to reduce the size of these bins, if
#                     running HMMRATAC on A desktop or a machine with limited memory.  Most
#                     desktops should handle A bin size of 1/10 the default.  A java heap
#                     space runtime error, is common When the bin size is too large for
#                     the machine (see [section 4](#troubleshooting)).
#   --model <File>
#                     This model (binary model generated by previous HMMRATAC run, suffixed
#                     with .model) will be used to decode the genome rather than building
#                     a new model using training regions (either provided by user with
#                      the -t option or created by HMMRATAC using the -u and -l
#                     options).  This can be used in conjunction with the --modelonly
#                     option (create the model, inspect it and run HMMRATAC with that model).
#   --modelonly <True || False>
#                     Boolean to determine if HMMRATAC should quit after generating the
#                     model.  This is a helpful option when trying to determine the best
#                     parameters for creating the model.  Default = false.
#   --maxTrain <int>
#                     Maximum number of training regions to use during model building.
#                     It is possible that the default number of training regions (1000) is 
#                     too many and will cause memory issues for smaller machines. Setting 
#                     this parameter to 1/2 of the default (ie. 500) could help with this 
#                     issue.
#
################################################
