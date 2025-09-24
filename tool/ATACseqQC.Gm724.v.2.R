# ATACseqQC.Gm724.v.2.R
# author: YCW
# update: 20250624
# update log: removing reduce redundant plots.

# date: 20231207
# Description: ATACseqQC 
## estimate and visualizing insert sizes distribution and TSSEScore



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ATACseqQC")

if (!require("ChIPseeker", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ChIPseeker")
}

if (!require("ATACseqQC", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("ATACseqQC")
}

rm(list=ls())
library(dplyr)
library(stringr)
library(tidyr)
library(GenomicFeatures)
library(ATACseqQC)
library(Rsamtools)
txdb <- loadDb("/bcst/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/TxDb/Gm724")
txs <- transcripts(txdb)

args = commandArgs(trailingOnly=TRUE)
cat("
          #USAGE:Rscript --vanilla /bcst/JYL/YCWang/Rscript/ATACseqQC.Gm724.v.1.R [arg1] \n
          #Arg1: input folder path, ex. ./$sampleID.trim.nochloro.nomt.nonclonal.unique.sort.bam \n
          ")

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("1 arguments must be supplied (input).n", call.=FALSE)
}

#----------------------
  bamfile<- args[1]
  print(paste("input:", bamfile))
  if (bamfile%>%as.character()%>%basename() %>%str_detect(".trim")) {
    sampleID<- str_extract(bamfile%>%as.character()%>%basename(), ".*(?=.trim)")
  }else{
    sampleID<- gsub(".bam", "", basename(bamfile))
  }
  
  print(paste("sampleID:", sampleID))
  
  bamfile.labels <- gsub(".bam", "", basename(bamfile))
  
  pdf(paste0(sampleID, ".ATACseq.QC.pdf"))
  #estimateLibComplexity(readsDupFreq(bamfile))
  
  
  # First, there should be a large proportion of reads with less than 100 bp, which represents the nucleosome-free region. Second, the fragment size distribution should have a clear periodicity, which is evident in the inset figure, indicative of nucleosome occupancy (present in integer multiples).
  fragSize <- fragSizeDist(bamfile, bamfile.labels)
  
  # Tn5 transposase has been shown to bind as a dimer and inserts two adapters into accessible DNA locations separated by 9 bp.2
  # 
  # Therefore, for downstream analysis, such as peak-calling and footprinting, all reads in input bamfile need to be shifted. The function shiftGAlignmentsList can be used to shift the reads. By default, all reads aligning to the positive strand are offset by +4bp, and all reads aligning to the negative strand are offset by -5bp.1
  # 
  # The adjusted reads will be written into a new bamfile for peak calling or footprinting
  

  
  ## bamfile tags to be read in
  possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                  "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                  "TC", "UQ"), 
                      "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                    "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                    "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                    "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                    "U2"))
  
  bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                       param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0]
  tags
  
  
  #seqlev <- "Gm01" ## subsample data for quick run
  seqinformation <- seqinfo(txdb)
  seqlev <- seqinformation%>%as.data.frame()%>%row.names()
  which <- as(seqinformation[seqlev], "GRanges")
  param <- ScanBamParam(which=which)
  gal <-  GenomicAlignments::readGAlignments(bamfile, use.names=TRUE, param=param)

  # gal <- ATACseqQC::readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
  # shiftedBamfile <- file.path(outPath, "shifted.bam")
  # gal1 <- shiftGAlignmentsList(gal)
  
  
  ## plot -----
  
  # Transcription Start Site (TSS) Enrichment Score  -------
  # TSS enrichment score is a raio between aggregate distribution of reads centered on TSSs and that flanking the corresponding TSSs. TSS score = the depth of TSS (each 100bp window within 1000 bp each side) / the depth of end flanks (100bp each end). TSSE score = max(mean(TSS score in each window)). TSS enrichment score is calculated according to the definition at https://www.encodeproject.org/data-standards/terms/#enrichment. Transcription start site (TSS) enrichment values are dependent on the reference files used; cutoff values for high quality data are listed in the following table from https://www.encodeproject.org/atac-seq/.
  
  tsse <- TSSEscore(gal, txs)
  tsse$TSSEscore %>% print()
  plot(100*(-9:10-.5), tsse$values, type="b", 
       main=paste(sampleID, ", TSSEscore:", tsse$TSSEscore%>%round(2)),
       xlab="distance to TSS",
       ylab="aggregate TSS score")

       

  
  

  dev.off()

