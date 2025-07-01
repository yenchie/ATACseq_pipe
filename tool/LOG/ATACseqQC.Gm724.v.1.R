# ATACseqQC.Gm724.v.1.R
# author: YCW
# date: 20231207
# Description: ATACseqQC 


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
library(ChIPseeker)
library(DOSE)
library(ggplot2)
library(ATACseqQC)
library(ChIPpeakAnno)
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
  
  
  ## files will be output into outPath
  outPath <- file.path("splited.bam")
  dir.create(outPath, recursive = T, mode = "0770")
  ## shift the coordinates of 5'ends of alignments in the bam file
  
  #seqlev <- "Gm01" ## subsample data for quick run
  seqinformation <- seqinfo(txdb)
  seqlev <- seqinformation%>%as.data.frame()%>%row.names()
  which <- as(seqinformation[seqlev], "GRanges")
  gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
  shiftedBamfile <- file.path(outPath, "shifted.bam")
  gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
  
  ## plot -----
  
  
  
  # Promoter/Transcript body (PT) score  -------
  # PT score is calculated as the coverage of promoter divided by the coverage of its transcript body. 
  # PT score will show if the signal is enriched in promoters.
  
  pt <- PTscore(gal1, txs)
  plot(pt$log2meanCoverage,
       pt$PT_score,
       xlab = "log2 mean coverage",
       ylab = "Promoter vs Transcript")
  
  
  # Nucleosome Free Regions (NFR) score  -------
  # NFR score is a ratio between cut signal adjacent to TSS and that flanking the corresponding TSS. 
  # Each TSS window of 400 bp is first divided into 3 sub-regions: the most upstream 150 bp (n1), the most downstream of 150 bp (n2), and the middle 100 bp (nf). Then the number of fragments with 5’ ends overlapping each region are calculated for each TSS. The NFR score for each TSS is calculated as NFR-score = log2(nf) - log2((n1+n2)/2). A plot can be generated with the NFR scores as Y-axis and the average signals of 400 bp window as X-axis, very like a MA plot for gene expression data.
  
  nfr <- NFRscore(gal1, txs)
  plot(nfr$log2meanCoverage, nfr$NFR_score, 
       main=paste(sampleID, "NFRscore for 200bp flanking TSSs"),
       xlab="log2 mean coverage",
       ylab="Nucleosome Free Regions score",
       xlim=c(-10, 0), ylim=c(-5, 5))
  
  # Transcription Start Site (TSS) Enrichment Score  -------
  # TSS enrichment score is a raio between aggregate distribution of reads centered on TSSs and that flanking the corresponding TSSs. TSS score = the depth of TSS (each 100bp window within 1000 bp each side) / the depth of end flanks (100bp each end). TSSE score = max(mean(TSS score in each window)). TSS enrichment score is calculated according to the definition at https://www.encodeproject.org/data-standards/terms/#enrichment. Transcription start site (TSS) enrichment values are dependent on the reference files used; cutoff values for high quality data are listed in the following table from https://www.encodeproject.org/atac-seq/.
  
  tsse <- TSSEscore(gal1, txs)
  tsse$TSSEscore
  plot(100*(-9:10-.5), tsse$values, type="b", 
       main=paste(sampleID, ", TSSEscore:", tsse$TSSEscore%>%round(2)),
       xlab="distance to TSS",
       ylab="aggregate TSS score")
  
  
  
  #txs <- txs[seqnames(txs) %in% "Gm01"]
  
  ## split reads by fragment length
  objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)
  ## list the files generated by splitGAlignmentsByCut.
  dir(outPath)
  
  
  bamfiles <- file.path(outPath,
                        c("NucleosomeFree.bam",
                          "mononucleosome.bam",
                          "dinucleosome.bam",
                          "trinucleosome.bam"))
  
  ## Plot the cumulative percentage of tag allocation in nucleosome-free 
  ## and mononucleosome bam files.
  #cumulativePercentage(bamfiles[1:2], as(seqinformation["Gm01"], "GRanges"))
  cumulativePercentage(bamfiles[1:2], as(seqinformation, "GRanges"))
  
  TSS <- promoters(txs, upstream=0, downstream=1)
  TSS <- unique(TSS)
  ## estimate the library size for normalization
  (librarySize <- estLibSize(bamfiles))
  
  ## calculate the signals around TSSs.
  NTILE <- 101
  dws <- ups <- 1010
  sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                       "mononucleosome",
                                       "dinucleosome",
                                       "trinucleosome")], 
                            TSS=TSS,
                            librarySize=librarySize,
                            seqlev=seqlev,
                            TSS.filter=0.5,
                            n.tile = NTILE,
                            upstream = ups,
                            downstream = dws)
  ## log2 transformed signals
  sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
  #plot heatmap
  featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                        zeroAt=.5, n.tile=NTILE)
  
  ## get signals normalized for nucleosome-free and nucleosome-bound regions.
  out <- featureAlignedDistribution(sigs, 
                                    reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt=.5, n.tile=NTILE, type="l", 
                                    ylab="Averaged coverage")
  
  ## rescale the nucleosome-free and nucleosome signals to 0~1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  out <- apply(out, 2, range01)
  matplot(out, type="l", xaxt="n", 
          xlab="Position (bp)", 
          ylab="Fraction of signal")
  axis(1, at=seq(0, 100, by=10)+1, 
       labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
  abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
  
  dev.off()

