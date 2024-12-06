# peakAnnt.Gm724.v.1.0.R
# author: YCW
# date: 20231222
# Description: run annotate to peak, this is for ATAC data set
# Version: GenomicFeatures_1.46.5; ChIPseeker_1.30.3
# Input:
#    1. path of peak files, accesapted format: bed, bedGraph, bedgraph(HMMRATAC outout) and broadPeak/gappedPeak
# Output:
#    1. annotation of each peak file
#    2. peak annotation comparision across all peak files 
# Reference: https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html

####################### install packages ######################
if (!require("GenomicFeatures", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GenomicFeatures")
}

if (!require("ChIPseeker", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ChIPseeker")
}

if (!require("DOSE", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("DOSE")
}

if (!require("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}

if (!require("ggupset", quietly = TRUE)) {
  install.packages("ggupset")
}

if (!require("forcats", quietly = TRUE)) {
  install.packages("forcats")
}

if (!require("Vennerable", quietly = TRUE)) {
  install.packages("Vennerable", repos="http://R-Forge.R-project.org", dependencies = T)
}
################################################


rm(list=ls())
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(ChIPseeker)
library(DOSE)
library(ggplot2)
library(Vennerable)
#library(ggvenn)
#library(VennDiagram)
library(ggpubr)
library(forcats)
library(tibble)

sessionInfo()

#----------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
cat("
          #USAGE:Rscript --vanilla /bcst/JYL/YCWang/Rscript/peakAnnt.Gm724.v.1.0.R [arg1] [arg2]\n
          #Arg1: input folder path, ex. /bcst/JYL/YCWang/testing/ATAC_test/peakannt/bed \n
          #Arg2: output file path ex. /bcst/JYL/YCWang/testing/ATAC_test/peakannt \n
          #Arg3: type of input files ex. bed \n
          ")

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("3 arguments must be supplied (input).n", call.=FALSE)
}

#----------------------------------------------------------------------------
datapath<- args[1]
outpath<- args[2]
filestype<- args[3]



####################### prepare Txdb ######################
# gffread /bcst/JYL/JYL_qnap/db/Gm/GmaxWm82ISU_01/v2.1/annotation/GmaxWm82ISU_01_724_v2.1.gene_exons.gff3 -T -o /bcst/JYL/JYL_qnap/db/Gm/GmaxWm82ISU_01/v2.1/annotation/GmaxWm82ISU_01_724_v2.1.gene_exons.gtf
# makeTxDbFromGFF(file,
#                 format=c("auto", "gff3", "gtf"),
#                 dataSource=NA,
#                 organism=NA,
#                 taxonomyId=NA,
#                 circ_seqs=NULL,
#                 chrominfo=NULL,
#                 miRBaseBuild=NA,
#                 metadata=NULL,
#                 dbxrefTag)

# seqlength<-read.delim("/bcst/JYL/JYL_qnap/db/Gm/GmaxWm82ISU_01/v2.1/assembly/Chrom_info_Wm82ISU01gnm2ann1", stringsAsFactors = F)
# txdb <- makeTxDbFromGFF(file="/bcst/JYL/JYL_qnap/db/Gm/GmaxWm82ISU_01/v2.1/annotation/GmaxWm82ISU_01_724_v2.1.gene_exons.gff3",
#                         format = "gff3", 
#                         dataSource="http://www.phytozome.org/",
#                         organism="Glycine max",
#                         chrominfo =Seqinfo(seqnames=seqlength$Chro,
#                                            seqlengths=seqlength$Length_bp,
#                                            genome="Gm724"))


#saveDb(txdb, "/bcst/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/TxDb/Gm724")

txdb <- loadDb("/bcst/JYL/JYL_qnap_2/YCWang/reference/Gm_ref/TxDb/Gm724")
#gene<-genes(txdb)%>%as.data.frame()


## run ----------------------------------------------------------------------------
setwd(outpath)

# input peak files -------
if(filestype %in%  c("bed")) {
  bedfiles = datapath %>% list.dirs() %>% list.files(".bed", full.names = T)
  cat("input file type: bed files")
  
} else if (filestype %in%  c("broadPeak")) {
  bedfiles = datapath %>% list.dirs() %>% list.files(".broadPeak$", full.names = T)
  cat("input file type: broadPeak files")
  
} else if (filestype %in%  c("gappedPeak")) {
  bedfiles = datapath %>% list.dirs() %>% list.files(".gappedPeak$", full.names = T)
  cat("input file type: gappedPeak files")
  

} else if (filestype %in%  c("bedGraph", "bdg")) {
  bedfiles = datapath %>% list.dirs() %>% list.files("bedGraph$", full.names = T)
  cat("input file type: bedGraph files")
  
} else if (filestype %in%  c("bedgraph", "bdg")) {
  bedfiles = datapath %>% list.dirs() %>% list.files("bedgraph$", full.names = T)
  cat("input file type: bedgraph files")  
  
}else{
  cat("input file type: can't be detected, please check. 
        accesapted format: bed, bedGraph, bedgraph(HMMRATAC outout) and broadPeak/gappedPeak")
}
cat("\n")
cat("total files read ")
cat(bedfiles)
cat("\n")

## run annotation for individul peak file ----------------------------------------------------------------------------
Bedfiles=NULL
for ( bed in 1:length(bedfiles)) {
sampleID<- bedfiles[bed]%>%as.character()%>%basename()%>%str_remove(".\\w+$")
cat(sampleID)

if(!dir.exists(file.path(sampleID))){
  dir.create(file.path(sampleID), mode = "0770", showWarnings = F)
}


peak <- readPeakFile(bedfiles[bed], as="GRanges", header=FALSE)
if(length(seqlevels(peak)%>%str_detect("Gm01")%>%which())==0){
  for (i in 1:9) {
    seqlevels(peak) <-
      sub("chr", "Gm", seqlevels(peak))
    seqlevels(peak) <-
      sub(paste0("Gm", i), paste0("Gm0", i), seqlevels(peak))
  }
}
 
head(peak)
print(length(peak))
# Assume peak is defined somewhere in your code

# Initialize peakaAnno as "failure"
peakanno <- "failure"

# Use tryCatch to catch errors
tryCatch({
  #
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), level = "gene",
                                        TxDb = txdb, verbose = F )
  print(peakAnno@annoStat)
  peakanno <- "done"
}, error = function(e) {
  # Code to handle the error
  print("An error occurred:")
  print(e)
})
print(peakanno)
# Check if peakaAnno is a failure
if (peakanno == "failure") {
  print("Annotation failed.")
} else if (peakanno == "done") {
  print("Annotation successful.")
  pdf(file.path(sampleID, paste0(sampleID, "ChIPseeker.plot.pdf")))
      
  
     
      # Visualize Genomic Annotation
      print("plotAnno plotting")
      print(plotAnnoPie(peakAnno))
      print(plotAnnoBar(peakAnno))
      #vennpie(peakAnno)
      
      print("upsetplot plotting")
      print(upsetplot(peakAnno))
      #covplot(peak, weightCol="V5", chrs=c("Gm01", "Gm02"))
      
      print("covplot plotting")
      if(filestype == "broadPeak"){
        print(covplot(peak, weightCol = "V5"))
      }else if(filestype == "bedGraph"){
        print(covplot(peak, weightCol = "V4"))
      }else if(filestype %in% c("bed", "bedgraph")){
        print(covplot(peak))
      }
      
      #peak=GenomicRanges::GRangesList(CBX6=readPeakFile(files[[4]]),CBX7=readPeakFile(files[[5]]))
      
      
      # Visualize distribution of TF-binding loci relative to TSS
      print("plotDistToTSS plotting")
      print(plotDistToTSS(peakAnno,
                    title="Distribution of transcription factor-binding loci\nrelative to TSS"))
      
      print("plot done")
  
  dev.off()
  
  # Output annotation
  dfGRanges = data.frame(as.GRanges(peakAnno))
  head(dfGRanges)
  write.table(dfGRanges,file = file.path(sampleID, paste0(sampleID, ".ACR.PeakAnnotation.txt")),
              sep="\t",row.names=F, quote = F)
  
  Bedfiles<-Bedfiles%>%bind_rows(data.frame(bed=bedfiles[bed]))
}
}



## run peak annotation comparision ----------------------------------------------------------------------------
# Read BED files into a GRangesList
grl <- GRangesList()
print(Bedfiles$bed)
for (bedfile in Bedfiles$bed) {
  # gr <- read.table(bedfile, header = FALSE, colClasses = c("character", "integer", "integer", "integer"), col.names = c("seqnames", "start", "end", "value"))
  gr <- read.table(bedfile, header = FALSE)
  colnames(gr)[1:3]<-c("seqnames", "start", "end")
  gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE)
  grl[[basename(bedfile)]] <- gr
}

head(grl)
length(grl)
grl1<-grl

if(length(seqlevels(grl1)%>%str_detect("Gm01")%>%which())==0){
  for (i in 1:9) {
    seqlevels(grl1) <-
      sub("chr", "Gm", seqlevels(grl1))
    seqlevels(grl1) <-
      sub(paste0("Gm", i), paste0("Gm0", i), seqlevels(grl1))
  }
}
seqlevels(grl1)%>%str_detect("Gm01")%>%which()
head(grl1)

# Profile of several ATAC peak data binding to TSS region
# Average profiles
groupname<-datapath%>%as.character()%>%basename()

names(grl1) <- names(grl1)%>%str_remove(".\\w+$")
pdf(file.path(paste0(groupname, ".ChIPseeker.plot.pdf")))
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(grl1, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), 
            xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency", 
            resample=500,
            facet="row")

# ATAC peak annotation comparision
peakAnnoList <- lapply(grl1, annotatePeak, TxDb=txdb, level = "gene",
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
p<-plotDistToTSS(peakAnnoList)

if(length(tagMatrixList)<=6){
  heatmap.rank=names(tagMatrixList)
  print(heatmap.rank)
}else{
  summed_list <- tagMatrixList[order(sapply(tagMatrixList, sum), decreasing = TRUE)][1:6]
  heatmap.rank<-names(summed_list)
  print(heatmap.rank)
}

# Use tryCatch to catch errors
tryCatch({
  #
  tagMatrixList_rank <- tagMatrixList[heatmap.rank]
  #print(tagMatrixList_rank)
  tagHeatmap(tagMatrixList_rank, xlim=c(-3000, 3000), color=rainbow(length(tagMatrixList_rank))
  )
}, error = function(e) {
  # Code to handle the error
  print("An error occurred:")
  print(e)
})


tryCatch({
  # 
  if(filestype == "broadPeak"){
    #covplot(grl1, weightCol = "V5") + facet_grid(.id ~ chr)
    covplot(grl1, weightCol = "V5") + facet_grid(chr ~ fct_rev(.id))
  }else if(filestype == "bedGraph"){
    covplot(grl1, weightCol = "V4") + facet_grid(chr ~ fct_rev(.id))
  }
  
}, error = function(e) {
  # Code to handle the error
  print("An error occurred:")
  print(e)
})



# ATAC gene annotation intersection ------------

genes <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)$geneId) 
# unique
genes <- lapply(genes, unique)

# Count the number of variables in the list
num_variables <- sapply(genes, function(x) length(unlist(x)))
cat(num_variables)


if(length(genes)>4){
  tophalf<- (num_variables%>%sort(decreasing = T))[length(genes)/2]
}else{
  tophalf<- (num_variables%>%sort(decreasing = T))[length(genes)]
}


# pick top4 of gene number
filtered_half_list <- Filter(function(x) length(unlist(x)) >= tophalf, genes)
num_variables_filtered <- sapply(filtered_half_list, function(x) length(unlist(x)))
cat(num_variables_filtered)

# Function to calculate Jaccard similarity between two vectors
jaccard_similarity <- function(vec1, vec2) {
  intersection_size <- length(intersect(vec1, vec2))
  union_size <- length(union(vec1, vec2))
  return((intersection_size / union_size)%>% round(4))
}

# Calculate Jaccard similarities for all pairs
input_list<-filtered_half_list
pairwise_sims <- combn(names(input_list), 2,
                       function(pair_names) {
                         vec1 <- input_list[[pair_names[1]]]
                         vec2 <- input_list[[pair_names[2]]]
                         similarity <- jaccard_similarity(vec1, vec2)
                         return(c(pair_names, similarity))
                       }, simplify = TRUE)

# Create a data frame with results
similarity_df <- data.frame(Vector1 = pairwise_sims[1,],
                            Vector2 = pairwise_sims[2,],
                            Similarity = pairwise_sims[3,]) %>%
  mutate(Similarity=Similarity%>%as.numeric())

write.table(similarity_df,file = file.path(paste0(groupname, "similarity.ACR.geneAnnotation.txt")),
            sep="\t",row.names=F, quote = F)

# Identify the four most similar pairs
top_4<-NULL
#start_from=8

find_top4 <- function(similarity_df, start_from) {
  cat(paste( "start_from searching number:", start_from))
  repeat {
    top_4<-NULL
    top_4_similar <- similarity_df[order(similarity_df$Similarity, decreasing = T), ][1:start_from, ]
    # cat the results
    # cat(similarity_df)
    cat("\n")
    cat("Top 4 Similar Pairs searching:")
    print(top_4_similar)
    top_4<- c(top_4_similar$Vector1, top_4_similar$Vector2)%>%unique()
    #cat(top_4)
    # check for success
    if (length(top_4) <= 4) break
    start_from=start_from-1
    cat(paste( "next searching number:", start_from))
    cat("\n")
  }
  cat("\n")
  cat("final return ------------------\n")
  cat(paste("when top pair-wise number is", start_from, "\n"))
  cat(paste("Top 4 Similar Pairs found. \n"))
  return(top_4)
  
}


start_from<-length(bedfiles)/2
if(length(genes)>4){
  top4<- find_top4(similarity_df, start_from)
}else{
  top4<- c(similarity_df$Vector1, similarity_df$Vector2)%>%unique()
}
print(top4)




# pick top4 of overlap gene number
filtered_fin_list <- genes[top4]
num_variables_filtered <- sapply(filtered_fin_list, function(x) length(unlist(x)))
cat(num_variables_filtered)

# venn plot
#pdf()
#vennplot(filtered_fin_list, by='Vennerable')
#dev.off()

#ggvenn(filtered_fin_list, set_name_size = 8, text_size = 6, show_outside= "auto") 
pcentFun <- function(x) {
  100 * (x /sum(x))
}
V <- Venn(filtered_fin_list)
if(length(filtered_fin_list)<=3){
  VennList <- compute.Venn(V, doWeights = F, type = "circles")
}else{
  VennList <- compute.Venn(V, doWeights = F, type = "ellipses")
}

Weight<- V@IndicatorWeight %>%as.data.frame()%>%rownames_to_column 
Weight_sorted <- Weight[order(match(Weight$rowname, VennList@FaceLabels$FaceName)), ]
areas<-Weight_sorted$.Weight
names(areas)<-Weight_sorted$rowname
areasPcent <- round(pcentFun(areas), digits=2)
VennList@FaceLabels$Signature <- paste0(areas, "\n", areasPcent%>%round(1), "%")

#plot(VennList, show = list(FaceText = c("signature"), DarkMatter = F))
gp <- VennThemes(VennList,increasingLineWidth=F, colourAlgorithm = "signature")
modify_vector <- function(element, option, value) {
  if (option %in% names(element)) {
    element[[option]] <- value
  }
  return(element)
}

gp[["SetText"]] <- lapply(gp[["SetText"]], function(x) modify_vector(x, "fontsize", 12))
gp[["FaceText"]] <- lapply(gp[["FaceText"]], function(x) modify_vector(x, "fontsize", 12))

gridExtra::grid.arrange(grid::grid.grabExpr(height=4, width=6, plot(
  VennList, show = list(FaceText = c("signature"),
                        DarkMatter = F),
  gp=gp
)),
top = textGrob(
  "Intersection of Gene Sets from Top 4 Similar Samples in Venn",
  gp = gpar( fontsize = 15)))


dev.off()


cat("job finish")



