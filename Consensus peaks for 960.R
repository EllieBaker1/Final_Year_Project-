library(DiffBind)
#Create SampleList#
Samples <- dir("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/", pattern="_peak2.narrowPeak")
SampleList <- data.frame(SampleID=substr(Samples,1,8),
                         Tissue=substr(Samples,5,6),
                         Factor="",
                         Condition=substr(Samples,5,6),
                         Treatment="5hmC",
                         Replicate="",
                         bamReads=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/", substr(Samples,1,12), "-chr19.bam"),
                         ControlID="",
                         bamControl="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Merged_Inp.bam",
                         Peaks=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/", substr(Samples,1,12) , "_filtered_peak2.narrowPeak"),
                         PeakCaller="narrow")
SampleList$Tissue[SampleList$Tissue %in% c("LT", "NL")] <- "Liver"
SampleList$Tissue[SampleList$Condition %in% c("NC", "TC")] <- "Colon"
SampleList$Condition[SampleList$Condition %in% c("LT", "TC")] <- "Tumour"
SampleList$Condition[SampleList$Condition %in% c("NL", "NC")] <- "Normal"
SampleList$Replicate <- c(1:4, 1:3, 1:3, 1:3)
write.csv(SampleList, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Filtered_peaks/SampleList_filteredpeaks4.csv",quote=F, row.names=F)

#Read SampleList into R#
SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Filtered_peaks/SampleList_filteredpeaks4.csv", header=T)#loaded in the samples filtered against the bed file for the 960 peaks#

#Create dba (an object in which data can be stored)#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Filtered_peaks/SampleList_filteredpeaks4.csv")
dba.plotHeatmap(myDBA)

NLLTNCTC <- dba(myDBA) #select for all samples in csv#
NLLTNCTC


####Consensus peak and occupancy analysis####
myDBA_consensus_NLLTNCTC <- dba.peakset(NLLTNCTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_NLLTNCTC <- dba(myDBA_consensus_NLLTNCTC, mask=myDBA_consensus_NLLTNCTC$masks$Consensus, minOverlap=1)
myDBA_consensus_NLLTNCTC #returns number of consensus peaks for all samples across all the 960 similiar peaks# 
consensus_peaks_NLLTNCTC <- dba.peakset(myDBA_consensus_NLLTNCTC, bRetrieve=TRUE) 
consensus_peaks_NLLTNCTC
myDBA_consensus_NLLTNCTC
dba.plotVenn(myDBA_consensus_NLLTNCTC, myDBA_consensus_NLLTNCTC$masks$Consensus)
NLLTNCTC.count <- dba.count(NLLTNCTC, consensus_peaks_NLLTNCTC, summits = 250) #counts number of consensus peaks for each sample. Summits makes the width of all the peaks the same (250 bp up and downstream of peak summit)
NLLTNCTC.count #displays number of consensus peaks across samples#
NLLTNCTC_normed <- dba.normalize(NLLTNCTC.count) #normalises number of peaks based on lib method - normalises by the mean number of reads across all samples being analysed#
NLLTNCTC_normed
hmap <- colorRampPalette(c("white", "grey", "black"))(n = 13)
dba.plotHeatmap(NLLTNCTC_normed, correlations = TRUE, colScheme = "Greens")#plot heatmaps of normalised 960 peak data#
dba.plotHeatmap(NLLTNCTC_normed, correlations = FALSE, colScheme = "hmap")#plot heatmaps of normalised 960 peak data#
