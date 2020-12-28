#####Load packages#####
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T) #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv") #create dba (object in which data can be stored)#
NLLTNCTC <- dba(myDBA) #select for all samples#
NLLTNCTC


####Consensus peak and occupancy analysis####
myDBA_consensus_NLLTNCTC <- dba.peakset(NLLTNCTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_NLLTNCTC <- dba(myDBA_consensus_NLLTNCTC, mask=myDBA_consensus_NLLTNCTC$masks$Consensus, minOverlap=1)
myDBA_consensus_NLLTNCTC #returns number of consensus peaks for all samples# 
consensus_peaks_NLLTNCTC <- dba.peakset(myDBA_consensus_NLLTNCTC, bRetrieve=TRUE) 
consensus_peaks_NLLTNCTC
myDBA_consensus_NLLTNCTC
dba.plotVenn(myDBA_consensus_NLLTNCTC, myDBA_consensus_NLLTNCTC$masks$Consensus)
NLLTNCTC.count <- dba.count(NLLTNCTC, consensus_peaks_NLLTNCTC, summits = 250) #counts number of consensus peaks for each sample. Summits makes the width of all the peaks the same (250 bp up and downstream of peak summit)
NLLTNCTC.count #displays number of consensus peaks across samples#
NLLTNCTC_normed <- dba.normalize(NLLTNCTC.count) #normalises number of peaks based on lib method - normalises by the mean number of reads across all samples being analysed#
NLLTNCTC_normed
dba.plotPCA(NLLTNCTC_normed, label=1)

dba.plotHeatmap(NLLTNCTC_normed, correlations = TRUE, colScheme = "custom")#plot heatmaps of normalised data#
dba.plotHeatmap(NLLTNCTC_normed,  contrast = 1, correlations = FALSE, colScheme = "Greens")#plot heatmaps of normalised data#

NLLTNCTC_contrast <- dba.contrast(NLLTNCTC_normed, categories=DBA_CONDITION) #Tells diffbind which comparison we are interested in - here it is all conditions#
NLLTNCTC_contrast #displays which group each sample is in for the differential analysis#
NLLTNCTC_analyzed <- dba.analyze(NLLTNCTC_contrast, bBlacklist = FALSE, bGreylist = FALSE) #Runs default differential binding analysis -DESeq2. Data has already been filtered against a "greylist#
NLLTNCTC_analyzed
NLLTNCTC.peaks <- dba.report(NLLTNCTC_analyzed, bCalled = TRUE, th=1) #Shows the differential binding analysis results and peaks called in group1 and group2. th=1 returns all peak results #
NLLTNCTC.peaks
data("TSS.human.GRCh38") #Load human genome dataset from org.Hs.eg.db package#
NLLTNCTC.peaks <- annotatePeakInBatch(NLLTNCTC.peaks, AnnotationData=TSS.human.GRCh38) # match the peaks found with genes in the human genome dataset TSS#
NLLTNCTC.peaks <- addGeneIDs(NLLTNCTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol')) #annotate the peaks found with gene names#
NLLTNCTC.peaks
