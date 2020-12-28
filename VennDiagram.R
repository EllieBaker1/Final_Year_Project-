#####Load packages#####
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T) #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv") #create dba (an object within which data can be stored)#
NLLTNCTC <- dba(myDBA) #select for only liver samples#
NLLTNCTC


####Consensus peak and occupancy analysis####
myDBA_consensus_NLLTNCTC <- dba.peakset(NLLTNCTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_NLLTNCTC <- dba(myDBA_consensus_NLLTNCTC, mask=myDBA_consensus_NLLTNCTC$masks$Consensus, minOverlap=1)
myDBA_consensus_NLLTNCTC #returns number of consensus peaks for Liver Tumour across all LT replicates and Normal Liver across all NL replicates# 
consensus_peaks_NLLTNCTC <- dba.peakset(myDBA_consensus_NLLTNCTC, bRetrieve=TRUE) 
consensus_peaks_NLLTNCTC
myDBA_consensus_NLLTNCTC
dba.plotVenn(myDBA_consensus_NLLTNCTC, myDBA_consensus_NLLTNCTC$masks$Consensus)


vennPeaks <- dba.plotVenn(myDBA_consensus_NLLTNCTC, myDBA_consensus_NLLTNCTC$masks$Consensus)
dba.plotVenn
peaksInAll = vennPeaks$inAll
bed <- data.frame(seqnames=seqnames(peaksInAll),
                  starts=start(peaksInAll) - 1,
                  ends=end(peaksInAll))
write.table(bed, file="E:/blacklist_file/peaksInAll.bed", quote=F, sep="\t", row.names=F, col.names=F)#produces a bed file that can be used to filter against on the command line#

