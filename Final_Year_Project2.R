library("DiffBind")

dirMacs2 = "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/"
Samples <- dir("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/", pattern="_peaks.narrowPeak")
SampleList <- data.frame(SampleID=substr(Samples,1,8),
                         Tissue=substr(Samples,5,6),
                         Factor="",
                         Condition=substr(Samples,5,6),
                         Treatment="5hmC",
                         Replicate="",
                         bamReads=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/", substr(Samples,1,12), "-chr19.bam"),
                         ControlID="",
                         bamControl="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Merged_Inp.bam",
                         Peaks=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/", substr(Samples,1,12) , "-chr19_macs2_peaks.narrowPeak"),
                         PeakCaller="narrow")
SampleList$Tissue[SampleList$Tissue %in% c("LT", "NL")] <- "Liver"
SampleList$Tissue[SampleList$Condition %in% c("NC", "TC")] <- "Colon"
SampleList$Condition[SampleList$Condition %in% c("LT", "TC")] <- "Tumour"
SampleList$Condition[SampleList$Condition %in% c("NL", "NC")] <- "Normal"
SampleList$Replicate <- c(1:5, 1:3, 1:3, 1:3)
write.csv(SampleList, file="SampleList1.csv",quote=F, row.names=F)
SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T)

myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv")
plot(myDBA)
dba.plotHeatmap(myDBA, correlations = TRUE)
myDBA
names(myDBA$masks)

LTTC <- dba(myDBA, mask=myDBA$masks$Tumour)
LTTC
myDBA_consensus <- dba.peakset(LTTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2)
myDBA_consensus <- dba(myDBA_consensus, mask=myDBA_consensus$masks$Consensus, minOverlap=1)
consensus_peaks <- dba.peakset(myDBA_consensus, bRetrieve=TRUE)

LTTC <- dba.count(LTTC, peaks=consensus_peaks, summits=250, bParallel = 30)
dba.contrast(LTTC)
dba.plotHeatmap(LTTC, correlations = FALSE, colScheme = "Blues")
LTTC.OL <- dba.overlap(myDBA_consensus, myDBA_consensus$masks$Consensus)
dba.plotPCA(LTTC, label=1)
dba.plotVenn(myDBA_consensus, myDBA_consensus$masks$Consensus)
dba.plotHeatmap(LTTC, correlations = TRUE, colScheme = "Blues")

LTTC <- dba.contrast(LTTC, categories=DBA_CONDITION)
LTTC <- dba.analyze(LTTC)
LTTC.peaks <- dba.report(LTTC, bCalled=TRUE, th=1)
LTTC.peaks

library(org.Hs.eg.db)
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
LTTC.peaks <- annotatePeakInBatch(LTTC.peaks, AnnotationData=TSS.human.GRCh38)
LTTC.peaks <- addGeneIDs(LTTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol'))
LTTC.peaks

write.csv(LTTC.peaks, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/LTTC.peaks1.csv")
UniqueLT <- LTTC.peaks[which(LTTC.peaks$Called1 >=2 & LTTC.peaks$Called2 < 2 & LTTC.peaks$FDR <= 0.1)]
UniqueTC <- LTTC.peaks[which(LTTC.peaks$Called1 < 2 & LTTC.peaks$Called2 >= 2 & LTTC.peaks$FDR <= 0.1)]
DiffPeak <- LTTC.peaks[which(LTTC.peaks$Called1 >=2 & LTTC.peaks$Called2 >= 2 & LTTC.peaks$FDR <= 0.1)]
write.csv(UniqueLT, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueLT_LTTC1.csv")
write.csv(UniqueTC, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueTC_LTTC1.csv")
write.csv(DiffPeak, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/DiffPeak_LTTC1.csv")


