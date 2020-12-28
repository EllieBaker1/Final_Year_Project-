##Load Package##
library("DiffBind")

#Create SampleList#
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

#Read SampleList into R#
SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T)

#Create dba (object in which data can be stored)#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv")
plot(myDBA)

#Select for liver samples#
NLLT <- dba(myDBA, mask=myDBA$masks$Liver)
NLLT #Returns number of enriched peaks called for in MACS2#


##Consensus peak and occupancy analysis##
myDBA_consensus <- dba.peakset(NLLT, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks#
myDBA_consensus <- dba(myDBA_consensus, mask=myDBA_consensus$masks$Consensus, minOverlap=1)
myDBA_consensus
consensus_peaks <- dba.peakset(myDBA_consensus, bRetrieve=TRUE) 
consensus_peaks
#Plot Venn diagram comparing peaks in each sample#
dba.plotVenn(myDBA_consensus, myDBA_consensus$masks$Consensus)

##Differential Binding Analysis##
NLLT <- dba.count(NLLT, summits=250, peaks=consensus_peaks) #counts number of consensus peaks (peaks that are found in at least two of the sample replicates)#
NLLT
NLLT <- dba.normalize(NLLT)#Data are normalized based on sequencing depth#
NLLT
dba.plotHeatmap(NLLT_normed, correlations = TRUE, colScheme = "Greens")#plot heatmaps of normalised data#
dba.plotHeatmap(NLLT_normed, correlations = FALSE, colScheme = "Greens")#plot heatmaps of normalised data#
NLLT <- dba.contrast(NLLT, categories=DBA_CONDITION)#tells DiffBind which comparisons we are interested in; here we are comparing Normal Liver and Liver tumour#
dba.contrast(NLLT)#confirms which group each sample is in#
NLLT <- dba.analyze(NLLT, bBlacklist=FALSE, bGreylist=FALSE)#This will run the default DESeq2 analysis on samples to produce an object that states which sites are differentially bound#
NLLT.peaks <- dba.report(NLLT, bCalled=TRUE, th=1)#This retrieves the differentially bound sites as GRanges object#
NLLT.peaks


#load library#
library(org.Hs.eg.db)
library(ChIPpeakAnno)

data(TSS.human.GRCh38)#Load human genome dataset from package stated#
NLLT.peaks <- annotatePeakInBatch(NLLT.peaks, AnnotationData=TSS.human.GRCh38)#match peaks from our dataset with genes in human genome dataset#
NLLT.peaks <- addGeneIDs(NLLT.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol')) #annotate peaks with gene names#
NLLT.peaks


#old method#

write.csv(NLLT.peaks, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/NLLT.peaks.csv")#Differential binding results (call1 and call2 and gene names)#

UniqueNL <- NLLT.peaks[which(NLLT.peaks$Called1 >=2 & NLLT.peaks$Called2 < 2 & NLLT.peaks$FDR <= 0.1)]#provides the peaks called in at least two of the replicates in group 1 but less than two in group 2 - unique to peak1#
UniqueLT <- NLLT.peaks[which(NLLT.peaks$Called1 < 2 & NLLT.peaks$Called2 >= 2 & NLLT.peaks$FDR <= 0.1)]#display peaks called in at least 2 replicates for group2 but less than 2 for group2 (unique to group2)#
DiffPeak <- NLLT.peaks[which(NLLT.peaks$Called1 >=2 & NLLT.peaks$Called2 >= 2 & NLLT.peaks$FDR <= 0.1)]#displays common peaks between the two groups#

write.csv(UniqueNL, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueNL_NLLT.csv")
write.csv(UniqueLT, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueLT_NLLT.csv")
write.csv(DiffPeak, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/DiffPeak_NLLT.csv")

##Heatmaps of NLLT normalised#
dba.plotHeatmap(NLLT_normed, correlations = TRUE, colScheme = "Greens")
dba.plotHeatmap(NLLT_normed, correlations = FALSE, colScheme = "Greens")

