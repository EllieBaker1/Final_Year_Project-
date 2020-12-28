#####Load packages#####
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T) #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv") #create dba (object in which data can be stored)#
dba.plotHeatmap(myDBA, correlations = TRUE)#plot heatmaps #
hmap <- colorRampPalette(c("white", "grey", "black"))(n = 13)
dba.plotHeatmap(myDBA, correlations = FALSE, colScheme = hmap )#plot heatmaps #

LTTC <- dba(myDBA, mask=myDBA$masks$Tumour) #select for only tumour samples#

LTTC #number of enriched peaks called for the tumour samples by MACS2#

####Consensus peak and occupancy analysis####
myDBA_consensus_LTTC <- dba.peakset(LTTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_LTTC <- dba(myDBA_consensus_LTTC, mask=myDBA_consensus_LTTC$masks$Consensus, minOverlap=1)
myDBA_consensus_LTTC #number of consensus peaks for Liver Tumour across all LT replicates and Tumour Colon across all TC replicates# 
consensus_peaks_LTTC <- dba.peakset(myDBA_consensus_LTTC, bRetrieve=TRUE) 
consensus_peaks_LTTC
dba.plotVenn(myDBA_consensus_LTTC, myDBA_consensus_LTTC$masks$Consensus)#plot venn diagram#


####Differential binding analysis####
LTTC.count <- dba.count(LTTC, consensus_peaks_LTTC, summits = 250) #counts number of consensus peaks for each sample. Summits makes the width of all the peaks the same - 250 bp up around peak summit
LTTC.count #number of consensus peaks across samples#
LTTC_normed <- dba.normalize(LTTC.count) #normalises number of peaks based on lib method#
LTTC_normed
dba.plotPCA(LTTC, label=1)
LTTC_contrast <- dba.contrast(LTTC_normed, categories=DBA_CONDITION) #Tells diffbind which comparison we are interested in. here we want to compare the 2 conditions for Tumour: Liver and tumour#
LTTC_contrast #displays which group each sample is in for the differential analysis#
LTTC_analyzed <- dba.analyze(LTTC_contrast, bBlacklist = FALSE, bGreylist = FALSE) #Runs default differential binding analysis - DESeq2. Our data has already been filtered against a "greylist#
LTTC.peaks <- dba.report(LTTC_analyzed, bCalled = TRUE, th=1) #Shows the differential binding analysis results and which peaks were called in each group. th=1 returns all peak results#
data("TSS.human.GRCh38") #Load human genome dataset from org.Hs.eg.db package#
LTTC.peaks <- annotatePeakInBatch(LTTC.peaks, AnnotationData=TSS.human.GRCh38) # match the peaks found with genes in the human genome dataset TSS#
LTTC.peaks <- addGeneIDs(LTTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol')) #annotate the peaks found with gene names#
LTTC.peaks

library(tidyverse) 
out <- as.data.frame(LTTC.peaks)
write.csv(out, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_NLLT/Normal.Liver_vs_Liver.Tumour_desq21.csv")#write data to file#

Tumour.Colon_enrich <- out %>%
  filter(FDR < 0.1 & Fold > 0) %>% #allows to create csv file keeping only the significant peaks for tumour colon (p < 0.05). Fold at > 0 pulls out all values that were called in group 1 at the contrast stage)
  select(seqnames, start, end, feature, FDR) #variables included in file# 
write.csv(Tumour.Colon_enrich, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_LTTC/Tumour.Colon_enriched_LTTC.csv") #write to file#

Liver.Tumour_enrich<- out %>%
  filter(FDR < 0.1 & Fold < 0) %>% #allows to create csv file keeping only the significant peaks forLiver Tumour (p < 0.05). Fold at < 0 pulls out all values that were called in group 2 at the contrast stage
  select(seqnames, start, end, feature, FDR) #variables included in file#
write.csv(Liver.Tumour_enrich, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_LTTC/Liver.Tumour_enriched_LTTC.csv")#write to file#

Similarties.LTTC <-out %>%
  filter(FDR > 0.1) %>%
  select(seqnames, start, end, feature, FDR) 
write.csv(Similarties.LTTC, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_LTTC/Similarities.LTTC.csv") #write to file# 
