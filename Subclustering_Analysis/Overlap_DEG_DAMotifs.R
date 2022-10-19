#Title: Overlap lists of DEG and differentially active motifs to identify high confidence transcription factors

library(Seurat)
library(Signac)
library(dplyr)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)
library(GenomicRanges)
library(tibble)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

# Overlap lists of upregulated genes and motifs with increased activity 
DEG <- read.csv("./Subclustering_Analysis/Supplemental_Table_4.csv",header=T)
DEG <- DEG[which(DEG$avg_log2FC>=0.25),]
colnames(DEG) <- c(paste0(rep("GEX.",5),colnames(DEG)[1:5]),"cluster","gene")
DEG$gene <- toupper(DEG$gene)

Motif <- read.csv("./Subclustering_Analysis/PTCSUB_chromvar_motif_logFC0_minpct0.1_padj0.01.csv",header=T)
colnames(Motif)
Motif <- Motif[,-1]
Motif$gene <- toupper(Motif$gene)
Motif <- Motif[which(Motif$motif.avg_diff.>0.1),]
OVERLAP <- data.frame("cluster"=c(),"gene"=c(),"motif"=c(),"motif.avg_diff."=c(), "motif.p_val"=c(),"motif.p_val_adj"=c(),"motif.pct.1"=c(),"motif.pct.2"=c(), 
                      "GEX.avg_log2FC"=c(),"GEX.p_val"=c(),"GEX.p_val_adj"=c(), "GEX.pct.1"=c(),"GEX.pct.2"=c())

Celltypes <- unique(Motif$cluster)

for(i in 1:length(Celltypes)){
  DEG_celltype <- DEG[which(DEG$cluster == Celltypes[[i]]),]
  Motif_celltype <- Motif[which(Motif$cluster== Celltypes[[i]]),]
  Motif_overlap <- Motif_celltype[which(Motif_celltype$gene%in%DEG_celltype$gene),]
  DEG_overlap <- DEG_celltype[which(DEG_celltype$gene%in%Motif_overlap$gene),]
  myDF <- merge(Motif_overlap,DEG_overlap,by=c("gene","cluster"))
  OVERLAP <- rbind(OVERLAP,myDF)
  }

dim(OVERLAP) 

# Overlap lists of downregulated genes and motifs with decreased activity 
DEG <- read.csv("./Subclustering_Analysis/Supplemental_Table_4.csv",header=T)
DEG <- DEG[which(DEG$avg_log2FC<0),]
colnames(DEG) <- c(paste0(rep("GEX.",5),colnames(DEG)[1:5]),"cluster","gene")
DEG$gene <- toupper(DEG$gene)

Motif <- read.csv("./Subclustering_Analysis/PTCSUB_chromvar_motif_logFC0_minpct0.1_padj0.01.csv",header=T)
colnames(Motif)
Motif <- Motif[,-1]
Motif$gene <- toupper(Motif$gene)
Motif <- Motif[which(Motif$motif.avg_diff.< (-0.1)),]
#Motif <- Motif[which(Motif$motif.avg_diff.>0.1),]
OVERLAP_NEG <- data.frame("cluster"=c(),"gene"=c(),"motif"=c(),"motif.avg_diff."=c(), "motif.p_val"=c(),"motif.p_val_adj"=c(),"motif.pct.1"=c(),"motif.pct.2"=c(), 
                      "GEX.avg_log2FC"=c(),"GEX.p_val"=c(),"GEX.p_val_adj"=c(), "GEX.pct.1"=c(),"GEX.pct.2"=c())

Celltypes <- unique(Motif$cluster)

for(i in 1:length(Celltypes)){
  DEG_celltype <- DEG[which(DEG$cluster == Celltypes[[i]]),]
  Motif_celltype <- Motif[which(Motif$cluster== Celltypes[[i]]),]
  Motif_overlap <- Motif_celltype[which(Motif_celltype$gene%in%DEG_celltype$gene),]
  DEG_overlap <- DEG_celltype[which(DEG_celltype$gene%in%Motif_overlap$gene),]
  myDF <- merge(Motif_overlap,DEG_overlap,by=c("gene","cluster"))
  OVERLAP_NEG <- rbind(OVERLAP_NEG,myDF)
}

dim(OVERLAP_NEG) 

#MERGE both lists
OVERLAP <- rbind(OVERLAP, OVERLAP_NEG)
dim(OVERLAP) 

OVERLAP$Motif_GEX_score <- OVERLAP$motif.avg_diff.* OVERLAP$GEX.avg_log2FC

#Order DF
OVERLAP <- OVERLAP[order(-OVERLAP$Motif_GEX_score),]
OVERLAP <- OVERLAP %>% arrange(factor(cluster, levels = Celltypes))
write.csv(OVERLAP, "PTC_SUB_TFs_concordant_motif_GEX_changes.csv",row.names = F)
