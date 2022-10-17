#Title: Plots for Figure 5
#Some parts from Parker Wilson's github: p4rkerw/Muto_Wilson_NComm_2020 (https://github.com/p4rkerw/Muto_Wilson_NComm_2020/blob/master/figures/figure3.R [last accessed 10/11/22])

library(Seurat)
library(Signac)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
set.seed(1234)

setwd()

#Read in the final R object of AKI_GFP+ ("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM") subclustering analysis 
PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")
PTC

#################################### Figure 5A - Heatmap DAR #################################### 
#Read in list of annotated differentially accessible regions
DAR <- read.csv("./Subclustering_Analysis/2022_5_31_peakAnno.1kb.allinfo.csv",header=T)

###Heatmap of DARs annotated as promoter (defined as 1000bp +/- TSS)
DAR_promoter <- DAR[which(DAR$annotation=="Promoter"),]

#Get average accessibility of DARs across clusters
dar_aver <- as.data.frame(AverageExpression(PTC, features = unique(DAR_promoter$DAR), assays = "peaks"))
dar_aver  <- dar_aver[do.call(order, c(dar_aver, list(decreasing=TRUE))),]
dar_aver$max <- max.col(dar_aver)
dar_aver  <- dar_aver[order(dar_aver$max), ]
dar_aver  <- dplyr::select(dar_aver, -max)
colnames(dar_aver) <- levels(Idents(PTC))
dar_aver <- t(dar_aver)

#Make heatmap
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
pdf("Figure_5A1.pdf",width=4.5,height=2)
pheatmap::pheatmap(dar_aver,scale = "column",
                   cluster_cols=F,cluster_rows = F,
                   color = col.pal,
                   show_rownames=T,show_colnames=F)
dev.off()

###Heatmap of DARs annotated as distal intergenic 
DAR_distal <- DAR[which(DAR$annotation=="Distal Intergenic"),]

#Get average accessibility of DARs across clusters
dar_aver <- as.data.frame(AverageExpression(PTC, features = unique(DAR_distal$DAR), assays = "peaks"))
dar_aver  <- dar_aver[do.call(order, c(dar_aver, list(decreasing=TRUE))),]
dar_aver$max <- max.col(dar_aver)
dar_aver  <- dar_aver[order(dar_aver$max), ]
dar_aver  <- dplyr::select(dar_aver, -max)
colnames(dar_aver) <- levels(Idents(PTC))
dar_aver <- t(dar_aver)

#Make heatmap
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
pdf("Figure_5A2.pdf",width=4.5,height=2)
pheatmap::pheatmap(dar_aver,scale = "column",
                   cluster_cols=F,cluster_rows = F,
                   color = col.pal,
                   show_rownames=T,show_colnames=F)
dev.off()

#################################### Figure 5B - Vcam1 Coverage plot #################################### 
#Read in matrix with all links
Links <- read.csv("./Subclustering_Analysis/PTC_SUB_allLinks.csv")

#Read in DAR
DAR <- read.csv("./Subclustering_Analysis/PTC_SUB_DAR_MACS2_logFC0.25_minpct0.05_padj0.01.csv")

#Get DAR in FR-PTC to highlight them in coverage plot
DAR_FRPTC <- DAR[which(DAR$cluster=="FR-PTC"),]
DAR_FRPTC.gr <- StringToGRanges(DAR_FRPTC$DAR, sep = c(":","-"))

#Get coordinates for DAR related to Vcam1 gene or links 
Vcam1_links <- Links[which(Links$gene =="Vcam1"),"peak"]
Vcam1_links <-StringToGRanges(Vcam1_links, sep = c(":","-"))

#Get overlap with DAR_FRPTC
overlap_peaks <- findOverlaps(DAR_FRPTC.gr,Vcam1_links)

#Get actual coordinates
Overlaps_Vcam1.gr <- granges(DAR_FRPTC.gr)[queryHits(overlap_peaks)]

#Make plot
pdf("Figure_5B.pdf", width=5, height=5)
CoveragePlot(
  object = PTC,
  region = "Vcam1",
  features = "Vcam1",
  expression.assay = "RNA",
  idents = levels(PTC),
  extend.upstream = 70000,
  extend.downstream = 80000,
  region.highlight = Overlaps_Vcam1.gr
)
dev.off()

#################################### Figure 5C & D - Heatmap motif activity & gene expression #################################### 

#Read in list with differentially expressed genes
DEG <- read.csv("./Subclustering_Analysis/PTCSUB_DEG_logFC0.25_minpct0.1_padj0.01.csv",header=T)
#Consider only up regulated genes
DEG <- DEG[which(DEG$avg_log2FC>=0.25),]
colnames(DEG) <- c(paste0(rep("GEX.",5),colnames(DEG)[1:5]),"cluster","gene")
DEG$gene <- toupper(DEG$gene)

#Read in list with differentially active motifs
Motif <- read.csv("./Subclustering_Analysis/PTCSUB_chromvar_motif_logFC0_minpct0.1_padj0.01.csv",header=T)
Motif$gene <- toupper(Motif$gene)

#Retain only differentially active motifs that are also differentially expressed
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

#Calculate a Motif_GEX_score
OVERLAP$Motif_GEX_score <- OVERLAP$motif.avg_diff.* OVERLAP$GEX.avg_log2FC

# Make heatmap with top 10 motifs per cluster (by Motif_GEX_score)
top10 <- OVERLAP  %>% group_by(cluster) %>% top_n(n = 10, wt = Motif_GEX_score)
top10$cluster <- factor(top10$cluster, levels= Celltypes)
top10 <- top10[order(top10$cluster),]
top10  <- top10  %>% group_by(cluster)
top10 <- top10 %>% arrange(desc(Motif_GEX_score), .by_group = TRUE)

#Get average expression per TF without exponentiation
PTC 
#The chromvar motif data stored in the data slot are z scores --> it should not be exponentiated (as is the default of AverageExpression())
#--> transfer the data into the counts slot and extract the average expression per cluster from there (no exponentiation will be performed when AverageExpression() is used on the counts slot)
PTC@assays$chromvar@counts <- PTC@assays$chromvar@data
aver_chromvar <- AverageExpression(PTC,assays = "chromvar", features = unique(top10$motif),slot="counts") %>% as.data.frame()
row.names(aver_chromvar) <- paste0(unique(top10$motif),sep="_",unique(top10$gene))
colnames(aver_chromvar) <- levels(Idents(PTC))

# Make heatmap
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

pdf("Figure_5C.pdf",width=5,height = 8)
pheatmap::pheatmap(aver_chromvar,scale = "row",
                   cluster_cols=F,cluster_rows = F,
                   color = col.pal,
                   show_rownames=T, angle_col=45)
dev.off()

# Make heatmap showing the gene expression of the top10 TFs across clusters
Genes_top10 <- paste0(substr(top10$gene,1,1),tolower(substr(top10$gene,2,10)))
aver <- AverageExpression(PTC,assays = "RNA", features = unique(Genes_top10),slot="data") %>% as.data.frame()
colnames(aver) <- levels(Idents(PTC))

pdf("Figure_5D.pdf",width=4,height = 8)
pheatmap::pheatmap(aver,scale = "row",
                   cluster_cols=F,cluster_rows = F,
                   color = col.pal,
                   show_rownames=T, angle_col=45)
dev.off()
