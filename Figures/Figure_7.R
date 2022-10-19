#Title: Plots for Figure 7

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)


setwd()

#Read in the final R object (contains RNA, chromVAR and peaks assay (peaks called using MACS2))
Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

#################################### Figure 7A - Umap plots #################################### 
# Umap plot by condition for control and AKI_GFP+
Idents(Ki67) <- "Group_time"
Control <- subset(Ki67, idents=c("Ctrl_4weeks","Ctrl_6months"))
AKI_GFPpos <- subset(Ki67, idents=c("AKI_GFP+_4weeks","AKI_GFP+_6months"))

pdf("./Figure_7A.pdf")
DimPlot(Control,reduction = "wnn.umap",label = F,cols=c("#b2b2b2","#333333")) + NoLegend()
levels(AKI_GFPpos) <- c("AKI_GFP+_4weeks","AKI_GFP+_6months")
DimPlot(AKI_GFPpos,reduction = "wnn.umap",label = F,cols=c("#8bd05a","#226400")) + NoLegend()
dev.off()

#################################### Figure 7B - Upset plot of DEG per cell type #################################### 
Celltype <- c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
              "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell")

#Remove Podocyte cluster (Comparison of AKI_GFP+ versus control cells not possible, because too few cells)
Celltype <- Celltype[-1]

#Make a list with all DEG AKI_GFP+ vs Control gene lists
DEG_list <- list()

for(i in 1:length(Celltype)){
  setwd("./GFPpos_vs_Ctrl/DEG")
  DEG <- read.csv(paste0("2022_3_31_DEG_GFPposvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  DEG_list[[Celltype[[i]]]] <- DEG[which(abs(DEG$avg_log2FC)>=0.25),"Gene"]
}

names(DEG_list) <-Celltype

#Make the combination matrix 
Mat <- make_comb_mat(DEG_list,mode = "distinct") 
Mat

#Show only sets with >= 10 genes in the plot
Mat1 <- Mat[comb_size(Mat) >= 10]

col_size = comb_size(Mat1)
row_size = set_size(Mat1)
ht = UpSet(Mat1, set_order = Celltype,
           top_annotation = upset_top_annotation(Mat1, ylim = c(0, max(col_size)*1.1)),
           right_annotation = upset_right_annotation(Mat1, ylim = c(0, max(row_size)*1.1)))

pdf("./Figure_7B.pdf")
ht <- draw(ht)
col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("intersection_size", {
  grid.text(col_size[col_od], x = seq_along(col_size), y = unit(col_size[col_od], "native") + unit(2, "pt"), default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

decorate_annotation("set_size", {
  grid.text(row_size[row_od], 
            unit(row_size[row_od], "native") + unit(2, "mm"), 
            rev(seq_len(length(row_size))), 
            default.units = "native", just = "bottom", rot = -90,
            gp = gpar(fontsize = 8))
})
dev.off()

#################################### Figure 7C - Violin plots of selected genes #################################### 
#subset only proximal tubule cells from Ki67 seurat object 
Idents(Ki67) <- "Celltype"

PT <- subset(Ki67,idents=c("PTS1_1","PTS1_2","PTS2","PTS1-S2","PTS3"))
PT

#Subset only AKI_GFP+ and control nuclei
Idents(PT) <- "Group"
PT <- subset(PT, idents=c("AKI_GFP+","Control"))

#Combine PTS1_1 and PTS1_2 for plotting
Idents(PT) <- "Celltype"
CELLTYPE <- as.vector(PT@meta.data$Celltype)
CELLTYPE[which(CELLTYPE%in%c("PTS1_2","PTS1_1"))] <- "PTS1"
PT@meta.data$Celltype_2 <- CELLTYPE

#Plot selected genes of interest
Genestoplot <- c("Nox4","Sash1", "Igfbp7","Slc34a1","Nrg1","Gclc")

pdf("Figure_7C.pdf",width=5,height=3)
Idents(PT) <- "Celltype_2"
levels(PT) <- c("PTS1","PTS1-S2","PTS2","PTS3")
for(i in 1:length(Genestoplot)){
  print(VlnPlot(PT,features = as.character(Genestoplot[[i]]),assay="RNA",split.by = "Group",pt.size = 0,cols = c("#128c26","#b8bdb5"))+ylab("Log normalized gene\n expression")+ggtitle(Genestoplot[[i]]))
}
dev.off()

#################################### Figure 7D - Volcano plot of DEG in PTS3: AKI_GFP+ versus Control #################################### 
#Read in list of differentially expressed genes between AKI_GFP+ PTS3 and Control PTS3
DEG <- read.csv("./GFPpos_vs_Ctrl/DEG/2022_3_31_DEG_GFPposvsCtrl_clusterPTS3_logFC0_minpct0.1_padj0.01.csv")

## Volcano plot
DEG <- DEG[which(DEG$p_val_adj<0.01),]
DEG  <- mutate(DEG , color = case_when(DEG$avg_log2FC > 0.25 & DEG$p_val_adj < 0.01 ~ "Increased expression in AKI_GFP+ PTS3",
                                       DEG$avg_log2FC < -0.25 & DEG$p_val_adj < 0.01 ~ "Increased expression in Control PTS3",
                                       abs(DEG$avg_log2FC) < 0.25 & DEG$p_val_adj < 0.01 ~ "nonsignificant",
                                       DEG$p_val_adj > 0.01 ~ "nonsignificant"))

#Label top 20 and some selected genes
#select genes to be labeled on the Volcanoplot
DEG <- DEG %>% mutate(threshold = p_val_adj < 0.01) %>% arrange(p_val_adj) %>% mutate(volcanolabels = "")

#Top20
DEG$volcanolabels[1:20] <- DEG$Gene[1:20]
DEG[which(DEG$p_val_adj==0),"p_val_adj"] <- 1e-300 #if p_val_adj==0 dot is cut off on volcano plot --> set to smallest possible number in R
GENES_toshow <- c("Baz2b","Nr1h4","Dab2","Gria3","Nrg1","Il34","Acy1","Slc22a13","Acox2","Ctnna2","Tgfbr1","Jarid2","Nr3c1","Nr3c2","Bach2","Slc9a8","Aqp1","Slc15a2","Acy3","Cdk6","Sash1","Igfbp7","Sema3c","Ppargc1b")
DEG[which(DEG$Gene%in%GENES_toshow),"volcanolabels"] <- DEG[which(DEG$Gene%in%GENES_toshow),"Gene"]

pdf("./Figure_7D.pdf")
ggplot(DEG, aes(x = avg_log2FC, y = -log10(p_val_adj), color=color,label = volcanolabels)) +
  geom_point(size = 1, alpha = 0.8, na.rm = T) + 
  geom_text_repel(aes(label=volcanolabels), max.overlaps = 1000,size=2,color="black") + 
  scale_color_manual(values =   c("#bc5f5f","#4aa4d5","darkgray"))+
  xlab("Average log2(Fold Change)") + 
  ylab("-log"[10]~"(p.adjust)") + theme_classic()+
  ggtitle("PTS3: AKI_GFP+ vs. Control")+
  theme(plot.title = element_text(hjust = 0.5, size=12),axis.text = element_text(size = 10))
dev.off()

#################################### Figure 7E & F - GO-term analysis of DEG in PTS3: AKI_GFP+ versus Control #################################### 
#Read in list of differentially expressed genes between AKI_GFP+ PTS3 and Control PTS3
DEG <- read.csv("./GFPpos_vs_Ctrl/DEG/2022_3_31_DEG_GFPposvsCtrl_clusterPTS3_logFC0_minpct0.1_padj0.01.csv")
DEG <- DEG[abs(DEG$avg_log2FC) >=0.25, ]

names(DEG)[names(DEG)=="Gene"] <- "mgi_symbol"

### Annotating Gene Symbol/Description with BiomaRt
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
DEG.anno <- getBM(values=DEG$mgi_symbol, filters = "mgi_symbol", attributes=c("ensembl_gene_id", "entrezgene_id", "description", 'mgi_symbol'), mart=mart)

#Make annotated dataframe unique by mgi.symbol
DEG.anno <- unique(DEG.anno)
DEG.anno <- subset(DEG.anno, !duplicated(subset(DEG.anno, select=c(mgi_symbol))))
DEG.anno <- merge(DEG,DEG.anno,by="mgi_symbol")

### GO term analysis with clusterprofiler
UP.anno <- DEG.anno[which(DEG.anno$avg_log2FC>=0.25),]
DOWN.anno <- DEG.anno[which(DEG.anno$avg_log2FC<=-0.25),]

GO.BP.up <- enrichGO(gene = UP.anno$entrezgene_id,OrgDb = "org.Mm.eg.db",ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10,readable=TRUE)
GO.BP.up.DF <-as.data.frame(simplify(GO.BP.up))
write.csv(GO.BP.up.DF, file="./GFPpos_vs_Ctrl/PTS3_GFPposvsCtrl_GOBP_UP_clusterprofiler.csv", row.names = F)

GO.BP.DWN <- enrichGO(gene = DOWN.anno$entrezgene_id,OrgDb = "org.Mm.eg.db",ont = "BP",pvalueCutoff  = 0.05, qvalueCutoff = 0.10, readable=TRUE)
GO.BP.DWN.DF <-as.data.frame(simplify(GO.BP.DWN))
write.csv(GO.BP.DWN.DF, file="./GFPpos_vs_Ctrl/PTS3_GFPposvsCtrl_GOBP_DOWN_clusterprofiler.csv", row.names = F)

###Make lollipop plot of top 10 GO terms
#Enriched GO terms
GO.BP.up.DF <- GO.BP.up.DF[c(1:10),c("Description","p.adjust")]
GO.BP.up.DF <- GO.BP.up.DF[order(GO.BP.up.DF$p.adjust), ]
x <- as.vector(GO.BP.up.DF$Description) # get levels
GO.BP.up.DF$Description <- factor(GO.BP.up.DF$Description, levels = rev(x)) # order GO terms by padjust order

pdf("./Figure_7E.pdf", height = 3, width=8)
ggplot(GO.BP.up.DF, aes(x = Description, y = -log10(p.adjust))) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = -log10(p.adjust)),
               color = "gray", lwd = 1.0) +
  geom_point(size = 4, color = "#bc5f5f",fill=alpha("#bc5f5f", 0.3)) +
  coord_flip() + theme_classic() + ylab("-log"[10]~"(p.adjust)") +theme(axis.title.y = element_blank(),plot.title = element_text(size=12),axis.text = element_text(size = 10)) +
  ggtitle("GO terms of 319 upregulated genes")
dev.off()

#Depleted GO terms
GO.BP.DWN.DF <- GO.BP.DWN.DF[c(1:10),c("Description","p.adjust")]
GO.BP.DWN.DF <- GO.BP.DWN.DF[order(GO.BP.DWN.DF$p.adjust), ]
x <- as.vector(GO.BP.DWN.DF$Description) # get levels
GO.BP.DWN.DF$Description <- factor(GO.BP.DWN.DF$Description, levels = rev(x)) # order GO terms by padjust order

pdf("./Figure_7F.pdf", height = 3, width=6)
ggplot(GO.BP.DWN.DF, aes(x = Description, y = -log10(p.adjust))) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = -log10(p.adjust)),
               color = "gray", lwd = 1.0) +
  geom_point(size = 4, color = "#4aa4d5",fill=alpha("#4aa4d5", 0.3)) +
  coord_flip() + theme_classic() + ylab("-log"[10]~"(p.adjust)") +theme(axis.title.y = element_blank(),plot.title = element_text(size=12),axis.text = element_text(size = 10)) +
  ggtitle("GO terms of 223 downregulated genes")
dev.off()

#################################### Figure 7G - Violin plots of differentially active motifs in PTS3: AKI_GFP+ versus Control #################################### 
#Subset PTS3 cells (AKI_GFP+ and Control) and calculate chromvar motif activity

Idents(Ki67) <- "Celltype"
DefaultAssay(Ki67) <- "peaks"
PTS3 <- subset(Ki67, idents="PTS3")
Idents(PTS3) <- "Group"
PTS3 <- subset(PTS3,idents=c("Control","AKI_GFP+"))

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object.
DefaultAssay(PTS3) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)) #this is using the core vertebrate database
motif.matrix <- CreateMotifMatrix(features = granges(PTS3), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
PTS3 <- SetAssayData(PTS3, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# calculate chromvar motif activity
PTS3 <- RunChromVAR(object = PTS3,genome = BSgenome.Mmusculus.UCSC.mm10)

# differentially active motifs between AKI_GFP+ PTS3 and Contorl PTS3
DefaultAssay(PTS3) <- "chromvar"
Idents(PTS3) <- "Group"
da_motifs <- FindMarkers(object = PTS3, ident.1 = "AKI_GFP+",ident.2 = "Control",only.pos = FALSE,test.use = 'LR',min.pct = 0.1,mean.fxn = rowMeans,fc.name = "avg_diff",logfc.threshold = 0)

#Add motif name
DefaultAssay(PTS3) <- "peaks"
da_motifs$gene <- row.names(da_motifs)
motif.names <- da_motifs$gene
da_motifs$gene <- ConvertMotifID(PTS3, id = motif.names)
da_motifs$motif <- motif.names 
da_motifs <- da_motifs[which(da_motifs$p_val_adj<0.01),]
da_motifs <- da_motifs[,c("motif","gene","pct.1","pct.2","avg_diff","p_val","p_val_adj")] #Reorder dataframe
colnames(da_motifs) <- c("motif","gene","pct.1","pct.2","motif.avg_diff","p_val","p_val_adj")
write.csv(da_motifs,"Supplemental_Table_10.csv", row.names = F)

#Plot Motif activity for Stat5a::Stat5b and NR3C1
Idents(Ki67) <- "Celltype"
levels(PTS3) <- c("AKI_GFP+","Control")

pdf("./Figure_7G.pdf",width=6)
VlnPlot(PTS3, features="MA0519.1",cols=c("#128c26","#b8bdb5"),pt.size = 0)+stat_summary(fun = median, geom='point', size = 25, colour = "darkblue", shape = 95)+ylab("chromVAR motif activity (z-score)")+
  ggtitle("MA0519.1_Stat5a::Stat5b")
VlnPlot(PTS3, features="MA0113.3",cols=c("#128c26","#b8bdb5"),pt.size = 0)+stat_summary(fun = median, geom='point', size = 25, colour = "darkblue", shape = 95)+ylab("chromVAR motif activity (z-score)")+
  ggtitle("MA0113.3_NR3C1")
dev.off()

