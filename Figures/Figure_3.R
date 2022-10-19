#Title: Plots for Figure 3

library(Seurat)
library(Signac)
library(ggplot2)
set.seed(1234)

setwd()

#Read in the final R object (contains RNA, chromVAR and peaks assay (peaks called using MACS2))
Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

#################################### Figure 3B - Umap plots #################################### 
##Generate dimplots - atac, rna and wnn
#by Celltype
p1 <- DimPlot(Ki67, reduction = "wnn.umap", group.by = "Celltype", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN") +NoLegend()
p2 <- DimPlot(Ki67, reduction = "umap.atac", group.by = "Celltype", label = FALSE, repel = TRUE) + ggtitle("ATAC") +NoLegend()
p3 <- DimPlot(Ki67, reduction = "umap.rna", group.by = "Celltype", label = FALSE, repel = TRUE) + ggtitle("RNA") +NoLegend()

pdf("Figure_3B.pdf")
print(p1)
print(p2)
print(p3)
dev.off()

#################################### Figure 3C - Dot plot of marker genes #################################### 
DefaultAssay(Ki67) <- "RNA"
Markergenes <- c("Nphs1","Nphs2","Hnf4a","Slc5a12","Cyp2e1","Slco1a1","Cyp7b1","Vcam1","Kcnip4","Jag1","Cdh6","Cdh13",
                  "Enox1","Slc12a1","Slc12a3","Calb1","Egfem1","Aqp2","Aqp4",
                  "Kit","Slc4a9","Slc26a4","Upk1b","Cfh","Pdgfrb","Flt1","Ptprc","Cd74","Wdfy4","Igha")

levels(Ki67) <- rev(levels(Ki67))
pdf("Figure_3C.pdf",width=8.5,height=6)
DotPlot(Ki67, features = Markergenes, col.min = 0, cols = c("darkgrey", "black"),dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90,vjust=0.5))
dev.off()

#################################### Figure 3D - Stacked bar plot normalized for nuclei number in group #################################### 
#Usage: Make a stacked bar plot showing the composition of clusters by group and time point, normalized for the total number of nuclei in each group

# Extract metadata from Seurat object
Meta <- Ki67@meta.data
Group_time <- unique(Meta$Group_time)
myDF_Group_time <- data.frame("Composition"=c(),"Celltype"=c(),"Group"=c())

# Normalize for total number of nuclei in each group
for(i in 1:length(Group_time)){
  CT <- Meta[which(Meta$Group_time==Group_time[[i]]),"Celltype"]
  CT <- table(CT)/length(CT)
  myDF_CT <- data.frame("Composition"=as.vector(CT),"Celltype"=names(CT),"Group"=rep(Group_time[[i]],length(CT)))
  myDF_Group_time <- rbind(myDF_Group_time, myDF_CT)
}

# Make Stacked bar plot
Celltype <- c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
              "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell")
myDF_Group_time$Celltype <- factor(myDF_Group_time$Celltype, rev(Celltype))
myDF_Group_time$Group <- factor(myDF_Group_time$Group, c("AKI_GFP+_4weeks","AKI_GFP+_6months","AKI_GFP-_4weeks","AKI_GFP-_6months","Ctrl_4weeks","Ctrl_6months"))

pdf("Figure_3D.pdf")
ggplot(myDF_Group_time, aes(fill=Group, x=Composition, y=Celltype)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity")+ scale_fill_manual(values=c("#61B329","#1B842C","#6394ec","#1e65e4","#D3D3D3","#bdbdbd"))+
  theme_classic()+ ggtitle("Composition of clusters\n(normalized for cell number in group)")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#################################### Figure 3E - Coverage plot of Hnf4a #################################### 

levels(Ki67) <- c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
                              "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell")

pdf("Figure_3E.pdf",width=4, height=5)
CoveragePlot(Ki67, region = "Hnf4a",assay = 'peaks', peaks = FALSE) 
dev.off()

