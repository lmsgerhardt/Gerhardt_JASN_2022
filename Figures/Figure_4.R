#Title: Plots for Figure 4

library(Seurat)
library(Signac)
library(ggplot2)
set.seed(1234)

setwd()

#Read in the final R object of AKI_GFP+ ("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM") subclustering analysis 
PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")
PTC

#################################### Figure 4A - Umap plot #################################### 
Idents(PTC) <- "Celltype"
levels(PTC) <- c("PTS1","PTS1-S2","PTS2","PTS3","Injured PTS1-S2","Injured PTS3","FR-PTC","LOH-TL-C","LOH-TL-JM")

pdf("Figure_4A.pdf")
DimPlot(PTC, reduction = "wnn.umap", label = TRUE, label.size = 7, repel = TRUE) + NoLegend()
dev.off()

#################################### Figure 4B - Dot plot of marker genes and stacked bar plot #################################### 
#Figure 4B.1 - Dot plot of marker genes
DefaultAssay(PTC) <- "RNA"

Markergenes <- c("Hnf4a","Lrp2","Slc34a1","Slc5a12","Nox4","Slc22a8","Cyp2e1","Slco1a1","Cyp7b1","Slc22a19","Bcat1","Grid1","Slc7a12","Kcnip4","Dock10","Vcam1","Sytl2","Spp1","Dcdc2a","Pdgfd","Havcr1","Ccl2","Krt20","Cp","Jag1","Slit3","Slc14a2","Cdh6","Cdh13","Apela")

levels(PTC) <- rev(levels(PTC))

pdf("Figure_4B1.pdf",height=4,width=8.1)
DotPlot(PTC, features = Markergenes, col.min = 0, cols = c("darkgrey", "black"),dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90,vjust=0.5))
dev.off()

#Figure 4B.2 - Stacked bar plot normalized by nuclei number in group
#Usage: Make a stacked bar plot showing the composition of clusters by group and time point (AKI_GFP+_4weeks, AKI_GF+_6months), normalized for the total number of nuclei in each group
#Extract meta data from the seurat object

Meta <- PTC@meta.data
Celltype <-c("PTS1","PTS1-S2","PTS2","PTS3","Injured PTS1-S2","Injured PTS3","FR-PTC","LOH-TL-C","LOH-TL-JM")
Group_time <- unique(Meta$Group_time)
myDF_Group_time <- data.frame("Composition"=c(),"Celltype"=c(),"Group"=c())

#Read in Ki67 seurat object to get the total number of nuclei per group and timepoint.
Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67
Meta_allclusters <- Ki67@meta.data

# Normalize for total number of nuclei in each group
for(i in 1:length(Group_time)){
  CT <- Meta[which(Meta$Group_time==Group_time[[i]]),"Celltype"]
  allclusters <- Meta_allclusters[which(Meta_allclusters$Group_time==Group_time[[i]]),"seurat_clusters"]
  CT <- table(CT)/length(allclusters)
  myDF_CT <- data.frame("Composition"=as.vector(CT),"Celltype"=names(CT),"Group"=rep(Group_time[[i]],length(CT)))
  myDF_Group_time <- rbind(myDF_Group_time, myDF_CT)
}

# make stacked bar plot
myDF_Group_time$Celltype <- factor(myDF_Group_time$Celltype, rev(Celltype))
myDF_Group_time$Group <- factor(myDF_Group_time$Group, c("AKI_GFP+_6months","AKI_GFP+_4weeks"))

pdf("Figure_4B2.pdf",height=5)
ggplot(myDF_Group_time, aes(fill=Group, x=Composition, y=Celltype)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity")+ scale_fill_manual(values=c("#1B842C","#61B329"))+
  theme_classic()+ ggtitle("Composition of clusters\n(normalized for cell number in group)")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()