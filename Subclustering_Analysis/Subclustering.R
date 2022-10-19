#Title: Subclustering analysis of AKI_GFP+ nuclei of the proximal tubule and the loop of Henle (thin limb) clusters 

library(Seurat)
library(Signac)
library(dplyr)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)
library(GenomicRanges)
library(tibble)
library(harmony)
set.seed(1234)

Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

# First subset AKI_GFP+ nuclei, then subset specific cell types

head(Ki67@meta.data)
Idents(Ki67) <- "Group"

Ki67 <- subset(Ki67, idents="AKI_GFP+")
Idents(Ki67) <- "Celltype"

PTC <- subset(Ki67,idents=c("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM"))
PTC

######################################## Re-analyze the subsetted object ######################################## 
rm(Ki67)
PTC@reductions$harmony.rna <- NULL
PTC@reductions$harmony.atac<- NULL
PTC@reductions$pca<- NULL
PTC@reductions$umap.rna<- NULL
PTC@reductions$umap.atac<- NULL
PTC@reductions$wnn.umap<- NULL
PTC@reductions$lsi<- NULL

########## Integrate RNA assay #############
DefaultAssay(PTC) <-"RNA"
PTC <- NormalizeData(PTC, normalization.method = "LogNormalize", scale.factor = 10000)
PTC <- FindVariableFeatures(PTC, selection.method = "vst", nfeatures = 3000)
PTC <- ScaleData(PTC,vars.to.regress = c("percent.mt","percent.ribo","nCount_RNA"))
PTC <- RunPCA(PTC,features = VariableFeatures(object = PTC),npcs =30)
PTC <- RunHarmony(object = PTC,"orig.ident", plot_convergence = TRUE,reduction = 'pca',assay.use = 'RNA',reduction.save = "harmony.rna",project.dim = FALSE,dims.use =1:30)
PTC <- RunUMAP(PTC,dims = 1:30,reduction = 'harmony.rna',reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
PTC

########## Integrate ATAC assay #############
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(PTC) <- "peaks"
PTC <- RunTFIDF(PTC)
PTC <- FindTopFeatures(PTC, min.cutoff = 'q0')
PTC <- RunSVD(PTC)
PTC <- RunUMAP(PTC, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
PTC <- RunHarmony(object = PTC,group.by.vars="orig.ident", plot_convergence = TRUE,reduction = 'lsi',assay.use = 'peaks',reduction.save = "harmony.atac",project.dim = FALSE,dims.use =2:30)
PTC <- RunUMAP(PTC, reduction = 'harmony.atac', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
PTC

########## WNN #############
#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
PTC <- FindMultiModalNeighbors(PTC, reduction.list = list("harmony.rna", "harmony.atac"), dims.list = list(1:30, 2:30))
PTC <- RunUMAP(PTC, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
PTC <- FindClusters(PTC, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution=0.5)

p1 <- DimPlot(PTC, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE,pt.size = 1.5) + ggtitle("WNN") +NoLegend()
p2 <- DimPlot(PTC, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE,label.size = 7, repel = TRUE,pt.size =  1.5) + ggtitle("ATAC") +NoLegend()
p3 <- DimPlot(PTC, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE,label.size = 7, repel = TRUE,pt.size =  1.5) + ggtitle("RNA") +NoLegend()
p4 <- DimPlot(PTC, reduction = "wnn.umap", group.by = "SampleID", label = TRUE, repel = TRUE,pt.size = 1.5) + ggtitle("WNN") +NoLegend()
p5 <- ggplot(PTC@meta.data, aes(x=seurat_clusters, fill=Group_time)) + geom_bar(position = "fill") + theme_classic()+ theme(axis.text.x = element_text(angle = 45,hjust = 1))
p6 <- ggplot(PTC@meta.data, aes(x=seurat_clusters, fill=SampleID)) + geom_bar(position = "fill") + theme_classic() +theme(axis.text.x = element_text(angle = 45,hjust = 1))


Idents(PTC) <-"seurat_clusters"
pdf("./Subclustering_Analysis/PTC_AKIGFP+_HarmonyInt.pdf")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
VlnPlot(PTC, features = "nFeature_RNA",pt.size = 0)
VlnPlot(PTC, features = "nCount_ATAC",pt.size = 0) 
VlnPlot(PTC, features = "nCount_RNA",pt.size = 0)
VlnPlot(PTC, features = "nCount_RNA",pt.size = 0)
VlnPlot(PTC, features = "percent.ribo",pt.size = 0)
VlnPlot(PTC, features = "percent.mt",pt.size = 0)
VlnPlot(PTC, features = "TSS.enrichment",pt.size = 0)
VlnPlot(PTC, features = "atac_peak_region_fragments",pt.size = 0)
markergenes <- c("Hnf4a","Lrp2","Slc34a1","Slc5a12","Cyp2e1","Slc7a13","Havcr1","Vcam1","Ccl2","Aqp1","Slc12a1","Enox1","Slc12a3","Aqp2","Pdgfrb","Myh11","Flt1","Ptprc")
print(FeaturePlot(object = PTC,reduction = "wnn.umap", features = markergenes,pt.size =  1.5,max.cutoff = 'q95',combine=F))
dev.off()

############################## Subcluster cluster 5 (Injured PT cluster) further ############################## 

PTC <- FindSubCluster(
  PTC,
  "5",
  graph.name = "wsnn",
  subcluster.name = "sub.cluster.injured.PT",
  resolution = 0.2,
  algorithm = 3
)

p1 <- DimPlot(PTC, reduction = "wnn.umap", group.by = "sub.cluster.injured.PT", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("WNN") +NoLegend()
p2 <- DimPlot(PTC, reduction = "umap.atac", group.by = "sub.cluster.injured.PT", label = TRUE,label.size = 7, repel = TRUE) + ggtitle("ATAC") +NoLegend()
p3 <- DimPlot(PTC, reduction = "umap.rna", group.by = "sub.cluster.injured.PT", label = TRUE,label.size = 7, repel = TRUE) + ggtitle("RNA") +NoLegend()
p4 <- DimPlot(PTC, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, repel = TRUE) + ggtitle("WNN") +NoLegend()
p5 <- ggplot(PTC@meta.data, aes(x=sub.cluster.injured.PT, fill=Group_time)) + geom_bar(position = "fill") + theme_classic()+ theme(axis.text.x = element_text(angle = 45,hjust = 1))
p6 <- ggplot(PTC@meta.data, aes(x=sub.cluster.injured.PT, fill=orig.ident)) + geom_bar(position = "fill") + theme_classic() +theme(axis.text.x = element_text(angle = 45,hjust = 1))

Idents(PTC) <-"sub.cluster.injured.PT"
pdf("./Subclustering_Analysis/PTC_AKIGFP+_HarmonyInt_SUB.pdf")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
VlnPlot(PTC, features = "nFeature_RNA",pt.size = 0)
VlnPlot(PTC, features = "nCount_ATAC",pt.size = 0) 
VlnPlot(PTC, features = "nCount_RNA",pt.size = 0)
VlnPlot(PTC, features = "nCount_RNA",pt.size = 0)
VlnPlot(PTC, features = "percent.ribo",pt.size = 0)
VlnPlot(PTC, features = "percent.mt",pt.size = 0)
VlnPlot(PTC, features = "TSS.enrichment",pt.size = 0)
VlnPlot(PTC, features = "atac_peak_region_fragments",pt.size = 0)
markergenes <- c("Hnf4a","Lrp2","Slc34a1","Slc5a12","Cyp2e1","Slc7a13","Havcr1","Vcam1","Ccl2","Aqp1","Slc12a1","Enox1","Slc12a3","Aqp2","Pdgfrb","Myh11","Flt1","Ptprc")
print(FeaturePlot(object = PTC,reduction = "wnn.umap", features = markergenes,max.cutoff = 'q95',combine=F))
dev.off()

############################## Remove 1 remaining heterotypic doublet cluster and 1 cluster containing low quality nuclei ############################## 
PTC2 <- subset(PTC, idents=c("0","1","2","3","5_0","5_1","5_2","4","6"))
PTC2

p1 <- DimPlot(PTC2, reduction = "wnn.umap", group.by = "sub.cluster.injured.PT", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("WNN") +NoLegend()
p2 <- DimPlot(PTC2, reduction = "umap.atac", group.by = "sub.cluster.injured.PT", label = TRUE,label.size = 7, repel = TRUE) + ggtitle("ATAC") +NoLegend()
p3 <- DimPlot(PTC2, reduction = "umap.rna", group.by = "sub.cluster.injured.PT", label = TRUE,label.size = 7, repel = TRUE) + ggtitle("RNA") +NoLegend()
p4 <- DimPlot(PTC2, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, repel = TRUE) + ggtitle("WNN") +NoLegend()
p5 <- ggplot(PTC2@meta.data, aes(x=sub.cluster.injured.PT, fill=Group_time)) + geom_bar(position = "fill") + theme_classic()+ theme(axis.text.x = element_text(angle = 45,hjust = 1))
p6 <- ggplot(PTC2@meta.data, aes(x=sub.cluster.injured.PT, fill=orig.ident)) + geom_bar(position = "fill") + theme_classic() +theme(axis.text.x = element_text(angle = 45,hjust = 1))

Idents(PTC2) <-"sub.cluster.injured.PT"
pdf("./Subclustering_Analysis/PTC_AKIGFP+_HarmonyInt_SUB_FINAL.pdf")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
VlnPlot(PTC2, features = "nFeature_RNA",pt.size = 0)
VlnPlot(PTC2, features = "nCount_ATAC",pt.size = 0) 
VlnPlot(PTC2, features = "nCount_RNA",pt.size = 0)
VlnPlot(PTC2, features = "nCount_RNA",pt.size = 0)
VlnPlot(PTC2, features = "percent.ribo",pt.size = 0)
VlnPlot(PTC2, features = "percent.mt",pt.size = 0)
VlnPlot(PTC2, features = "TSS.enrichment",pt.size = 0)
VlnPlot(PTC2, features = "atac_peak_region_fragments",pt.size = 0)
markergenes <- c("Hnf4a","Lrp2","Slc34a1","Slc5a12","Cyp2e1","Slc7a13","Havcr1","Vcam1","Ccl2","Aqp1","Slc12a1","Enox1","Slc12a3","Aqp2","Pdgfrb","Myh11","Flt1","Ptprc")
print(FeaturePlot(object = PTC2,reduction = "wnn.umap", features = markergenes,max.cutoff = 'q95',combine=F))
dev.off()

#Rename the clusters
levels(PTC2) <- c("2","3","0","1","5_1","5_0","5_2","4","6")
new.cluster.ids <- c("PTS1","PTS1-S2","PTS2","PTS3","Injured PTS1-S2","Injured PTS3","FR-PTC","LOH-TL-C","LOH-TL-JM")
names(new.cluster.ids) <- levels(PTC2)
PTC2 <- RenameIdents(PTC2, new.cluster.ids)
PTC2@meta.data$Celltype <- PTC2@active.ident

p1 <- DimPlot(PTC2, reduction = "wnn.umap", group.by = "Celltype", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") +NoLegend()
p2 <- DimPlot(PTC2, reduction = "umap.atac", group.by = "Celltype", label = TRUE,label.size =5, repel = TRUE) + ggtitle("ATAC") +NoLegend()
p3 <- DimPlot(PTC2, reduction = "umap.rna", group.by = "Celltype", label = TRUE,label.size = 5, repel = TRUE) + ggtitle("RNA") +NoLegend()
p5 <- ggplot(PTC2@meta.data, aes(x=Celltype, fill=Group_time)) + geom_bar(position = "fill") + theme_classic()+ theme(axis.text.x = element_text(angle = 45,hjust = 1))
p6 <- ggplot(PTC2@meta.data, aes(x=Celltype, fill=orig.ident)) + geom_bar(position = "fill") + theme_classic() +theme(axis.text.x = element_text(angle = 45,hjust = 1))


Idents(PTC2) <-"Celltype"
pdf("./Subclustering_Analysis/PTC_AKIGFP+_HarmonyInt_SUB_FINAL_Celltype.pdf")
print(p1)
print(p2)
print(p3)
print(p5)
print(p6)
VlnPlot(PTC2, features = "nFeature_RNA",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "nCount_ATAC",pt.size = 0) +theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "nCount_RNA",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "nCount_RNA",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "percent.ribo",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "percent.mt",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "TSS.enrichment",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
VlnPlot(PTC2, features = "atac_peak_region_fragments",pt.size = 0)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

#Add new cell type names to meta data
PTC2@meta.data$Group_time_Celltype <- paste0(PTC2@meta.data$Group_time,sep="_",PTC2@meta.data$Celltype)
PTC2@meta.data$Group_Celltype <- paste0(PTC2@meta.data$Group,sep="_",PTC2@meta.data$Celltype)

saveRDS(PTC2,"./Subclustering_Analysis/PTC_FINAL_SUB.Rds")
saveRDS(PTC,"./Subclustering_Analysis/PTC_all_SUB.Rds")