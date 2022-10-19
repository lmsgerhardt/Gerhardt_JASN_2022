#Title: Clustering

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)
library(tibble)
library(harmony)
set.seed(1234)

Ki67 <- readRDS("./Ki67_filtered.Rds")
Ki67

########## Integrate RNA assay #############
DefaultAssay(Ki67) <-"RNA"
Ki67 <- NormalizeData(Ki67, normalization.method = "LogNormalize", scale.factor = 10000)
Ki67 <- FindVariableFeatures(Ki67, selection.method = "vst", nfeatures = 3000)
Ki67 <- ScaleData(Ki67,vars.to.regress = c("percent.mt","percent.ribo","nCount_RNA"))
Ki67 <- RunPCA(Ki67,features = VariableFeatures(object = Ki67),npcs =30)
Ki67 <- RunHarmony(
  object = Ki67,
  "orig.ident", 
  plot_convergence = TRUE,
  reduction = 'pca',
  assay.use = 'RNA',
  reduction.save = "harmony.rna",
  project.dim = FALSE,dims.use =1:30 
)

Ki67 <- RunUMAP(Ki67,dims = 1:30,reduction = 'harmony.rna', reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Ki67

########## Integrate ATAC assay #############
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(Ki67) <- "ATAC"
Ki67 <- RunTFIDF(Ki67)
Ki67 <- FindTopFeatures(Ki67, min.cutoff = 'q0')
Ki67 <- RunSVD(Ki67)
str(Ki67)
Ki67 <- RunHarmony(
  object = Ki67,
  group.by.vars=c("orig.ident"), 
  plot_convergence = TRUE,
  reduction = 'lsi',
  assay.use = 'ATAC',
  reduction.save = "harmony.atac",
  project.dim = FALSE,dims.use =2:30
)

Ki67 <- RunUMAP(Ki67, reduction = 'harmony.atac', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

########## WNN #############
#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

Ki67 <- FindMultiModalNeighbors(Ki67, reduction.list = list("harmony.rna", "harmony.atac"), dims.list = list(1:30, 2:30))
Ki67 <- RunUMAP(Ki67, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Ki67 <- FindClusters(Ki67, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution=0.8)

p1 <- DimPlot(Ki67, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p12 <- DimPlot(Ki67, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p13 <- DimPlot(Ki67, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(Ki67, reduction = "wnn.umap", group.by = "Group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p3 <- ggplot(Ki67@meta.data, aes(x=seurat_clusters, fill=Group)) + geom_bar(position = "fill") + theme_classic() 
p4 <- ggplot(Ki67@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position = "fill") + theme_classic()
p5 <- ggplot(Ki67@meta.data, aes(x=seurat_clusters, fill=Group_time)) + geom_bar(position = "fill") + theme_classic()

table(Idents(Ki67))

pdf("Plots_Clustering_wnn_30PCs_res0.8.pdf",width=8)
print(p1)
print(p12)
print(p13)
print(p2)
print(p3)
print(p4)
print(p5)
Idents(Ki67) <- "seurat_clusters"
VlnPlot(Ki67, features = "nFeature_RNA", pt.size = 0)
VlnPlot(Ki67, features = "nCount_ATAC", pt.size = 0) 
VlnPlot(Ki67, features = "nCount_RNA", pt.size = 0)
VlnPlot(Ki67, features = "percent.ribo", pt.size = 0)
VlnPlot(Ki67, features = "percent.mt", pt.size = 0)
VlnPlot(Ki67, features = "TSS.enrichment", pt.size = 0)
VlnPlot(Ki67, features = "atac_peak_region_fragments", pt.size = 0)
markergenes <- c("Nphs2","Hnf4a","Lrp2","Slc34a1","Slc5a12","Cyp2e1","Slc7a13","Havcr1","Vcam1","Aqp1","Slc12a1","Enox1","Slc12a3","Calb1","Scnn1g","Aqp2","Atp6v1g3","Kit","Slc26a4","Eln","Pdgfrb","Myh11","Flt1","Ptprc")
print(FeaturePlot(object = Ki67,reduction = "wnn.umap", features = markergenes,pt.size = 0.1,max.cutoff = 'q95',combine=F))
dev.off()

saveRDS(Ki67,"Ki67_clustered.Rds")

########## Name cell types #############
Idents(Ki67) <- "seurat_clusters"
my_levels <- c("25","2","3","4","0","1","16","15","17","10","5","6","18","9","8","7","12","19","23","14","13","11","22","20",
               "21","24","26")
levels(Ki67) <- my_levels

new.cluster.ids <- c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
                     "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell", "Doublets_1","Doublets_2","Small unknown")

names(new.cluster.ids) <- levels(Ki67)
Ki67 <- RenameIdents(Ki67, new.cluster.ids)
Ki67@meta.data$Celltype <- Ki67@active.ident

########## Remove the two remaining doublet clusters and the small cluster containing < 10 nuclei  #############
Idents(Ki67) <- "Celltype"
Ki67_FINAL <- subset(Ki67, idents=c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","FR-PTC","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
                                    "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell"))

p1 <- DimPlot(Ki67, reduction = "wnn.umap", group.by = "Celltype", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN") +NoLegend()
p2 <- DimPlot(Ki67, reduction = "umap.atac", group.by = "Celltype", label = FALSE, repel = TRUE) + ggtitle("ATAC") +NoLegend()
p3 <- DimPlot(Ki67, reduction = "umap.rna", group.by = "Celltype", label = FALSE, repel = TRUE) + ggtitle("RNA") +NoLegend()

pdf("Umap_plots_celltype_FINAL.pdf")
print(p1)
print(p2)
print(p3)
dev.off()

########## Add metadata on Celltype per group and time point #############
Ki67@meta.data$Group_time_Celltype <- paste0(Ki67@meta.data$Group_time, sep="_",Ki67@meta.data$Celltype)

saveRDS(Ki67_FINAL, "Ki67_FINAL.Rds")

