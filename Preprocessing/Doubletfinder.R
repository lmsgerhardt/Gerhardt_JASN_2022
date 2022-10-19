#Title: Doubletfinder to remove heterotypic doublets
#Some parts adapted from Parker Wilson's github: p4rkerw/Muto_Wilson_NComm_2020 (https://github.com/p4rkerw/Muto_Wilson_NComm_2020/blob/master/snRNA_prep/seurat_rna_process.R [last accessed 10/17/22])

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
library(DoubletFinder)
library(tibble)
set.seed(1234)

#The following estimated doublet rates will be used:
#10% in samples with <11000 nuclei, i.e. in samples Ctrl_6months, AKI_GFP+_4weeks_1, AKI_GFP+_4weeks_2, AKI_GFP+_4weeks_3, AKI_GFP+_6months_1, AKI_GFP+_6months_2, AKI_GFP+_6months_3, AKI_GFP+_6months_4, AKI_GFP-_6months
#15% in samples with 11000-13000, i.e. in samples Ctrl_4weeks_1, AKI_GFP-_4weeks
#25% in one sample with 16093 nuclei (sample Ctrl_4weeks_2)

setwd("/Ki67_analysis/Doubletfinder")

# Identify doublets separately for all samples with an estimated doublet rate of 10%
SAMPLES <- c("Ctrl_6months", "AKI_GFP+_4weeks_1", "AKI_GFP+_4weeks_2", "AKI_GFP+_4weeks_3", "AKI_GFP+_6months_1", "AKI_GFP+_6months_2", "AKI_GFP+_6months_3", "AKI_GFP+_6months_4", "AKI_GFP-_6months")

for(i in 1:length(SAMPLES)){
  dir.create(SAMPLES[[i]])
  setwd(SAMPLES[[i]])
  
  ## Read in the data after ambient RNA removal with SoupX and create Seurat object with basic filtering of low quality nuclei
  SObj.data <- Read10X(data.dir = paste0("/Ki67_analysis/SoupX/",SAMPLES[[i]],"/strainedCounts"))
  SObj <-  CreateSeuratObject(counts = SObj.data, min.cells = 3, min.features = 200)
  head(SObj@meta.data)
  
  ## Perform PCA and clustering
  SObj <- NormalizeData(SObj)
  SObj <- ScaleData(SObj, features=row.names(SObj))
  SObj <- FindVariableFeatures(SObj, selection.method = "vst", nfeatures = 3000)
  SObj  <- RunPCA(SObj)
  ElbowPlot(SObj,ndims = 50, reduction = "pca")
  SObj <- FindNeighbors(SObj, dims = 1:20)
  SObj <- RunUMAP(SObj, dims = 1:20)
  DimPlot(SObj)
  SObj <- FindClusters(SObj, algorithm = 3, verbose = FALSE, resolution=0.8)
  DimPlot(SObj)
  
  ## Doublet removal with the assumption that doublets represent 10% of cells.
  # pK 
  sweep.res.list_kidney <- paramSweep_v3(SObj, PCs = 1:20, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(SObj, PCs = 1:20, pN = 0.25, pK = pK,
                                      nExp = round(0.1*nrow(SObj@meta.data)), 
                                      reuse.pANN = FALSE, sct = F)
  
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  pdf(paste0(SAMPLES[[i]],"_10perct.pdf"))
  print(p1) 
  print(p2)
  print(p3)
  FeaturePlot(SObj,features = "Slc12a1")
  FeaturePlot(SObj,features = "Slc12a3")
  FeaturePlot(SObj,features = "Lrp2")
  FeaturePlot(SObj,features = "Flt1")
  FeaturePlot(SObj,features = "Ptprc")
  dev.off()
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes_10perct <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  table(df_doublet_barcodes_10perct$V2) 
  
  #save
  write.csv(df_doublet_barcodes_10perct,paste0(SAMPLES[[i]],"_df_doublet_barcodes_10perct.csv"),row.names = F)
  
  setwd('..')
}

# Identify doublets separately for all samples with an estimated doublet rate of 15%
SAMPLES <- c("Ctrl_4weeks_1", "AKI_GFP-_4weeks")

for(i in 1:length(SAMPLES)){
  dir.create(SAMPLES[[i]])
  setwd(SAMPLES[[i]])
  
  ## Read in the data after ambient RNA removal with SoupX and create Seurat object with basic filtering of low quality nuclei
  SObj.data <- Read10X(data.dir = paste0("/Ki67_analysis/SoupX/",SAMPLES[[i]],"/strainedCounts"))
  SObj <-  CreateSeuratObject(counts = SObj.data, min.cells = 3, min.features = 200)
  head(SObj@meta.data)
  
  ## Perform PCA and clustering
  SObj <- NormalizeData(SObj)
  SObj <- ScaleData(SObj, features=row.names(SObj))
  SObj <- FindVariableFeatures(SObj, selection.method = "vst", nfeatures = 3000)
  SObj  <- RunPCA(SObj)
  ElbowPlot(SObj,ndims = 50, reduction = "pca")
  SObj <- FindNeighbors(SObj, dims = 1:20)
  SObj <- RunUMAP(SObj, dims = 1:20)
  DimPlot(SObj)
  SObj <- FindClusters(SObj, algorithm = 3, verbose = FALSE, resolution=0.8)
  DimPlot(SObj)
  
  ## Doublet removal with the assumption that doublets represent 15% of cells.
  # pK 
  sweep.res.list_kidney <- paramSweep_v3(SObj, PCs = 1:20, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(SObj, PCs = 1:20, pN = 0.25, pK = pK,
                                      nExp = round(0.15*nrow(SObj@meta.data)), 
                                      reuse.pANN = FALSE, sct = F)
  
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  pdf(paste0(SAMPLES[[i]],"_15perct.pdf"))
  print(p1) 
  print(p2)
  print(p3)
  FeaturePlot(SObj,features = "Slc12a1")
  FeaturePlot(SObj,features = "Slc12a3")
  FeaturePlot(SObj,features = "Lrp2")
  FeaturePlot(SObj,features = "Flt1")
  FeaturePlot(SObj,features = "Ptprc")
  dev.off()
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes_15perct <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  table(df_doublet_barcodes_15perct$V2) 
  
  #save
  write.csv(df_doublet_barcodes_15perct,paste0(SAMPLES[[i]],"_df_doublet_barcodes_15perct.csv"),row.names = F)
  
  setwd('..')
}

# Identify doublets separately for the sample with an estimated doublet rate of 25%
SAMPLES <- "Ctrl_4weeks_2"

for(i in 1:length(SAMPLES)){
  dir.create(SAMPLES[[i]])
  setwd(SAMPLES[[i]])
  
  ## Read in the data after ambient RNA removal with SoupX and create Seurat object with basic filtering of low quality nuclei
  SObj.data <- Read10X(data.dir = paste0("/Ki67_analysis/SoupX/",SAMPLES[[i]],"/strainedCounts"))
  SObj <-  CreateSeuratObject(counts = SObj.data, min.cells = 3, min.features = 200)
  head(SObj@meta.data)
  
  ## Perform PCA and clustering
  SObj <- NormalizeData(SObj)
  SObj <- ScaleData(SObj, features=row.names(SObj))
  SObj <- FindVariableFeatures(SObj, selection.method = "vst", nfeatures = 3000)
  SObj  <- RunPCA(SObj)
  ElbowPlot(SObj,ndims = 50, reduction = "pca")
  SObj <- FindNeighbors(SObj, dims = 1:20)
  SObj <- RunUMAP(SObj, dims = 1:20)
  DimPlot(SObj)
  SObj <- FindClusters(SObj, algorithm = 3, verbose = FALSE, resolution=0.8)
  DimPlot(SObj)
  
  ## Doublet removal with the assumption that doublets represent 25% of cells.
  # pK 
  sweep.res.list_kidney <- paramSweep_v3(SObj, PCs = 1:20, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(SObj, PCs = 1:20, pN = 0.25, pK = pK,
                                      nExp = round(0.25*nrow(SObj@meta.data)), 
                                      reuse.pANN = FALSE, sct = F)
  
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  pdf(paste0(SAMPLES[[i]],"_25perct.pdf"))
  print(p1) 
  print(p2)
  print(p3)
  FeaturePlot(SObj,features = "Slc12a1")
  FeaturePlot(SObj,features = "Slc12a3")
  FeaturePlot(SObj,features = "Lrp2")
  FeaturePlot(SObj,features = "Flt1")
  FeaturePlot(SObj,features = "Ptprc")
  dev.off()
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes_25perct <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  table(df_doublet_barcodes_25perct$V2) 
  
  #save
  write.csv(df_doublet_barcodes_25perct,paste0(SAMPLES[[i]],"_df_doublet_barcodes_25perct.csv"),row.names = F)
  
  setwd('..')
}

########################################### Aggregate all information into one singlet / doublet matrix ###########################################

Ctrl_4weeks_1 <- read.csv("./Ctrl_4weeks_1/Ctrl_4weeks_1_df_doublet_barcodes_15perct.csv",header=T)
AKI_GFPpos_4weeks_1 <- read.csv("./AKI_GFP+_4weeks_1/AKI_GFP+_4weeks_1_df_doublet_barcodes_10perct.csv",header=T)
Ctrl_4weeks_2 <- read.csv("./Ctrl_4weeks_2/Ctrl_4weeks_2_df_doublet_barcodes_25perct.csv",header=T)
AKI_GFPpos_4weeks_2 <- read.csv("./AKI_GFP+_4weeks_2/AKI_GFP+_4weeks_2_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPpos_6months_1 <- read.csv("./AKI_GFP+_6months_1/AKI_GFP+_6months_1_df_doublet_barcodes_10perct.csv",header=T)
Ctrl_6months <- read.csv("./Ctrl_6months/Ctrl_6months_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPpos_6months_2 <- read.csv("./AKI_GFP+_6months_2/AKI_GFP+_6months_2_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPpos_6months_3 <- read.csv("./AKI_GFP+_6months_3/AKI_GFP+_6months_3_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPneg_6months <- read.csv("./AKI_GFP-_6months/AKI_GFP-_6months_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPpos_6months_4 <- read.csv("./AKI_GFP+_6months_4/AKI_GFP+_6months_4_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPpos_4weeks_3 <- read.csv("./AKI_GFP+_4weeks_3/AKI_GFP+_4weeks_3_df_doublet_barcodes_10perct.csv",header=T)
AKI_GFPneg_4weeks <- read.csv("./AKI_GFP-_4weeks/AKI_GFP-_4weeks_df_doublet_barcodes_15perct.csv",header=T)

#Rename V1 to fit naming in the aggregated seurat object
AKI_GFPpos_4weeks_1$V1 <- paste0(substr(AKI_GFPpos_4weeks_1$V1,1,17),rep("2",length(AKI_GFPpos_4weeks_1$V1)))
Ctrl_4weeks_2$V1 <- paste0(substr(Ctrl_4weeks_2$V1,1,17),rep("3",length(Ctrl_4weeks_2$V1)))
AKI_GFPpos_4weeks_2$V1 <- paste0(substr(AKI_GFPpos_4weeks_2$V1,1,17),rep("4",length(AKI_GFPpos_4weeks_2$V1)))
AKI_GFPpos_6months_1$V1 <- paste0(substr(AKI_GFPpos_6months_1$V1,1,17),rep("5",length(AKI_GFPpos_6months_1$V1)))
Ctrl_6months$V1 <- paste0(substr(Ctrl_6months$V1,1,17),rep("6",length(Ctrl_6months$V1)))
AKI_GFPpos_6months_2$V1 <- paste0(substr(AKI_GFPpos_6months_2$V1,1,17),rep("7",length(AKI_GFPpos_6months_2$V1)))
AKI_GFPpos_6months_3$V1 <- paste0(substr(AKI_GFPpos_6months_3$V1,1,17),rep("8",length(AKI_GFPpos_6months_3$V1)))
AKI_GFPneg_6months$V1 <- paste0(substr(AKI_GFPneg_6months$V1,1,17),rep("9",length(AKI_GFPneg_6months$V1)))
AKI_GFPpos_6months_4$V1 <- paste0(substr(AKI_GFPpos_6months_4$V1,1,17),rep("10",length(AKI_GFPpos_6months_4$V1)))
AKI_GFPpos_4weeks_3$V1 <- paste0(substr(AKI_GFPpos_4weeks_3$V1,1,17),rep("11",length(AKI_GFPpos_4weeks_3$V1)))
AKI_GFPneg_4weeks$V1 <- paste0(substr(AKI_GFPneg_4weeks$V1,1,17),rep("12",length(AKI_GFPneg_4weeks$V1)))

mymatrix <- rbind(Ctrl_4weeks_1,AKI_GFPpos_4weeks_1,Ctrl_4weeks_2,AKI_GFPpos_4weeks_2,AKI_GFPpos_6months_1,Ctrl_6months,AKI_GFPpos_6months_2,
                  AKI_GFPpos_6months_3,AKI_GFPneg_6months,AKI_GFPpos_6months_4,AKI_GFPpos_4weeks_3,AKI_GFPneg_4weeks)
dim(mymatrix) 
table(mymatrix$V2)

write.csv(mymatrix,"Ki67_singlets_doublets_matrix.csv",row.names = F)