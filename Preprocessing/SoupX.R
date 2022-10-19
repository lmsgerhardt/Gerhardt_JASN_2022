#Title: remove ambient RNA on each sample separately using SoupX

library(Seurat)
library(SoupX)
library(ggplot2)
library(DropletUtils)

#Make vector with all sample names; then run for loop across all samples (Cellranger out put files are stored in Cellranger_outs directory)
SAMPLES <- c("Ctrl_4weeks_1","Ctrl_4weeks_2","AKI_GFP+_4weeks_1","AKI_GFP+_4weeks_2","AKI_GFP+_4weeks_3","AKI_GFP-_4weeks",
             "Ctrl_6months","AKI_GFP+_6months_1","AKI_GFP+_6months_2","AKI_GFP+_6months_3","AKI_GFP+_6months_4",
             "AKI_GFP-_6months")

for(i in 1:length(SAMPLES)){
  setwd("./SoupX")
  
  #Read in the filtered and raw gene expression matrix (10X output files)
  #the 10x hdf5 file contains both data types --> only use gene expression 
  filtered.10X <- Read10X_h5(paste0("./Cellranger_outs/",SAMPLES[[i]],"/filtered_feature_bc_matrix.h5"))
  rna_counts_filt <- filtered.10X$`Gene Expression`
  rm(filtered.10X)
  
  raw.10X <- Read10X_h5(paste0("./Cellranger_outs/",SAMPLES[[i]],"/raw_feature_bc_matrix.h5"))
  rna_counts_raw <- raw.10X$`Gene Expression`
  rm(raw.10X)
  
  #################################################### BASIC SEURAT WORKFLOW TO GENERATE CLUSTERING #################################################### 
  dir.create(SAMPLES[[i]])
  setwd(SAMPLES[[i]])
  
  #Create Seurat object without filtering
  SObj <-  CreateSeuratObject(counts = rna_counts_filt, min.cells = 0, min.features = 0)
  head(SObj@meta.data)
  SObj[["percent.mt"]] <- PercentageFeatureSet(SObj, pattern = "^mt-")
  C<-GetAssayData(object = SObj, slot = "counts")
  rb.genes <- rownames(SObj)[grep("^Rp[sl]",rownames(SObj))]
  percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
  SObj <- AddMetaData(SObj, percent.ribo, col.name = "percent.ribo")
  
  # Visualize QC metrics as a violin plot
  pdf(paste0(SAMPLES[[i]],"_QCplot.pdf"))
  print(VlnPlot(SObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), pt.size=0, ncol = 4))
  dev.off()
  
  SObj <- NormalizeData(SObj)
  SObj <- ScaleData(SObj, features=row.names(SObj))
  SObj <- FindVariableFeatures(SObj, selection.method = "vst", nfeatures = 3000)
  SObj  <- RunPCA(SObj)
  SObj <- FindNeighbors(SObj, dims = 1:20)
  SObj <- RunUMAP(SObj, dims = 1:20)
  SObj <- FindClusters(SObj, algorithm = 3, verbose = FALSE, resolution=0.5)
  
  #save plots
  pdf(paste0(SAMPLES[[i]],"_Plots_beforeSoupX.pdf"))
  print(ElbowPlot(SObj,ndims = 50, reduction = "pca"))
  print(DimPlot(SObj))
  print(FeaturePlot(SObj,features=c("Slc34a1","Slc5a12","Slc12a1","Aqp1","Aqp2","Slc12a3") #feature plots of some abundant genes to compare distribution before and after ambient RNA removal
              ,order=T, combine=F))
  dev.off()
  
  #################################################### RUN SOUPX TO REMOVE AMBIENT RNA #################################################### 
  Clusternames <- SObj@meta.data$seurat_clusters
  names(Clusternames) <- row.names(SObj@meta.data)
  #Create Soup Channel object
  toc <- rna_counts_filt
  dim(toc) #
  tod <- rna_counts_raw
  sc = SoupChannel(tod, toc)
  
  #Profiling the soup
  sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc = estimateSoup(sc)
  
  sc = setClusters(sc, setNames(Clusternames, names(Clusternames)))
  
  #Measuring the contamination fraction
  sc = autoEstCont(sc)
  
  #Correcting expression profile
  out = adjustCounts(sc)
  
  #investigating changes in expression
  #We can get a sense for what has been the most strongly decreased by looking at the fraction of cells that were non-zero 
  #now set to zero after correction.
  cntSoggy = rowSums(sc$toc > 0)
  cntStrained = rowSums(out > 0)
  mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
  mostZeroed
  
  #If on the other hand we focus on genes for which there is a quantitative difference,we find genes associated with metabolism and translation. 
  tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
  
  #Save the cleaned up matrix
  DropletUtils:::write10xCounts("./strainedCounts", out)
  
  ########################################### Remake Seurat object from cleaned matrix to see if there is a difference ###########################################
  SObj.data <- Read10X(data.dir = "./strainedCounts")
  SObj <-  CreateSeuratObject(counts = SObj.data, min.cells = 0, min.features = 0)
  head(SObj@meta.data)
  SObj[["percent.mt"]] <- PercentageFeatureSet(SObj, pattern = "^mt-")
  C<-GetAssayData(object = SObj, slot = "counts")
  rb.genes <- rownames(SObj)[grep("^Rp[sl]",rownames(SObj))]
  percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
  SObj <- AddMetaData(SObj, percent.ribo, col.name = "percent.ribo")
  
  SObj <- NormalizeData(SObj)
  SObj <- ScaleData(SObj, features=row.names(SObj))
  SObj <- FindVariableFeatures(SObj, selection.method = "vst", nfeatures = 3000)
  SObj  <- RunPCA(SObj)
  SObj <- FindNeighbors(SObj, dims = 1:20)
  SObj <- RunUMAP(SObj, dims = 1:20)
  SObj <- FindClusters(SObj, algorithm = 3, verbose = FALSE, resolution=0.5)
  pdf(paste0(SAMPLES[[i]],"_Plots_afterSoupX.pdf"))
  print(DimPlot(SObj))
  print(FeaturePlot(SObj,features=c("Slc34a1","Slc5a12","Slc12a1","Aqp1","Aqp2","Slc12a3") 
              ,order=T, combine=F))
  dev.off()
  
  setwd('..')
}

########################################### Aggregate all matrices in which ambient RNA was removed into one  matrix ###########################################

setwd("./SoupX")

# Load each data set and create a merged Seurat object (order and numbering as in Ki67_libraries for cellranger aggr)
Ctrl_4weeks_1 <- Read10X("./Ctrl_4weeks_1/strainedCounts")
AKI_GFPpos_4weeks_1 <- Read10X("./AKI_GFP+_4weeks_1/strainedCounts")
Ctrl_4weeks_2 <- Read10X("./Ctrl_4weeks_2/strainedCounts")
AKI_GFPpos_4weeks_2 <- Read10X("./AKI_GFP+_4weeks_2/strainedCounts")
AKI_GFPpos_6months_1 <- Read10X("./AKI_GFP+_6months_1/strainedCounts")
Ctrl_6months <- Read10X("./Ctrl_6months/strainedCounts")
AKI_GFPpos_6months_2 <- Read10X("./AKI_GFP+_6months_2/strainedCounts")
AKI_GFPpos_6months_3 <- Read10X("./AKI_GFP+_6months_3/strainedCounts")
AKI_GFPneg_6months <- Read10X("./AKI_GFP-_6months/strainedCounts")
AKI_GFPpos_6months_4 <- Read10X("./AKI_GFP+_6months_4/strainedCounts")
AKI_GFPpos_4weeks_3 <- Read10X("./AKI_GFP+_4weeks_3/strainedCounts")
AKI_GFPneg_4weeks <- Read10X("./AKI_GFP-_4weeks/strainedCounts")

Ctrl_4weeks_1 <- CreateSeuratObject(counts = Ctrl_4weeks_1, min.cells = 0, min.features = 0,project="Ctrl_4weeks_1")
Ctrl_4weeks_1
head(Ctrl_4weeks_1@meta.data)

AKI_GFPpos_4weeks_1 <- CreateSeuratObject(counts = AKI_GFPpos_4weeks_1, min.cells = 0, min.features = 0,project="AKI_GFPpos_4weeks_1")
AKI_GFPpos_4weeks_1
AKI_GFPpos_4weeks_1 <- RenameCells(AKI_GFPpos_4weeks_1, new.names = c(paste0(substr(row.names(AKI_GFPpos_4weeks_1@meta.data),1,17),rep("2",length(row.names(AKI_GFPpos_4weeks_1@meta.data))))))
head(AKI_GFPpos_4weeks_1@meta.data)

Ctrl_4weeks_2 <- CreateSeuratObject(counts = Ctrl_4weeks_2, min.cells = 0, min.features = 0,project="Ctrl_4weeks_2")
Ctrl_4weeks_2
Ctrl_4weeks_2 <- RenameCells(Ctrl_4weeks_2, new.names = c(paste0(substr(row.names(Ctrl_4weeks_2@meta.data),1,17),rep("3",length(row.names(Ctrl_4weeks_2@meta.data))))))
head(Ctrl_4weeks_2@meta.data)

AKI_GFPpos_4weeks_2 <- CreateSeuratObject(counts = AKI_GFPpos_4weeks_2, min.cells = 0, min.features = 0,project="AKI_GFPpos_4weeks_2")
AKI_GFPpos_4weeks_2
AKI_GFPpos_4weeks_2 <- RenameCells(AKI_GFPpos_4weeks_2, new.names = c(paste0(substr(row.names(AKI_GFPpos_4weeks_2@meta.data),1,17),rep("4",length(row.names(AKI_GFPpos_4weeks_2@meta.data))))))
head(AKI_GFPpos_4weeks_2@meta.data)

AKI_GFPpos_6months_1 <- CreateSeuratObject(counts = AKI_GFPpos_6months_1, min.cells = 0, min.features = 0,project="AKI_GFPpos_6months_1")
AKI_GFPpos_6months_1
AKI_GFPpos_6months_1 <- RenameCells(AKI_GFPpos_6months_1, new.names = c(paste0(substr(row.names(AKI_GFPpos_6months_1@meta.data),1,17),rep("5",length(row.names(AKI_GFPpos_6months_1@meta.data))))))
head(AKI_GFPpos_6months_1@meta.data)

Ctrl_6months <- CreateSeuratObject(counts = Ctrl_6months, min.cells = 0, min.features = 0,project="Ctrl_6months")
Ctrl_6months
Ctrl_6months <- RenameCells(Ctrl_6months, new.names = c(paste0(substr(row.names(Ctrl_6months@meta.data),1,17),rep("6",length(row.names(Ctrl_6months@meta.data))))))
head(Ctrl_6months@meta.data)

AKI_GFPpos_6months_2 <- CreateSeuratObject(counts = AKI_GFPpos_6months_2, min.cells = 0, min.features = 0,project="AKI_GFPpos_6months_2")
AKI_GFPpos_6months_2
AKI_GFPpos_6months_2 <- RenameCells(AKI_GFPpos_6months_2, new.names = c(paste0(substr(row.names(AKI_GFPpos_6months_2@meta.data),1,17),rep("7",length(row.names(AKI_GFPpos_6months_2@meta.data))))))
head(AKI_GFPpos_6months_2@meta.data)

AKI_GFPpos_6months_3 <- CreateSeuratObject(counts = AKI_GFPpos_6months_3, min.cells = 0, min.features = 0,project="AKI_GFPpos_6months_3")
AKI_GFPpos_6months_3
AKI_GFPpos_6months_3 <- RenameCells(AKI_GFPpos_6months_3, new.names = c(paste0(substr(row.names(AKI_GFPpos_6months_3@meta.data),1,17),rep("8",length(row.names(AKI_GFPpos_6months_3@meta.data))))))
head(AKI_GFPpos_6months_3@meta.data)

AKI_GFPneg_6months <- CreateSeuratObject(counts = AKI_GFPneg_6months, min.cells = 0, min.features = 0,project="AKI_GFPneg_6months")
AKI_GFPneg_6months
AKI_GFPneg_6months <- RenameCells(AKI_GFPneg_6months, new.names = c(paste0(substr(row.names(AKI_GFPneg_6months@meta.data),1,17),rep("9",length(row.names(AKI_GFPneg_6months@meta.data))))))
head(AKI_GFPneg_6months@meta.data)

AKI_GFPpos_6months_4 <- CreateSeuratObject(counts = AKI_GFPpos_6months_4, min.cells = 0, min.features = 0,project="AKI_GFPpos_6months_4")
AKI_GFPpos_6months_4
AKI_GFPpos_6months_4 <- RenameCells(AKI_GFPpos_6months_4, new.names = c(paste0(substr(row.names(AKI_GFPpos_6months_4@meta.data),1,17),rep("10",length(row.names(AKI_GFPpos_6months_4@meta.data))))))
head(AKI_GFPpos_6months_4@meta.data)

AKI_GFPpos_4weeks_3 <- CreateSeuratObject(counts = AKI_GFPpos_4weeks_3, min.cells = 0, min.features = 0,project="AKI_GFPpos_4weeks_3")
Ki67_91P_4wAKI_GFPpos_4weeks_3
AKI_GFPpos_4weeks_3 <- RenameCells(AKI_GFPpos_4weeks_3, new.names = c(paste0(substr(row.names(AKI_GFPpos_4weeks_3@meta.data),1,17),rep("11",length(row.names(AKI_GFPpos_4weeks_3@meta.data))))))
head(AKI_GFPpos_4weeks_3@meta.data)

AKI_GFPneg_4weeks <- CreateSeuratObject(counts = AKI_GFPneg_4weeks, min.cells = 0, min.features = 0,project="AKI_GFPneg_4weeks")
AKI_GFPneg_4weeks
AKI_GFPneg_4weeks <- RenameCells(AKI_GFPneg_4weeks, new.names = c(paste0(substr(row.names(AKI_GFPneg_4weeks@meta.data),1,17),rep("12",length(row.names(AKI_GFPneg_4weeks@meta.data))))))
head(AKI_GFPneg_4weeks@meta.data)

Ki67 <- merge(x = Ctrl_4weeks_1, y = c(AKI_GFPpos_4weeks_1,Ctrl_4weeks_2,AKI_GFPpos_4weeks_2,AKI_GFPpos_6months_1,Ctrl_6months,
                                       AKI_GFPpos_6months_2,AKI_GFPpos_6months_3,AKI_GFPneg_6months,AKI_GFPpos_6months_4,AKI_GFPpos_4weeks_3,AKI_GFPneg_4weeks))

rm(Ctrl_4weeks_1,AKI_GFPpos_4weeks_1,Ctrl_4weeks_2,AKI_GFPpos_4weeks_2,AKI_GFPpos_6months_1,Ctrl_6months,
   AKI_GFPpos_6months_2,AKI_GFPpos_6months_3,AKI_GFPneg_6months,AKI_GFPpos_6months_4,AKI_GFPpos_4weeks_3,AKI_GFPneg_4weeks)

head(Ki67@meta.data)
tail(Ki67@meta.data)

GEX <- Ki67@assays$RNA@counts
head(GEX) 
rm(Ki67)

saveRDS(GEX,"./Ki67_aggr_GEX_SoupX.Rds")


