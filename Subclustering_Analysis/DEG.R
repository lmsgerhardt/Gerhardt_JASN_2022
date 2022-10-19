#Title: Differentially expressed genes between AKI_GFP+ PTC and LOH subclusters as identified in Figure 4A

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

#Read in the final R object of AKI_GFP+ ("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM") subclustering analysis 
PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")
PTC

Idents(PTC) <- "Celltype"

DefaultAssay(PTC) <- "RNA"
PTC <- NormalizeData(PTC, normalization.method = "LogNormalize", scale.factor = 10000)
DEG <- FindAllMarkers(PTC,min.pct = 0.1, logfc.threshold = 0.25,only.pos = FALSE)
DEG <- DEG[DEG$p_val_adj < 0.01,]

DEG <- DEG[,c("cluster","gene","pct.1","pct.2","avg_log2FC","p_val","p_val_adj")] #Reorder

write.csv(da_peaks,"./Subclustering_Analysis/Supplemental_Table_4.csv",row.names = F)

