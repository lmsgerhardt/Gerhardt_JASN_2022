#Title: Differentially accessible regions between AKI_GFP+ PTC and LOH subclusters as identified in Figure 4A

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
DefaultAssay(PTC) <- 'peaks'

#Find differentially accessible peaks
da_peaks <- FindAllMarkers(object = PTC,min.pct = 0.05,test.use = 'LR',latent.vars = 'atac_peak_region_fragments',logfc.threshold = 0.25)
colnames(da_peaks)[which(colnames(da_peaks)=="gene")] <- "DAR"

#Identify closest features
cf <- ClosestFeature(PTC,rownames(da_peaks))
da_peaks <- cbind(da_peaks, Closest_gene=cf$gene_name, gene_id=cf$gene_id,distance=cf$distance)
da_peaks <- da_peaks[which(da_peaks$p_val_adj<0.01),]
da_peaks <- da_peaks[,c("cluster","DAR","Closest_gene","gene_id","distance","pct.1","pct.2","avg_log2FC","p_val","p_val_adj")] #Reorder
write.csv(da_peaks,"./Subclustering_Analysis/Supplemental_Table_5.csv",row.names = F)

