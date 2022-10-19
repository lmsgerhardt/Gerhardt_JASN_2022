#Title: Identify differentially accessible regions 

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
set.seed(1234)

Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

DefaultAssay(Ki67) <- 'peaks'

da_peaks <- FindAllMarkers(object = Ki67,min.pct = 0.05,test.use = 'LR',latent.vars = 'atac_peak_region_fragments',logfc.threshold = 0.25)
da_peaks$DAR <- row.names(da_peaks)

# Identify closest feature
cf <- ClosestFeature(Ki67,rownames(da_peaks))
da_peaks <- cbind(da_peaks, gene=cf$gene_name, gene_id=cf$gene_id,distance=cf$distance)

da_peaks <- da_peaks[which(da_peaks$p_val_adj<0.01),]
write.csv(da_peaks,"./Ki67_DAR_MACS2_logFC0.25_minpct0.05_padj0.01.csv",row.names = T)