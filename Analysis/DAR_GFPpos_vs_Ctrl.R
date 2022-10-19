#Title: Differentially accessible regions between AKI_GFP+ and control cells per cell type

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
set.seed(1234)

#Read in the final R object (contains RNA, chromVAR and peaks assay (peaks called using MACS2))
Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

Idents(Ki67) <- "Celltype"

DefaultAssay(Ki67) <- "peaks"

clusters <- levels(Idents(Ki67))
clusters <- clusters[2:length(clusters)] #exclude podocytes (too few cells for meaningful comparison)

### 4 weeks and 6 months post AKI samples pooled
for(i in 1:length(clusters)){
  SUB <- subset(Ki67, idents=clusters[[i]])
  Idents(SUB) <- "Group"
  table(Idents(SUB))
  DefaultAssay(SUB) <- "peaks"
  DAR_SUB <- FindMarkers(SUB, ident.1 = "AKI_GFP+", ident.2="Control", 
                         min.pct = 0.05, logfc.threshold = 0,only.pos = FALSE,
                         test.use = 'LR',latent.vars = 'atac_peak_region_fragments')
  DAR_SUB$DAR <- row.names(DAR_SUB)
  cf <- ClosestFeature(SUB,rownames(DAR_SUB))
  DAR_SUB <- cbind(DAR_SUB, gene=cf$gene_name, gene_id=cf$gene_id,distance=cf$distance)
  DAR_SUB <- DAR_SUB[which(DAR_SUB$p_val_adj < 0.01),]
  write.csv(DAR_SUB,paste0("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_cluster",clusters[[i]],"_logFC0_minpct0.05_padj0.01.csv"),row.names = F)
  rm(DAR_SUB,SUB)
}

### only 4 weeks post AKI samples
for(i in 1:length(clusters)){
  SUB <- subset(Ki67, idents=clusters[[i]])
  Idents(SUB) <- "Group_time"
  table(Idents(SUB))
  DefaultAssay(SUB) <- "peaks"
  DAR_SUB <- FindMarkers(SUB, ident.1 = "AKI_GFP+_4weeks", ident.2="Ctrl_4weeks", 
                         min.pct = 0.05, logfc.threshold = 0,only.pos = FALSE,
                         test.use = 'LR',latent.vars = 'atac_peak_region_fragments')
  DAR_SUB$DAR <- row.names(DAR_SUB)
  cf <- ClosestFeature(SUB,rownames(DAR_SUB))
  DAR_SUB <- cbind(DAR_SUB, gene=cf$gene_name, gene_id=cf$gene_id,distance=cf$distance)
  DAR_SUB <- DAR_SUB[which(DAR_SUB$p_val_adj < 0.01),]
  write.csv(DAR_SUB,paste0("./GFPpos_vs_Ctrl/DAR/4weeks/2022_9_1_DAR_GFPpos4wvsCtrl4w_cluster",clusters[[i]],"_logFC0_minpct0.05_padj0.01.csv"),row.names = F)
  rm(DAR_SUB,SUB)
}

### only 6 months post AKI samples
for(i in 1:length(clusters)){
  SUB <- subset(Ki67, idents=clusters[[i]])
  Idents(SUB) <- "Group_time"
  table(Idents(SUB))
  DefaultAssay(SUB) <- "peaks"
  DAR_SUB <- FindMarkers(SUB, ident.1 = "AKI_GFP+_6months", ident.2="Ctrl_6months", 
                         min.pct = 0.05, logfc.threshold = 0,only.pos = FALSE,
                         test.use = 'LR',latent.vars = 'atac_peak_region_fragments')
  DAR_SUB$DAR <- row.names(DAR_SUB)
  cf <- ClosestFeature(SUB,rownames(DAR_SUB))
  DAR_SUB <- cbind(DAR_SUB, gene=cf$gene_name, gene_id=cf$gene_id,distance=cf$distance)
  write.csv(DAR_SUB,paste0("./GFPpos_vs_Ctrl/DAR/6months/2022_9_1_DAR_GFPpos6mvsCtrl6m_cluster",clusters[[i]],"_logFC0_minpct0.1_padj0.05.csv"),row.names = F)
  rm(DAR_SUB,SUB)
}