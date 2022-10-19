#Title: Differentially expressed genes between AKI_GFP- and control cells per cell type

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)
library(tibble)
set.seed(1234)

#Read in the final R object (contains RNA, chromVAR and peaks assay (peaks called using MACS2))
Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

Idents(Ki67) <- "Celltype"

DefaultAssay(Ki67) <- "RNA"

Celltype <- levels(Idents(Ki67))
Celltype <- Celltype[2:length(Celltype)] #exclude podocytes (too few cells for meaningful comparison)

for(i in 1:length(Celltype)){
  SUB <- subset(Ki67, idents=Celltype[[i]])
  Idents(SUB) <- "Group"
  table(Idents(SUB))
  DefaultAssay(SUB) <- "RNA"
  SUB <- NormalizeData(SUB)
  DEG_SUB <- FindMarkers(SUB, ident.1 = "AKI_GFP-", ident.2="Control", min.pct = 0.1, logfc.threshold = 0,only.pos = FALSE)
  DEG_SUB$Gene <- row.names(DEG_SUB)
  DEG_SUB <- DEG_SUB[which(DEG_SUB$p_val_adj < 0.01),]
  write.csv(DEG_SUB,paste0("./GFPneg_vs_Ctrl/DEG/2022_3_31_DEG_GFPnegvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),row.names = F)
  rm(DEG_SUB,SUB)
}

