#Title: Differential gene expression analysis

library(Seurat)
library(Signac)
set.seed(1234)

Ki67 <- readRDS("./Ki67_FINAL.Rds")
Ki67

Idents(Ki67) <- "Celltype"

DefaultAssay(Ki67) <- "RNA"
Ki67 <- NormalizeData(Ki67, normalization.method = "LogNormalize", scale.factor = 10000)
DEG <- FindAllMarkers(Ki67,min.pct = 0.1, logfc.threshold = 0.25,only.pos = FALSE)
DEG <- DEG[DEG$p_val_adj < 0.01,]
DEG <- DEG[,c("cluster","gene","pct.1","pct.2","avg_log2FC","p_val","p_val_adj")] #reorder
write.csv(DEG,"Supplemental_Table_3.csv",row.names = F)