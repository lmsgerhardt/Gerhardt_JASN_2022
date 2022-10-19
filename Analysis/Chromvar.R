#Title: Motif analysis using ChromVAR on peaks called with MACS2

library(Seurat)
library(Signac)
library(dplyr)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)
library(GenomicRanges)
library(tibble)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

Ki67 <- readRDS("./Ki67_FINAL_macs2.Rds")
Ki67

Idents(Ki67) <- "Celltype"

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object. Use the peaks created by MACS2 for that.
DefaultAssay(Ki67) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)) #this is using the core vertebrate database
motif.matrix <- CreateMotifMatrix(features = granges(Ki67), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
Ki67 <- SetAssayData(Ki67, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# Run chromVAR
Ki67 <- RunChromVAR(
  object = Ki67,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(Ki67,"./Ki67_SObj_rna_peak_motif_assay.Rds")

# Identify differentially active motifs
DefaultAssay(Ki67) <- "chromvar"
Ki67

# Make a matrix with all the enriched motifs per celltype
da_motifs <- FindAllMarkers(
  object = Ki67,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.1,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0 #find all cluster-specific motifs
)

# Add motif name
DefaultAssay(Ki67) <- "peaks"
motif.names <- da_motifs$gene
da_motifs$gene <- ConvertMotifID(Ki67, id = motif.names)
da_motifs$motif <- motif.names 
colnames(da_motifs) <- c("motif.p_val","motif.avg_diff.","motif.pct.1","motif.pct.2","motif.p_val_adj","cluster","gene","motif")

#Filter to keep only the motifs with padj < 0.01
da_motifs <- da_motifs[which(da_motifs$motif.p_val_adj<0.01),]
write.csv(da_motifs,"./Ki67_chromvar_motif_logFC0_minpct0.1_padj0.01.csv")

sessionInfo()