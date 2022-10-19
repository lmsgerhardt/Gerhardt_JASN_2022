#Title: ChromVAR motif analysis 

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

PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")
PTC

Idents(PTC) <- "Celltype"

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object. Use the peaks created by MACS2 for that.
DefaultAssay(PTC) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)) #this is using the core vertebrate database
motif.matrix <- CreateMotifMatrix(features = granges(PTC), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
PTC <- SetAssayData(PTC, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# Run chromVAR
PTC <- RunChromVAR(
  object = PTC,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(PTC,"./Subclustering_Analysis/PTC_FINAL_SUB.Rds")

# Identify differentially active motifs
DefaultAssay(PTC) <- "chromvar"
PTC

#Make a matrix with all the enriched motifs per cluster
da_motifs <- FindAllMarkers(
  object = PTC,
  only.pos = FALSE,
  test.use = 'LR',
  min.pct = 0.1,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0 #find all cluster-specific motifs
)

#Add motif name
DefaultAssay(PTC) <- "peaks"
motif.names <- da_motifs$gene
da_motifs$gene <- ConvertMotifID(PTC, id = motif.names)
da_motifs$motif <- motif.names 
colnames(da_motifs) <- c("motif.p_val","motif.avg_diff.","motif.pct.1","motif.pct.2","motif.p_val_adj","cluster","gene","motif")

#Filter to keep only the motifs with padj < 0.01
da_motifs <- da_motifs[which(da_motifs$motif.p_val_adj<0.01),]
write.csv(da_motifs,"./Subclustering_Analysis/PTCSUB_chromvar_motif_logFC0_minpct0.1_padj0.01.csv")
