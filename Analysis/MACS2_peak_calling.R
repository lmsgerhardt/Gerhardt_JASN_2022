#Title: peak calling with MACS2 

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(Matrix)
set.seed(1234)

Ki67 <- readRDS("./Ki67_FINAL.Rds")
Ki67

Idents(Ki67) <- "Group_time_Celltype"

table(Idents(Ki67))

DefaultAssay(Ki67) <- "ATAC"

peaks <- CallPeaks(
  object = Ki67,
  group.by = "Group_time_Celltype",
  macs2.path = "/lgerhard/miniconda3_2/envs/Seurat_Signac_MACS2/bin/macs2"
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks,"2022_04_01_peaks.Rds")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(Ki67),
  features = peaks,
  cells = colnames(Ki67)
)

saveRDS(macs2_counts,"2022_04_01_macs2_counts.Rds")

frag.file <- "Ki67_analysis/Cellranger_outs/Aggregated/outs/Aggregated_atac_fragments.tsv.gz" #replace with full path to fragment file
annotations <- readRDS("./annotations.Rds")

# create a new assay using the MACS2 peak set and add it to the Seurat object
Ki67[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 3,
  annotation = annotations
)
DefaultAssay(Ki67) <- "ATAC"
Ki67
saveRDS(Ki67,"./Ki67_FINAL_macs2_ATAC.Rds")

#Save Copy of Rds without ATAC assay to reduce size
DefaultAssay(Ki67) <- "peaks"
Ki67[["ATAC"]] <- NULL
saveRDS(Ki67,"./Ki67_FINAL_macs2.Rds")

sessionInfo()