#Title: Link genes to peaks 

library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)

PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")
PTC

DefaultAssay(PTC) <- "peaks"
Idents(PTC) <- "Celltype"

# first compute the GC content for each peak
PTC <- RegionStats(PTC, genome = BSgenome.Mmusculus.UCSC.mm10)

# link peaks to genes
PTC <- LinkPeaks(
  object = PTC,
  peak.assay = "peaks",
  expression.assay = "RNA",
  distance = 100000,
  min.distance = NULL,
)

Links(PTC)

saveRDS(PTC,"./Subclustering_Analysis/PTC_FINAL_SUB.Rds")

# save matrix with link information
Links_Matrix <- as.data.frame(Links(PTC))
write.csv(Links_Matrix,"./Subclustering_Analysis/PTC_allLinks.csv")