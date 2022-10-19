#Title: Making of Seurat object and quality control / filtering

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
set.seed(1234)

# Create a Seurat object based on the gene expression data, and then add in the ATAC-seq data as a second assay. 
# load the 10x output file (filtered_feature_bc_matrix.h5), which contains RNA and ATAC data.
inputdata.10x <- Read10X_h5("./Cellranger_outs/Aggregated/outs/Aggregated_filtered_feature_bc_matrix.h5")

#Read in the aggr.csv containing the names of the aggregated samples
aggcsv <- read.csv("./Cellranger_outs/Aggregated/outs/aggr.csv", header = TRUE, row.names = 1)

#Read in SoupX corrected counts
SoupX_counts <- readRDS("./SoupX/Ki67_aggr_GEX_SoupX.Rds")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
Ki67 <- CreateSeuratObject(counts = rna_counts)
range(Ki67@assays$RNA@counts) 
head(Ki67@meta.data)

#Order SoupX_counts according to order of Ki67 object
SoupX_counts_ordered <- SoupX_counts[,row.names(Ki67@meta.data)]
SoupX_counts_ordered[1:5,1:5]
unique(row.names(SoupX_counts_ordered)==row.names(Ki67@assays$RNA@counts)) 

#Replace counts with SoupX corrected counts (ambient RNA removed)
Ki67@assays$RNA@counts <- SoupX_counts_ordered
Ki67@assays$RNA@data <- SoupX_counts_ordered
range(Ki67@assays$RNA@counts) 
Ki67 <- NormalizeData(Ki67)
range(Ki67@assays$RNA@data) 

# Now add in the ATAC-seq data 
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
atac_counts[1:5,1:5]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) 
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "Ki67_analysis/Cellranger_outs/Aggregated/outs/Aggregated_atac_fragments.tsv.gz" #replace with full path to fragment file
chrom_assay <- CreateChromatinAssay(counts = atac_counts,sep = c(":", "-"),genome = 'mm10',fragments = frag.file,min.cells = 3,annotation = annotations)

Ki67[["ATAC"]] <- chrom_assay
Ki67

saveRDS(Ki67,"Ki67raw.Rds")

################################### ADD METADATA ########################################
# Add metadata to the Seurat object
# Some parts adapted from Parker Wilson's github: p4rkerw/Muto_Wilson_NComm_2020 (https://github.com/p4rkerw/Muto_Wilson_NComm_2020/blob/master/snATAC_prep/signac_atac_process.R [last accessed 10/17/22])
# Add orig. ident.
samples <- sapply(strsplit(rownames(Ki67@meta.data), split="-"), "[[", 2) 
label <- seq(length(rownames(aggcsv))) 
orig.ident <- rownames(aggcsv)
orig.ident

sampleID <- plyr::mapvalues(samples, from = label, to = orig.ident)
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(orig.ident=sampleID, row.names=rownames(Ki67@meta.data)))
head(Ki67@meta.data)

# Add metadata on group and group_time
ORIGIDENT <- Ki67@meta.data$orig.ident
ORIGIDENT[which(ORIGIDENT %in% c("Ctrl_4weeks_1","Ctrl_4weeks_2","Ctrl_6months"))] <- "Control"
ORIGIDENT[which(ORIGIDENT %in% c("AKI_GFP+_4weeks_1","AKI_GFP+_4weeks_2","AKI_GFP+_4weeks_3","AKI_GFP+_6months_1","AKI_GFP+_6months_2","AKI_GFP+_6months_3","AKI_GFP+_6months_4"))] <- "AKI_GFP+"
ORIGIDENT[which(ORIGIDENT %in% c("AKI_GFP-_4weeks","AKI_GFP-_6months"))] <- "AKI_GFP-"
names(ORIGIDENT) <- rownames(Ki67@meta.data)
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(Group=ORIGIDENT, row.names=rownames(Ki67@meta.data)))
head(Ki67@meta.data)

ORIGIDENT <- Ki67@meta.data$orig.ident
ORIGIDENT[which(ORIGIDENT %in% c("Ctrl_4weeks_1","Ctrl_4weeks_2"))] <- "Ctrl_4weeks"
ORIGIDENT[which(ORIGIDENT %in% c("Ctrl_6months"))] <- "Ctrl_6months"
ORIGIDENT[which(ORIGIDENT %in% c("AKI_GFP+_4weeks_1","AKI_GFP+_4weeks_2","AKI_GFP+_4weeks_3"))] <- "AKI_GFP+_4weeks"
ORIGIDENT[which(ORIGIDENT %in% c("AKI_GFP+_6months_1","AKI_GFP+_6months_2","AKI_GFP+_6months_3","AKI_GFP+_6months_4"))] <- "AKI_GFP+_6months"
names(ORIGIDENT) <- rownames(Ki67@meta.data)
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(Group_time=ORIGIDENT, row.names=names(ORIGIDENT)))

# Add metadata on mitochondrial and ribosomal count
Ki67[["percent.mt"]] <- PercentageFeatureSet(Ki67, pattern = "^mt-")
range(Ki67@meta.data$percent.mt) 
C<-GetAssayData(object = Ki67, slot = "counts")
rb.genes <- rownames(Ki67)[grep("^Rp[sl]",rownames(Ki67))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
Ki67 <- AddMetaData(Ki67, percent.ribo, col.name = "percent.ribo")
range(Ki67@meta.data$percent.ribo)
head(Ki67@meta.data)

# Calculate TSSEnrichment and nucleosome signal
DefaultAssay(Ki67) <- "ATAC"
Ki67 <- TSSEnrichment(Ki67, fast = TRUE)
Ki67 <- NucleosomeSignal(object = Ki67)

#Add information on peak_region_fragments and TSS fragments. This information is not provided by cellranger aggr, but only by cellranger count
Ctrl_4weeks_1_1 <- read.csv("./Cellranger_outs/Ctrl_4weeks_1/outs/per_barcode_metrics.csv",header=T)
#subset only the true cells
Ctrl_4weeks_1_1 <- Ctrl_4weeks_1_1[which(Ctrl_4weeks_1_1$is_cell ==1),]
#Check if the barcode column represents the barcodes actually used in Ki67 object after subsetting only the true cells
unique(Ctrl_4weeks_1_1$barcode %in%row.names(Ki67@meta.data[which(Ki67@meta.data$orig.ident == "Ctrl_4weeks_1"),]))

ORIGIDENT[which(ORIGIDENT %in% c("Ctrl_4weeks_1","Ctrl_4weeks_2"))] <- "Ctrl_4weeks"
ORIGIDENT[which(ORIGIDENT %in% c("Ctrl_6months"))] <- "Ctrl_6months"
ORIGIDENT[which(ORIGIDENT %in% c("AKI_GFP+_4weeks_1","AKI_GFP+_4weeks_2","AKI_GFP+_4weeks_3"))] <- "AKI_GFP+_4weeks"
ORIGIDENT[which(ORIGIDENT %in% c("AKI_GFP+_6months_1","AKI_GFP+_6months_2","AKI_GFP+_6months_3","AKI_GFP+_6months_4"))] <- "AKI_GFP+_6months"

#Read barcode information in for all samples and subset only true cells
AKI_GFPpos_4weeks_1_2 <- read.csv("./Cellranger_outs/AKI_GFPpos_4weeks_1/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_4weeks_1_2 <- AKI_GFPpos_4weeks_1_2[which(AKI_GFPpos_4weeks_1_2$is_cell ==1),]
Ctrl_4weeks_2_3 <- read.csv("./Cellranger_outs/Ctrl_4weeks_2/outs/per_barcode_metrics.csv",header=T)
Ctrl_4weeks_2_3 <- Ctrl_4weeks_2_3[which(Ctrl_4weeks_2_3$is_cell ==1),]
AKI_GFPpos_4weeks_2_4 <- read.csv("./Cellranger_outs/AKI_GFPpos_4weeks_2/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_4weeks_2_4 <- AKI_GFPpos_4weeks_2_4[which(AKI_GFPpos_4weeks_2_4$is_cell ==1),]
AKI_GFPpos_6months_1_5 <- read.csv("./Cellranger_outs/AKI_GFPpos_6months_1/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_6months_1_5 <- AKI_GFPpos_6months_1_5[which(AKI_GFPpos_6months_1_5$is_cell ==1),]
Ctrl_6months_6 <- read.csv("./Cellranger_outs/Ctrl_6months_6/outs/per_barcode_metrics.csv",header=T)
Ctrl_6months_6 <- Ctrl_6months_6[which(Ctrl_6months_6$is_cell ==1),]
AKI_GFPpos_6months_2_7 <- read.csv("./Cellranger_outs/AKI_GFPpos_6months_2/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_6months_2_7 <- AKI_GFPpos_6months_2_7[which(AKI_GFPpos_6months_2_7$is_cell ==1),]
AKI_GFPpos_6months_3_8 <- read.csv(".Cellranger_outs/AKI_GFPpos_6months_3/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_6months_3_8 <- AKI_GFPpos_6months_3_8[which(AKI_GFPpos_6months_3_8$is_cell ==1),]
AKI_GFPneg_6months_9 <- read.csv("./Cellranger_outs/AKI_GFPneg_6months/outs/per_barcode_metrics.csv",header=T)
AKI_GFPneg_6months_9 <- AKI_GFPneg_6months_9[which(AKI_GFPneg_6months_9$is_cell ==1),]
AKI_GFPpos_6months_4_10 <- read.csv("./Cellranger_outs/AKI_GFPpos_6months_4/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_6months_4_10 <- AKI_GFPpos_6months_4_10[which(AKI_GFPpos_6months_4_10$is_cell ==1),]
AKI_GFPpos_4weeks_3_11 <- read.csv("./Cellranger_outs/AKI_GFPpos_4weeks_3/outs/per_barcode_metrics.csv",header=T)
AKI_GFPpos_4weeks_3_11 <- AKI_GFPpos_4weeks_3_11[which(AKI_GFPpos_4weeks_3_11$is_cell ==1),]
AKI_GFPneg_4weeks_12 <- read.csv("./Cellranger_outs/AKI_GFPneg_4weeks_12/outs/per_barcode_metrics.csv",header=T)
AKI_GFPneg_4weeks_12 <- AKI_GFPneg_4weeks_12[which(AKI_GFPneg_4weeks_12$is_cell ==1),]

# Format barcodes to fit naming in the aggregated seurat object (only barcode column)
AKI_GFPpos_4weeks_1_2$barcode <- paste(substr(AKI_GFPpos_4weeks_1_2$barcode,1,17),rep(2,length(AKI_GFPpos_4weeks_1_2$barcode)),sep="")
Ctrl_4weeks_2_3$barcode <- paste(substr(Ctrl_4weeks_2_3$barcode,1,17),rep(3,length(Ctrl_4weeks_2_3$barcode)),sep="")
AKI_GFPpos_4weeks_2_4$barcode <- paste(substr(AKI_GFPpos_4weeks_2_4$barcode,1,17),rep(4,length(AKI_GFPpos_4weeks_2_4$barcode)),sep="")
AKI_GFPpos_6months_1_5$barcode <- paste(substr(AKI_GFPpos_6months_1_5$barcode,1,17),rep(5,length(AKI_GFPpos_6months_1_5$barcode)),sep="")
Ctrl_6months_6$barcode <- paste(substr(Ctrl_6months_6$barcode,1,17),rep(6,length(Ctrl_6months_6$barcode)),sep="")
AKI_GFPpos_6months_2_7$barcode <- paste(substr(AKI_GFPpos_6months_2_7$barcode,1,17),rep(7,length(AKI_GFPpos_6months_2_7$barcode)),sep="")
AKI_GFPpos_6months_3_8$barcode <- paste(substr(AKI_GFPpos_6months_3_8$barcode,1,17),rep(8,length(AKI_GFPpos_6months_3_8$barcode)),sep="")
AKI_GFPneg_6months_9$barcode <- paste(substr(AKI_GFPneg_6months_9$barcode,1,17),rep(9,length(AKI_GFPneg_6months_9$barcode)),sep="")
AKI_GFPpos_6months_4_10$barcode <- paste(substr(AKI_GFPpos_6months_4_10$barcode,1,17),rep(10,length(AKI_GFPpos_6months_4_10$barcode)),sep="")
AKI_GFPpos_4weeks_3_11$barcode <- paste(substr(AKI_GFPpos_4weeks_3_11$barcode,1,17),rep(11,length(AKI_GFPpos_4weeks_3_11$barcode)),sep="")
AKI_GFPneg_4weeks_12$barcode <- paste(substr(AKI_GFPneg_4weeks_12$barcode,1,17),rep(12,length(AKI_GFPneg_4weeks_12$barcode)),sep="")

myDF <- rbind(Ctrl_4weeks_1_1,AKI_GFPpos_4weeks_1_2,Ctrl_4weeks_2_3,AKI_GFPpos_4weeks_2_4,AKI_GFPpos_6months_1_5,Ctrl_6months_6,AKI_GFPpos_6months_2_7,
              AKI_GFPpos_6months_3_8,AKI_GFPneg_6months_9,AKI_GFPpos_6months_4_10,AKI_GFPpos_4weeks_3_11,AKI_GFPneg_4weeks_12)

# Add metadata on peak_region_fragments and TSS fragments
atac_TSS_fragments <- myDF$atac_TSS_fragments
names(atac_TSS_fragments) <- myDF$barcode
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(atac_TSS_fragments=atac_TSS_fragments, 
                                                     row.names=names(atac_TSS_fragments)))
atac_peak_region_fragments<- myDF$atac_peak_region_fragments
names(atac_peak_region_fragments) <- myDF$barcode
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(atac_peak_region_fragments=atac_peak_region_fragments,
                                                     row.names=names(atac_peak_region_fragments)))
atac_mitochondrial_reads<- myDF$atac_mitochondrial_reads
names(atac_mitochondrial_reads) <- myDF$barcode
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(atac_mitochondrial_reads=atac_mitochondrial_reads,
                                                     row.names=names(atac_mitochondrial_reads)))
Ki67@meta.data$atac_pct_fragments_in_peaks <- (Ki67@meta.data$atac_peak_region_fragments/Ki67@meta.data$nCount_ATAC)*100

pdf("Ki67_QCplots.pdf",width=45)
Idents(Ki67) <- "orig.ident"
levels(Ki67) <- c("Ctrl_4weeks_1","Ctrl_4weeks_2","Ctrl_6months", "AKI_GFP+_4weeks_1", "AKI_GFP+_4weeks_2", "AKI_GFP+_4weeks_3", "AKI_GFP+_6months_1", "AKI_GFP+_6months_2", "AKI_GFP+_6months_3", "AKI_GFP+_6months_4","AKI_GFP-_4weeks" ,"AKI_GFP-_6months")
VlnPlot(Ki67, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_ATAC", "atac_peak_region_fragments","TSS.enrichment",'nucleosome_signal'), ncol=8,
        log = TRUE, pt.size = 0) + NoLegend()
dev.off()

saveRDS(Ki67,"./Ki67raw.Rds")

################################### Remove doublets and low quality nuclei ########################################
# Filter nFeature_RNA >= 200 as was done for the identification of heterotypic doublets using DoubletFinder
Ki67 <- subset(x = Ki67, nFeature_RNA >= 200) 
Ki67
range(Ki67@meta.data$nFeature_RNA) # somehow 3 nuclei with nFeature_RNA =< 200 remain --> repeated filtering removes these 3 nuclei
Ki67 <- subset(x = Ki67, nFeature_RNA >= 200) 
Ki67
range(Ki67@meta.data$nFeature_RNA)

#Add heterotypic doubletinformation and filter them
doublets <- read.csv("./Doubletfinder/Ki67_singlets_doublets_matrix.csv",header=T,stringsAsFactors = F)
length(unique(doublets$V1))

#Check that all nuclei in the Ki67 object are in the doubletinformation df
table(row.names(Ki67@meta.data)%in%doublets$V1) 

#Add singlet/doublet metadata
Ki67 <- AddMetaData(object=Ki67, metadata=data.frame(Singlet_Doublet=doublets$V2,
                                                     row.names=doublets$V1))
table(Ki67@meta.data$Singlet_Doublet)

Idents(Ki67) <- "Singlet_Doublet"
Ki67 <- subset(Ki67,idents="Singlet")
Ki67 

#Filter low quality nuclei
Ki67_filtered <- subset(
  x = Ki67,
  subset = 
    nFeature_RNA > 350&
    nFeature_RNA < 3500 &
    nCount_RNA < 8000 &
    nCount_RNA > 500 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    TSS.enrichment > 2&
    nucleosome_signal < 3 &
    percent.mt < 1 &
    percent.ribo < 3
)
Ki67_filtered

saveRDS(Ki67_filtered,"Ki67_filtered.Rds")

sessionInfo()