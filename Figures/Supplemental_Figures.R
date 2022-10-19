#Title: Supplemental Figures

library(Signac) 
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(ggrepel)
library(tibble)
library(dplyr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene # TxDb object is an object form to store transcript metdata
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(Matrix)
library(GenomicRanges)
set.seed(1234)

setwd()

#################################### Supplemental Figure 6 - QC Violin plots #################################### 
# Read in the unfiltered Seurat object
Ki67raw <- readRDS("./Ki67raw.Rds")
Ki67raw

pdf("Supplemental_Figure_6.pdf",width=45)
VlnPlot(Ki67raw, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_ATAC", "atac_peak_region_fragments","TSS.enrichment",'nucleosome_signal'), ncol=8,
        log = TRUE, pt.size = 0) + NoLegend()
dev.off()

#################################### Supplemental Figure 8 ####################################
### Supplemental Figure 8A - annotation of DARs in subclustering analysis (clusters as in Figure 4A)
#Read in the final R object of AKI_GFP+ ("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM") subclustering analysis 
PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")

# Read in DARs between subclusters 
DAR <- read.csv("./Subclustering_Analysis/Supplemental_Table_5.csv")
DAR <- unique(DAR$DAR)

# convert the DAR to GRanges objects to annotate
DAR.gr <- StringToGRanges(DAR, sep = c(":","-"))
head(DAR.gr)

# annotate the list of GRanges DAR for all clusters and make pie plot
all.peakAnno.1kb <- annotatePeak(DAR.gr, TxDb = txdb,level = "transcript",
                                 tssRegion = c(-1000, 1000), verbose = FALSE) 

all.peakAnno.1kb 

pdf("Supplemental_Figure_8A.pdf")
plotAnnoPie(all.peakAnno.1kb)
dev.off()

# Make data frame and save annotations
#Get dataframe
all_peakAnno_1kb_df <- as.data.frame(all.peakAnno.1kb)

#Convert entrez into gene ID using Biomart
# Get the entrez IDs
entrez <-all_peakAnno_1kb_df$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Mmusculus.v79,keys = entrez,columns = c("GENENAME"),keytype = "ENTREZID")
head(annotations_edb)
length(unique(annotations_edb$ENTREZID)) 
annotations_edb_nodup <- annotations_edb[!duplicated(annotations_edb$ENTREZID),]

# Change IDs to character type to merge
annotations_edb_nodup$ENTREZID <- as.character(annotations_edb_nodup$ENTREZID)
dim(annotations_edb_nodup)

# Write to file
all_peakAnno_1kb_df <- all_peakAnno_1kb_df %>% left_join(annotations_edb_nodup, by=c("geneId"="ENTREZID")) 

#Reorder and annotate dataframe
colnames(all_peakAnno_1kb_df)
all_peakAnno_1kb_df <- all_peakAnno_1kb_df[,c("cluster","DAR","seqnames","start","end","width","strand",
                                              "avg_log2FC","pct.1","pct.2","p_val_adj","p_val","annotation",
                                              "distanceToTSS","geneId","transcriptId","GENENAME","geneChr",
                                              "geneStart","geneEnd","geneLength","geneStrand","Closest_Gene","Closest_Gene_id","Closest_Gene_distance")]
colnames(all_peakAnno_1kb_df) <- c("cluster","DAR","seqnames","start","end","width","strand",
                                   "avg_log2FC","pct.1","pct.2","p_val_adj","p_val","annotation",
                                   "distanceToTSS","EntrezID_Chipseeker","TranscriptId_Chipseeker","Gene_Chipseeker","geneChr",
                                   "geneStart","geneEnd","geneLength","geneStrand","Closest_Gene","Closest_Gene_id","Closest_Gene_distance")

write.csv(all_peakAnno_1kb_df, "./Subclustering_Analysis/2022_5_31_peakAnno.1kb.allinfo.csv", row.names = F)

### Supplemental Figure 8B - Coverage plot for Ccl2 and links between peaks and gene
#Read in matrix with all links
Links <- read.csv("./Subclustering_Analysis/PTC_allLinks.csv")

#Read in DAR
DAR <- read.csv("./Subclustering_Analysis/Supplemental_Table_5.csv")

#Get DAR in FR-PTC to highlight them in coverage plot
DAR_FRPTC <- DAR[which(DAR$cluster=="FR-PTC"),]
DAR_FRPTC.gr <- StringToGRanges(DAR_FRPTC$DAR, sep = c(":","-"))

#Get coordinates for DAR related to Ccl2 gene or links 
Ccl2_links <- Links[which(Links$gene =="Ccl2"),"peak"]
Ccl2_links <-StringToGRanges(Ccl2_links, sep = c(":","-"))

#Get overlap with DAR_FRPTC
overlap_peaks <- findOverlaps(DAR_FRPTC.gr,Ccl2_links)

#Get actual coordinates
Overlaps_Ccl2.gr <- granges(DAR_FRPTC.gr)[queryHits(overlap_peaks)]

#Make plot
pdf("Supplemental_Figure_8B.pdf", width=5, height=5)
CoveragePlot(
  object = PTC,
  region = "Ccl2",
  features = "Ccl2",
  expression.assay = "RNA",
  idents = levels(PTC),
  extend.upstream = 60000,
  extend.downstream = 15000,
  region.highlight = Overlaps_Ccl2.gr
)
dev.off()

### Supplemental Figure 8C - Correlation of between cluster differences in TF GEX (log2FC) & motif activity (deviation score)
#Read in matrix with differentially expressed genes and differentially active motifs between clusters
df <- read.csv("./Subclustering_Analysis/Correlation_chromVar_GEX.csv")

# filter the df for pairs with significant changes in chromvar activity
df.f <- dplyr::filter(df, chrom_pval < 0.01, gene_pval < 0.01)
df.f <- na.omit(df.f)

# create a column for each unique gene-motif combo
df.f <- dplyr::mutate(df.f, gene_motif = paste0(gene,"_",motif))

# compute a pearson r2 and pval for each motif-gene combo across all celltypes
combos <- unique(df.f$gene_motif)
df.pearson <- lapply(combos, function(combo) {
  require(ggrepel)
  
  df <- dplyr::filter(df.f, gene_motif == combo)
  pearson <- tryCatch(cor.test(df$chromvar, df$rna, method="pearson", conf.level=0.95), error=function(e) NULL)
  df$corr <- tryCatch(signif(pearson$estimate,2), error=function(e) NULL)
  df$pval <- tryCatch(signif(pearson$p.value,2), error=function(e) NULL)
  df$num_celltypes <- length(unique(df$celltype))
  return(df)
})

# aggregate a final df with all the available pearson R2 coefficients and pvals
df.final <- bind_rows(df.pearson)
df.final <- dplyr::distinct(df.final, celltype, gene_motif, .keep_all = TRUE)
dim(df.final) 

#Take only one value per TF (because the corr coeff is identical across all celltypes)
df.final.unique <- df.final[!duplicated(df.final$motif),]
dim(df.final.unique)

#Add information of significance
df.final$Group <- df.final$pval
df.final$Group[which(df.final$Group < 0.05)] <- "Significant"
df.final$Group[which(!df.final$Group=="Significant")] <- "Not significant"
df.final$Group <- factor(df.final$Group,levels=rev(c("Significant","Not significant")))
# Change histogram plot line colors by groups
ggplot(df.final, aes(x=corr, color=Group,fill=Group)) +
  geom_histogram(alpha=0.5, position="identity",binwidth=0.05)+theme_classic()+
  ggtitle("Correlation of between cluster differences in TF gene \nexpression (log2FC) & in motif activity (deviation score)")+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+xlab("Pearson Coefficient")+ylab("Frequency")

#Add information of significance
df.final.unique$Group <- df.final.unique$pval
df.final.unique$Group[which(df.final.unique$Group < 0.05)] <- "Significant"
df.final.unique$Group[which(!df.final.unique$Group=="Significant")] <- "Not significant"
df.final.unique$Group <- factor(df.final.unique$Group,levels=rev(c("Significant","Not significant")))

# Make plot
pdf("Supplemental_Figure_8C.pdf", width=5)
ggplot(df.final.unique, aes(x=corr, color=Group,fill=Group)) +
  geom_histogram(alpha=0.5, position="identity",binwidth=0.05)+theme_classic()+
  ggtitle("Correlation of between cluster differences\nin TF gene expression (log2FC) & in motif activity\n(deviation score)")+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+xlab("Pearson Coefficient")+ylab("Frequency")
dev.off()

### Supplemental Figure 8D - Dot plot showing expression of BCAA and FAO genes
DefaultAssay(PTC) <- "RNA"
levels(PTC) <- rev(levels(PTC))

BCAA <- c("Bckdha","Bckdhb","Acadm","Mmut","Hibch","Ivd","Mccc1","Mccc2")
FAO <- c("Cd36","Cpt1a","Cpt2","Acox1","Acox2","Ppargc1a","Ppara")
BCAA_FAO <- c(BCAA, FAO)

pdf("Supplemental_Figure_8D.pdf", height=5)
DotPlot(PTC, features = BCAA_FAO, col.min = 0, cols = c("darkgrey", "black"),dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90,vjust=0.5),plot.title = element_text(hjust = 0.5))
dev.off()

#################################### Supplemental Figure 9 - TF-target gene interactions using DoRothEA ####################################
#Source: https://saezlab.github.io/dorothea/
#Download TF regulon dataset (dorothea_mm.rda) from DoRothEA (https://github.com/saezlab/dorothea/tree/master/data)
load("./Subclustering_Analysis/dorothea_mm.rda")
dorothea_mm[1:5,]

#Only keep confidence level A 
dorothea_mmA <- dorothea_mm[which(dorothea_mm$confidence=="A"),]
length(unique(dorothea_mmA$tf))
length(unique(dorothea_mmA$target))

# Get top10 differentially active transcription factors
#Read in list with differentially expressed genes
DEG <- read.csv("./Subclustering_Analysis/Supplemental_Table_4.csv",header=T)
#Consider only up regulated genes
DEG <- DEG[which(DEG$avg_log2FC>=0.25),]
colnames(DEG) <- c(paste0(rep("GEX.",5),colnames(DEG)[1:5]),"cluster","gene")
DEG$gene <- toupper(DEG$gene)

#Read in list with differentially active motifs
Motif <- read.csv("./Subclustering_Analysis/PTCSUB_chromvar_motif_logFC0_minpct0.1_padj0.01.csv",header=T)
Motif$gene <- toupper(Motif$gene)

#Retain only differentially active motifs that are also differentially expressed
OVERLAP <- data.frame("cluster"=c(),"gene"=c(),"motif"=c(),"motif.avg_diff."=c(), "motif.p_val"=c(),"motif.p_val_adj"=c(),"motif.pct.1"=c(),"motif.pct.2"=c(), 
                      "GEX.avg_log2FC"=c(),"GEX.p_val"=c(),"GEX.p_val_adj"=c(), "GEX.pct.1"=c(),"GEX.pct.2"=c())

Celltypes <- unique(Motif$cluster)

for(i in 1:length(Celltypes)){
  DEG_celltype <- DEG[which(DEG$cluster == Celltypes[[i]]),]
  Motif_celltype <- Motif[which(Motif$cluster== Celltypes[[i]]),]
  Motif_overlap <- Motif_celltype[which(Motif_celltype$gene%in%DEG_celltype$gene),]
  DEG_overlap <- DEG_celltype[which(DEG_celltype$gene%in%Motif_overlap$gene),]
  myDF <- merge(Motif_overlap,DEG_overlap,by=c("gene","cluster"))
  OVERLAP <- rbind(OVERLAP,myDF)
}

#Calculate a Motif_GEX_score and get top10 TF per cluster
OVERLAP$Motif_GEX_score <- OVERLAP$motif.avg_diff.* OVERLAP$GEX.avg_log2FC
top10 <- OVERLAP  %>% group_by(cluster) %>% top_n(n = 10, wt = Motif_GEX_score)

#Check intersect of top10 TF with DoRothEA data base
length(intersect(top10$gene,unique(toupper(dorothea_mmA$tf)))) 

# subset the dorothea_mmA dataframe for the identified TFs
dorothea_mmA <- dorothea_mmA[which(dorothea_mmA$tf%in%paste0(substr(top10$gene,1,1),tolower(substr(top10$gene,2,10)))),]

# subset the dorothea_mmA dataframe further to keep only the target genes that are differentially expressed between PT and LOH clusters as in Figure 4A 
DEG <- read.csv("./Subclustering_Analysis/Supplemental_Table_4.csv",header=T)
dorothea_mmA <- dorothea_mmA[which(dorothea_mmA$target%in%DEG$gene),]
dim(dorothea_mmA)

# subset the dorothea_mmA dataframe further to keep only the activating TF-target genes interactions 
dorothea_mmA <- dorothea_mmA[which(dorothea_mmA$mor=="1"),]

#Save list
write.csv(dorothea_mmA,"./Subclustering_Analysis/Dorothea_mm_A_mor1_top10TFs.csv",row.names = F)

#change column names for Gephi and save edges, i.e. connection between TFs ("Source") and target genes, to load into Gephi
colnames(dorothea_mmA)
colnames(dorothea_mmA) <- c("Source","confidence","Target","mor")
Edges <- dorothea_mmA[,c(1,3)]
write.csv(Edges,"./Subclustering_Analysis/Supplemental_Table_8.csv",row.names = F)

#save list of al genes in the data set with an "Id" and a "Label" column to load into Gephi
Nodes <- unique(c(dorothea_mmA$Source, dorothea_mmA$Target))
Nodes <- data.frame("Id" =Nodes, "Label"=Nodes)
write.csv(Nodes,"./Subclustering_Analysis/Dorothea_mm_A_mor1_top10TFs_Nodes.csv",row.names = F)

# Get z-scores for expression of all genes in the list across all PT (all segments pooled), injured PT (all segments pooled), FR-PTC, LOH-TL (LOH-TL of juxtamedullary and cortical nephrons pooled)
Idents(PTC) <- "Celltype"
levels(PTC) <- c("PTS1","PTS1-S2","PTS2","PTS3","Injured PTS1-S2","Injured PTS3","FR-PTC","LOH-TL-C","LOH-TL-JM")
        
#Pool clusters into 4 groups
new.cluster.ids <- c("PT-S1/S2/S3","PT-S1/S2/S3","PT-S1/S2/S3","PT-S1/S2/S3","Injured PT","Injured PT","FR-PTC","LOH-TL","LOH-TL")
names(new.cluster.ids) <- levels(PTC)
PTC <- RenameIdents(PTC, new.cluster.ids)
PTC@meta.data$GRN <- PTC@active.ident

DefaultAssay(PTC) <- "RNA"
PTC <- NormalizeData(PTC)

aver <- AverageExpression(PTC,features = unique(c(dorothea_mmA$Target, dorothea_mmA$Source)),assays = "RNA",slot = "data") %>% as.data.frame()

# Calculate z.score of gene expression across the clusters
Gene_names <- row.names(aver)
all_celltypes <- colnames(aver)
cluster.z_score <- aver

for(i in Gene_names){
  myD <- mean(as.numeric(aver[i,]))
  mySD <- sd(as.numeric(aver[i,]))
  for(x in all_celltypes)
  {
    cluster.z_score[i,x] <- (as.numeric(aver[i,x]) -myD)/mySD
  }
  rm(myD,mySD)
}

colnames(cluster.z_score) <- levels(PTC)
cluster.z_score$Id <- row.names(cluster.z_score)

PT <- cluster.z_score[,c(5,1)]
InjPT <- cluster.z_score[,c(5,2)]
FR <- cluster.z_score[,c(5,3)]
LOH <- cluster.z_score[,c(5,4)]

write.csv(PT,"./Subclustering_Analysis/PT_zscore_Avexp_dorotheammA_top10TFs.csv",row.names = F)
write.csv(InjPT,"./Subclustering_Analysis/InjPT_zscore_Avexp_dorotheammA_top10TFs.csv",row.names = F)
write.csv(FR,"./Subclustering_Analysis/FRPTC_zscore_Avexp_dorotheammA_top10TFs.csv",row.names = F)
write.csv(LOH,"./Subclustering_Analysis/LOH_zscore_Avexp_dorotheammA_top10TFs.csv",row.names = F)

#################################### Supplemental Figure 10 - Footprint analysis ####################################
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)

#Calculate footprint for all TFs with concordant GEX and Motif accessibility as identified in OVERLAP dataframe above
DefaultAssay(PTC) <- "peaks"
Idents(PTC) <- "Celltype"
levels(PTC)

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# add motif information
PTC <- AddMotifs(PTC, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pwm)

# gather the footprinting information for sets of motifs
PTC <- Footprint(object = PTC,motif.name = unique(OVERLAP$motif),genome = BSgenome.Mmusculus.UCSC.mm10,in.peaks = T)

MyM <- GetFootprintData(
  PTC,
  unique(OVERLAP$motif),
  assay = "peaks",
  group.by = "celltype",
  idents = NULL
)

write.csv(MyM,"./Subclustering_Analysis/Footprint_Matrix.csv")

# Plot footprint information for Hnf4a, Nfkb1 and Rela
pdf("Supplemental_Figure_10.pdf",height=10)
p1 <- PlotFootprint(PTC, features = c("MA0114.4","MA0105.4","MA0107.1"),label=FALSE)
p1 + patchwork::plot_layout(ncol = 1)
dev.off()

#################################### Supplemental Figure 11 - DAR between AKI_GFP+ and control cells per cell type ####################################
### Supplemental Figure 11A - Upset plot of DAR per cell type
library(ComplexHeatmap)

Celltype <- c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
              "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell")

Celltype <- Celltype[-1] #exclude podocytes (too few cells for meaningful comparison) 

#Make a list with all datasets of differentially accessible regions between AKI_GFP+ and control cells per cell types
DAR_list <- list()

for(i in 1:length(Celltype)){
  DAR <- read.csv(paste0("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  DAR_list[[Celltype[[i]]]] <- DAR[which(abs(DAR$avg_log2FC)>=0.1),"DAR"]
}

#Make the combination matrix (mode = distinct)
Mat <- make_comb_mat(DAR_list,mode = "distinct") 
Mat

#Make upset plots showing DAR sets >= 5 regions
Mat1 <- Mat[comb_size(Mat) >= 5]

col_size = comb_size(Mat1)
row_size = set_size(Mat1)
ht = UpSet(Mat1, set_order = Celltype,
           top_annotation = upset_top_annotation(Mat1, ylim = c(0, max(col_size)*1.1)),
           right_annotation = upset_right_annotation(Mat1, ylim = c(0, max(row_size)*1.1)))

pdf("Supplemental_Figure_11A.pdf",height = 6.5)
ht <- draw(ht)
col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("intersection_size", {
  grid.text(col_size[col_od], x = seq_along(col_size), y = unit(col_size[col_od], "native") + unit(2, "pt"), default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

decorate_annotation("set_size", {
  grid.text(row_size[row_od], 
            unit(row_size[row_od], "native") + unit(2, "mm"), 
            rev(seq_len(length(row_size))), 
            default.units = "native", just = "bottom", rot = -90,
            gp = gpar(fontsize = 8))
})
dev.off()

### Supplemental Figure 11B - Coverage plot for Keg1 across AKI_GFP+ and Control proximal tubule cells
#Read in the final R object (contains RNA, chromVAR and peaks assay (peaks called using MACS2))
Ki67 <- readRDS("./Ki67_SObj_rna_peak_motif_assay.Rds")
Ki67

# subset AKI_GFP+ and control proximal tubule cells for plotting, pool PTS1_1 and PTS1_2 for plotting
Idents(Ki67) <- "Group"
Ki67 <- subset(Ki67, idents=c("Control","AKI_GFP+"))
Idents(Ki67) <- "Celltype"
PT <- subset(Ki67, idents=c("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3"))
Celltype <- as.vector(PT@meta.data$Celltype)
Celltype[which(Celltype %in% c("PTS1_1","PTS1_2"))] <- "PTS1"
PT@meta.data$Celltype <- Celltype
PT@meta.data$Group_Celltype <- paste0(PT@meta.data$Group,sep="_",PT@meta.data$Celltype)

#read in DAR between AKI_GFP+ and control cells across proximal tubule clusters
DARPTS1_1 <- read.csv("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_clusterPTS1_1_logFC0_minpct0.05_padj0.01.csv")
DARPTS1_2 <- read.csv("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_clusterPTS1_2_logFC0_minpct0.05_padj0.01.csv")
DARPTS12 <- read.csv("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_clusterPTS1-S2_logFC0_minpct0.05_padj0.01.csv")
DARPTS2 <- read.csv("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_clusterPTS2_logFC0_minpct0.05_padj0.01.csv")
DARPTS3 <- read.csv("./GFPpos_vs_Ctrl/DAR/2022_3_31_DAR_GFPposvsCtrl_clusterPTS3_logFC0_minpct0.05_padj0.01.csv")
DAR <- rbind(DARPTS1_1,DARPTS1_2,DARPTS12,DARPTS2,DARPTS3)

DefaultAssay(PT) <- "peaks"
Idents(PT) <- "Group_Celltype"
levels(PT) <- c("AKI_GFP+_PTS1","Control_PTS1","AKI_GFP+_PTS1-S2","Control_PTS1-S2","AKI_GFP+_PTS2","Control_PTS2","AKI_GFP+_PTS3", "Control_PTS3")  

#highlight all DAR (min pct.1 0.05, avg_log2FC 0, padj < 0.01) for Keg1
Keg1 <- DAR[which(DAR$gene=="Keg1"&DAR$p_val_adj<0.01),]
Keg1 
Keg1 <- Keg1$DAR
ranges.show <- StringToGRanges(Keg1)
ranges.show$color <- "lightgrey"

pdf("Supplemental_Figure_11B.pdf", width=8)
CoveragePlot(
  object = PT,
  region.highlight = ranges.show,
  region = "Keg1",
  features = "Keg1",
  expression.assay = "RNA",
  idents = levels(PT),
  extend.downstream = 5000,
  extend.upstream = 5000)
dev.off()

#################################### Supplemental Figure 12 - DAR and DEG between AKI_GFP+ and control cells per cell type at 4 weeks and 6months post AKI ####################################
### Supplemental Figure 12A - DEG between AKI_GFP+ and control cells per cell type at 4 weeks and 6 months post AKI
Celltype <- c("Podo","PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM","CTAL","MTAL","DT", "DT-CNT","CNT",
              "PC_1","PC_2","IC-A","IC-B","Uro","Interst","Vasc","Macroph","Leuk other", "Bcell")

Celltype <- Celltype[-1] #exclude podocytes (too few cells for meaningful comparison) 

#Make a list with all datasets
#Read in all DEG IRIGFP vs Ctrl gene lists into one list
MyDF <- data.frame("Celltype"=c(),"Shared"=c(), "Concordant" = c(), "Notconcordant"=c(),
                   "4w_vs_4wCtrl"=c(),"6m_vs_6mCtrl"=c(), "Concordant_percent_of_6m"=c(),
                   "Notconcordant_percent_of_6m"=c(),"Concordant_percent_of_4w"=c(),
                   "Notconcordant_percent_of_4w"=c())

for(i in 1:length(Celltype)){
  AKI4w_vs_Ctrl4w <- read.csv(paste0("./GFPpos_vs_Ctrl/DEG/4weeks/2022_9_1_DEG_GFPpos4wvsCtrl4w_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  AKI4w_vs_Ctrl4wall <- AKI4w_vs_Ctrl4w[which(abs(AKI4w_vs_Ctrl4w$avg_log2FC)>=0.25),"Gene"]
  AKI4w_vs_Ctrl4wup <- AKI4w_vs_Ctrl4w[which(AKI4w_vs_Ctrl4w$avg_log2FC>=0.25),"Gene"]
  AKI4w_vs_Ctrl4wdown <- AKI4w_vs_Ctrl4w[which(AKI4w_vs_Ctrl4w$avg_log2FC<=-0.25),"Gene"]
  AKI6m_vs_Ctrl6m <- read.csv(paste0("./GFPpos_vs_Ctrl/DEG/6months/2022_9_1_DEG_GFPpos6mvsCtrl6m_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  AKI6m_vs_Ctrl6mall <- AKI6m_vs_Ctrl6m[which(abs(AKI6m_vs_Ctrl6m$avg_log2FC)>=0.25),"Gene"]
  AKI6m_vs_Ctrl6mup <- AKI6m_vs_Ctrl6m[which(AKI6m_vs_Ctrl6m$avg_log2FC>=0.25),"Gene"]
  AKI6m_vs_Ctrl6mdown <- AKI6m_vs_Ctrl6m[which(AKI6m_vs_Ctrl6m$avg_log2FC<=-0.25),"Gene"]
  #General overlap
  Shared <- length(intersect(AKI4w_vs_Ctrl4wall,AKI6m_vs_Ctrl6mall))
  Shared 
  Sharedup <- length(intersect(AKI4w_vs_Ctrl4wup,AKI6m_vs_Ctrl6mup))
  Shareddown <- length(intersect(AKI4w_vs_Ctrl4wdown,AKI6m_vs_Ctrl6mdown))
  Shareddown
  Concordant <- Sharedup + Shareddown
  Concordant
  Notconcordant <- Shared - Concordant
  Notconcordant
  
  myDF <- data.frame("Celltype"=Celltype[[i]],"Shared"=Shared, "Concordant" = Concordant, "Notconcordant"=Notconcordant,
                     "4w_vs_4wCtrl"=length(AKI4w_vs_Ctrl4wall),"6m_vs_6mCtrl"=length(AKI6m_vs_Ctrl6mall), "Concordant_percent_of_6m"=(Concordant/length(AKI6m_vs_Ctrl6mall))*100,
                     "Notconcordant_percent_of_6m"=(Notconcordant/length(AKI6m_vs_Ctrl6mall))*100,"Concordant_percent_of_4w"=(Concordant/length(AKI4w_vs_Ctrl4wall))*100,
                     "Notconcordant_percent_of_4w"=(Notconcordant/length(AKI4w_vs_Ctrl4wall))*100)
  MyDF <- rbind(MyDF, myDF)
}

#Add unique values
MyDF$Unique_6m <- MyDF$X6m_vs_6mCtrl - (MyDF$Concordant+MyDF$Notconcordant)
MyDF$Unique_4w <- MyDF$X4w_vs_4wCtrl - (MyDF$Concordant+MyDF$Notconcordant)

#Replace NA with 0
MyDF[is.na(MyDF)] <- 0

#Make DF to plot 
DF_toplot <- data.frame("Cell type"=rep(Celltype,3), "Group" = c(rep("Shared (concordant)",2*length(Celltype)),rep("Shared (nonconcordant)",2*length(Celltype)),
                                                                 rep("Unique",2*length(Celltype))) ,"Time point" = c(rep("4 weeks",length(Celltype)),rep("6 months",length(Celltype)),rep("4 weeks",length(Celltype)),
                                                                                                                     rep("6 months",length(Celltype)),rep("4 weeks",length(Celltype)),rep("6 months",length(Celltype))),
                        "Number of genes"=c(MyDF$Concordant,MyDF$Concordant,MyDF$Notconcordant,MyDF$Notconcordant,MyDF$Unique_4w,MyDF$Unique_6m))


#Make plot 
DF_toplot$Cell.type <- factor(DF_toplot$Cell.type,levels =Celltype)
DF_toplot$Group <- factor(DF_toplot$Group,levels=c("Unique","Shared (concordant)","Shared (nonconcordant)"))

pdf("Supplemental_Figure_12A.pdf",width=15, height=5)
ggplot(DF_toplot,                        
       aes(x = Time.point,
           y = Number.of.genes,
           fill = Group)) + 
  geom_bar(stat = "identity",
           position = "stack") +theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Number of differentially expressed genes")+
  facet_grid(~ Cell.type) + scale_fill_manual(values = c("#a7adba","#65737e","black"))
dev.off()

### Supplemental Figure 12B - DAR between AKI_GFP+ and control cells per cell type at 4 weeks and 6 months post AKI
#Make a list with all datasets
MyDF <- data.frame("Celltype"=c(),"Shared"=c(), "Concordant" = c(), "Notconcordant"=c(),
                   "4w_vs_4wCtrl"=c(),"6m_vs_6mCtrl"=c(), "Concordant_percent_of_6m"=c(),
                   "Notconcordant_percent_of_6m"=c(),"Concordant_percent_of_4w"=c(),
                   "Notconcordant_percent_of_4w"=c())

for(i in 1:length(Celltype)){
  AKI4w_vs_Ctrl4w <- read.csv(paste0("./GFPpos_vs_Ctrl/DAR/4weeks/2022_9_1_DAR_GFPpos4wvsCtrl4w_cluster",Celltype[[i]],"_logFC0_minpct0.05_padj0.01.csv"),header = T)
  AKI4w_vs_Ctrl4wall <- AKI4w_vs_Ctrl4w[which(abs(AKI4w_vs_Ctrl4w$avg_log2FC)>=0.1),"DAR"]
  AKI4w_vs_Ctrl4wup <- AKI4w_vs_Ctrl4w[which(AKI4w_vs_Ctrl4w$avg_log2FC>=0.1),"DAR"]
  AKI4w_vs_Ctrl4wdown <- AKI4w_vs_Ctrl4w[which(AKI4w_vs_Ctrl4w$avg_log2FC<=-0.1),"DAR"]
  AKI6m_vs_Ctrl6m <- read.csv(paste0("./GFPpos_vs_Ctrl/DAR/6months/2022_9_1_DAR_GFPpos6mvsCtrl6m_cluster",Celltype[[i]],"_logFC0_minpct0.05_padj0.01.csv"),header = T)
  AKI6m_vs_Ctrl6mall <- AKI6m_vs_Ctrl6m[which(abs(AKI6m_vs_Ctrl6m$avg_log2FC)>=0.1),"DAR"]
  AKI6m_vs_Ctrl6mup <- AKI6m_vs_Ctrl6m[which(AKI6m_vs_Ctrl6m$avg_log2FC>=0.1),"DAR"]
  AKI6m_vs_Ctrl6mdown <- AKI6m_vs_Ctrl6m[which(AKI6m_vs_Ctrl6m$avg_log2FC<=-0.1),"DAR"]
  #General overlap
  Shared <- length(intersect(AKI4w_vs_Ctrl4wall,AKI6m_vs_Ctrl6mall))
  Shared 
  Sharedup <- length(intersect(AKI4w_vs_Ctrl4wup,AKI6m_vs_Ctrl6mup))
  Shareddown <- length(intersect(AKI4w_vs_Ctrl4wdown,AKI6m_vs_Ctrl6mdown))
  Shareddown
  Concordant <- Sharedup + Shareddown
  Concordant
  Notconcordant <- Shared - Concordant
  Notconcordant
  
  myDF <- data.frame("Celltype"=Celltype[[i]],"Shared"=Shared, "Concordant" = Concordant, "Notconcordant"=Notconcordant,
                     "4w_vs_4wCtrl"=length(AKI4w_vs_Ctrl4wall),"6m_vs_6mCtrl"=length(AKI6m_vs_Ctrl6mall), "Concordant_percent_of_6m"=(Concordant/length(AKI6m_vs_Ctrl6mall))*100,
                     "Notconcordant_percent_of_6m"=(Notconcordant/length(AKI6m_vs_Ctrl6mall))*100,"Concordant_percent_of_4w"=(Concordant/length(AKI4w_vs_Ctrl4wall))*100,
                     "Notconcordant_percent_of_4w"=(Notconcordant/length(AKI4w_vs_Ctrl4wall))*100)
  MyDF <- rbind(MyDF, myDF)
}

#Add unique values
MyDF$Unique_6m <- MyDF$X6m_vs_6mCtrl - (MyDF$Concordant+MyDF$Notconcordant)
MyDF$Unique_4w <- MyDF$X4w_vs_4wCtrl - (MyDF$Concordant+MyDF$Notconcordant)

#Replace NA with 0
MyDF[is.na(MyDF)] <- 0

#Make DF to plot 
DF_toplot <- data.frame("Cell type"=rep(Celltype,3), "Group" = c(rep("Shared (concordant)",2*length(Celltype)),rep("Shared (nonconcordant)",2*length(Celltype)),
                                                                 rep("Unique",2*length(Celltype))) ,"Time point" = c(rep("4 weeks",length(Celltype)),rep("6 months",length(Celltype)),rep("4 weeks",length(Celltype)),
                                                                                                                     rep("6 months",length(Celltype)),rep("4 weeks",length(Celltype)),rep("6 months",length(Celltype))),
                        "Number of genes"=c(MyDF$Concordant,MyDF$Concordant,MyDF$Notconcordant,MyDF$Notconcordant,MyDF$Unique_4w,MyDF$Unique_6m))

#Make plot
DF_toplot$Cell.type <- factor(DF_toplot$Cell.type,levels =Celltype)
DF_toplot$Group <- factor(DF_toplot$Group,levels=c("Unique","Shared (concordant)","Shared (nonconcordant)"))

pdf("Supplemental_Figure_12B.pdf",width=15, height=5)
ggplot(DF_toplot,                         # Draw barplot with grouping & stacking
       aes(x = Time.point,
           y = Number.of.genes,
           fill = Group)) + 
  geom_bar(stat = "identity",
           position = "stack") +theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+
  ylab("Number of differentially accessible regions")+
  facet_grid(~ Cell.type) + scale_fill_manual(values = c("#a7adba","#65737e","black"))
dev.off()

### Supplemental Figure 12C - Volcano plot of DEG of AKI_GFP+ PTS3 versus Control PTS3 at 4 weeks and 6 months post AKI
## 4 weeks post AKI
# Read in list of differentially expressed genes of AKI_GFP+ PTS3 versus Control PTS3 at 4 weeks
DEG <- read.csv("./GFPpos_vs_Ctrl/DEG/4weeks/2022_9_1_DEG_GFPpos4wvsCtrl4w_clusterPTS3_logFC0_minpct0.1_padj0.01.csv")

# Volcano plot
DEG <- DEG[which(DEG$p_val_adj<0.01),]
DEG  <- mutate(DEG , color = case_when(DEG$avg_log2FC > 0.25 & DEG$p_val_adj < 0.01 ~ "Increased expression in AKI_GFP+ PTS3",
                                       DEG$avg_log2FC < -0.25 & DEG$p_val_adj < 0.01 ~ "Increased expression in Control PTS3",
                                       abs(DEG$avg_log2FC) < 0.25 & DEG$p_val_adj < 0.01 ~ "nonsignificant",
                                       DEG$p_val_adj > 0.01 ~ "nonsignificant"))

#Only label top 30 and selected genes
#select genes to be labeled on the Volcanoplot
DEG <- DEG %>% mutate(threshold = p_val_adj < 0.01) %>% arrange(p_val_adj) %>% mutate(volcanolabels = "")
DEG$volcanolabels[1:30] <- DEG$Gene[1:30]
DEG[which(DEG$p_val_adj==0),"p_val_adj"] <- 1e-300
GENES_toshow <- c("Baz2b","Nr1h4","Dab2","Gria3","Nrg1","Il34","Acy1","Slc22a13","Acox2","Ctnna2","Tgfbr1","Jarid2","Nr3c1","Nr3c2","Bach2","Slc9a8","Aqp1","Slc15a2","Acy3","Cdk6","Sash1","Igfbp7","Sema3c","Ppargc1b","Ccl28")
DEG[which(DEG$Gene%in%GENES_toshow),"volcanolabels"] <- DEG[which(DEG$Gene%in%GENES_toshow),"Gene"]

pdf("Supplemental_Figure_12C1.pdf")
ggplot(DEG, aes(x = avg_log2FC, y = -log10(p_val_adj), color=color,label = volcanolabels)) +
  geom_point(size = 1, alpha = 0.8, na.rm = T) + 
  geom_text_repel(aes(label=volcanolabels), max.overlaps = 1000,size=2,color="black") + 
  scale_color_manual(values =   c("#bc5f5f","#4aa4d5","darkgray"))+
  xlab("Average log2(Fold Change)") + 
  ylab("-log"[10]~"(p.adjust)") + theme_classic()+
  ggtitle("PTS3: AKI_GFP+_4weeks vs. Control_4weeks")+
  theme(plot.title = element_text(hjust = 0.5, size=12),axis.text = element_text(size = 10))
dev.off()

## 6 months post AKI
DEG <- read.csv("./GFPpos_vs_Ctrl/DEG/6months/2022_9_1_DEG_GFPpos6mvsCtrl6m_clusterPTS3_logFC0_minpct0.1_padj0.01.csv")

# Volcano plot
DEG <- DEG[which(DEG$p_val_adj<0.01),]
DEG  <- mutate(DEG , color = case_when(DEG$avg_log2FC > 0.25 & DEG$p_val_adj < 0.01 ~ "Increased expression in AKI_GFP+ PTS3",
                                       DEG$avg_log2FC < -0.25 & DEG$p_val_adj < 0.01 ~ "Increased expression in Control PTS3",
                                       abs(DEG$avg_log2FC) < 0.25 & DEG$p_val_adj < 0.01 ~ "nonsignificant",
                                       DEG$p_val_adj > 0.01 ~ "nonsignificant"))

#Label top 30 and selected genes
DEG <- DEG %>% mutate(threshold = p_val_adj < 0.01) %>% arrange(p_val_adj) %>% mutate(volcanolabels = "")
DEG$volcanolabels[1:30] <- DEG$Gene[1:30]
DEG[which(DEG$p_val_adj==0),"p_val_adj"] <- 1e-300
GENES_toshow <- c("Baz2b","Nr1h4","Dab2","Gria3","Nrg1","Il34","Acy1","Slc22a13","Acox2","Ctnna2","Tgfbr1","Jarid2","Nr3c1","Nr3c2","Bach2","Slc9a8","Aqp1","Slc15a2","Acy3","Cdk6","Sash1","Igfbp7","Sema3c","Ppargc1b","Ccl28")
DEG[which(DEG$Gene%in%GENES_toshow),"volcanolabels"] <- DEG[which(DEG$Gene%in%GENES_toshow),"Gene"]

pdf("Supplemental_Figure_12C2.pdf")
ggplot(DEG, aes(x = avg_log2FC, y = -log10(p_val_adj), color=color,label = volcanolabels)) +
  geom_point(size = 1, alpha = 0.8, na.rm = T) + 
  geom_text_repel(aes(label=volcanolabels), max.overlaps = 1000,size=2,color="black") + 
  scale_color_manual(values =   c("#bc5f5f","#4aa4d5","darkgray"))+
  xlab("Average log2(Fold Change)") + 
  ylab("-log"[10]~"(p.adjust)") + theme_classic()+
  ggtitle("PTS3: AKI_GFP+_6months vs. Control_6months")+
  theme(plot.title = element_text(hjust = 0.5, size=12),axis.text = element_text(size = 10))
dev.off()

#################################### Supplemental Figure 15 #################################### 
##  Supplemental Figure 15A - overlap of upregulated between AKI_GFP+ vs. Control and AKI_GFP- vs. Control cells

#Make a list with all datasets
#Read in all DEG IRIGFP vs Ctrl gene lists into one list
MyDF_Genes <- data.frame("Celltype"=c(),"Overlap_genes"=c())
MyDF <- data.frame("Celltype"=c(),"Overlap"=c(),"Percent_of_GFPpos"=c(),"GFPpos_vs_Ctrl"=c(),"GFPneg_vs_Ctrl"=c())

for(i in 1:length(Celltype)){
  DEGpos <- read.csv(paste0("./GFPpos_vs_Ctrl/DEG/2022_3_31_DEG_GFPposvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  DEGpos <- DEGpos[which(DEGpos$avg_log2FC>=0.25),"Gene"]
  DEGneg <- read.csv(paste0("./GFPneg_vs_Ctrl/DEG/2022_3_31_DEG_GFPnegvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  DEGneg <- DEGneg[which(DEGneg$avg_log2FC>=0.25),"Gene"]
  Shared <- intersect(DEGpos,DEGneg)
  Percent_of_GFP <- (length(Shared)/length(DEGpos))*100
  myDF <- data.frame("Celltype"=Celltype[[i]],"Overlap"=length(Shared),"Percent_of_GFPpos"=Percent_of_GFP,"GFPpos_vs_Ctrl"=length(DEGpos),"GFPneg_vs_Ctrl"=length(DEGneg))
  MyDF <- rbind(MyDF, myDF)
  MyGenes <- data.frame("Celltype"=rep(Celltype[[i]],length(Shared)),"Overlap_genes"=Shared)
  MyDF_Genes <- rbind(MyDF_Genes,MyGenes)
}

MyDF$GFPpos_only <- MyDF$GFPpos_vs_Ctrl-MyDF$Overlap

#Make dataframe to plot
DF_toplot <- data.frame("Number of upregulated genes"=c(MyDF$GFPpos_only,MyDF$Overlap) ,"Celltype"=c(MyDF$Celltype,MyDF$Celltype),
                        "Group"=c(rep("Up in AKI_GFP+ only",length(Celltype)),rep("Up in AKI_GFP+ & AKI_GFP- vs. Control",length(Celltype))))
MyDF$Celltype <- factor(MyDF$Celltype, levels=rev(Celltype))
DF_toplot$Group <- factor(DF_toplot$Group, levels=rev(c("Up in AKI_GFP+ & AKI_GFP- vs. Control","Up in AKI_GFP+ only")))

pdf("Supplemental_Figure_15A.pdf")
ggplot(DF_toplot, aes(fill=Group, x=Number.of.upregulated.genes, y=Celltype))+geom_bar(position="stack", stat="identity")+ scale_fill_manual(values = c("#c6c9cb",	"#71777e"))+
  xlab("No. of upregulated genes")+theme_classic()+ylab("Cell type")
dev.off()

##  Supplemental Figure 15B - overlap of downregulated between AKI_GFP+ vs. Control and AKI_GFP- vs. Control cells

MyDF_Genes <- data.frame("Celltype"=c(),"Overlap_genes"=c())
MyDF <- data.frame("Celltype"=c(),"Overlap"=c(),"Percent_of_GFPpos"=c(),"GFPpos_vs_Ctrl"=c(),"GFPneg_vs_Ctrl"=c())

for(i in 1:length(Celltype)){
  DEGpos <- read.csv(paste0("./GFPpos_vs_Ctrl/DEG/2022_3_31_DEG_GFPposvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  DEGpos <- DEGpos[which(DEGpos$avg_log2FC<=-0.25),"Gene"]
  DEGneg <- read.csv(paste0("./GFPneg_vs_Ctrl/DEG/2022_3_31_DEG_GFPnegvsCtrl_cluster",Celltype[[i]],"_logFC0_minpct0.1_padj0.01.csv"),header = T)
  DEGneg <- DEGneg[which(DEGneg$avg_log2FC<=-0.25),"Gene"]
  Shared <- intersect(DEGpos,DEGneg)
  Percent_of_GFP <- (length(Shared)/length(DEGpos))*100
  myDF <- data.frame("Celltype"=Celltype[[i]],"Overlap"=length(Shared),"Percent_of_GFPpos"=Percent_of_GFP,"GFPpos_vs_Ctrl"=length(DEGpos),"GFPneg_vs_Ctrl"=length(DEGneg))
  MyDF <- rbind(MyDF, myDF)
  MyGenes <- data.frame("Celltype"=rep(Celltype[[i]],length(Shared)),"Overlap_genes"=Shared)
  MyDF_Genes <- rbind(MyDF_Genes,MyGenes)
}

MyDF$GFPpos_only <- MyDF$GFPpos_vs_Ctrl-MyDF$Overlap
MyDF$Celltype <- factor(MyDF$Celltype, levels=rev(Celltype))

#Make dataframe to plot
DF_toplot <- data.frame("Number of downregulated genes"=c(MyDF$GFPpos_only,MyDF$Overlap) ,"Celltype"=c(MyDF$Celltype,MyDF$Celltype),
                        "Group"=c(rep("Down in AKI_GFP+ only",length(Celltype)),rep("Down in AKI_GFP+ & AKI_GFP- vs. Control",length(Celltype))))
DF_toplot$Group <- factor(DF_toplot$Group, levels=rev(c("Down in AKI_GFP+ & AKI_GFP- vs. Control","Down in AKI_GFP+ only")))

pdf("Supplemental_Figure_15B.pdf")
ggplot(DF_toplot, aes(fill=Group, x=Number.of.downregulated.genes, y=Celltype))+geom_bar(position="stack", stat="identity")+ scale_fill_manual(values = c("#c6c9cb",	"#71777e"))+
  xlab("No. of downregulated genes")+theme_classic()+ylab("Cell type")
dev.off()