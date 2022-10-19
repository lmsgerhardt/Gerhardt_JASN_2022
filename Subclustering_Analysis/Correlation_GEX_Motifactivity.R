#Title: Correlation of GEX and Motif activity of differentially expressed genes and differentially active motifs
#Code adapted from Muto et al., Nature Communications, 2021, PMID: 33850129 (https://github.com/p4rkerw/Muto_Wilson_NComm_2020/blob/master/analysis/corr_chromVar_TF_exp.R)

library(Seurat)
library(Signac)
library(dplyr)
library(tidyverse)
library(ggplot2)
set.seed(1234)

#Read in the final R object of AKI_GFP+ ("PTS1_1","PTS1_2","PTS1-S2","PTS2","PTS3","Injured PT","LOH-TL-C","LOH-TL-JM") subclustering analysis 
PTC <- readRDS("./Subclustering_Analysis/PTC_FINAL_SUB.Rds")

DefaultAssay(PTC) <- "chromvar"

# create a data frame of motifs and corresponding gene names
motifs <- PTC$peaks@motifs@motif.names

motif_names <- data.frame(genes=unlist(motifs)) %>% rownames_to_column(var="motif") %>% arrange(genes) 

# format the motif_names df so each motif corresponds to a single gene
# the JASPAR database has multiple genes associated with each motif separated by "::"
motif_names <- distinct(motif_names) %>% #subset unique rows
  arrange(genes) %>%
  dplyr::filter(!str_detect(genes, pattern="::"))%>%
  dplyr::filter(!str_detect(genes, pattern="var.")) %>% # remove any genes associated with multiple variants
  as.data.frame()

# fix the one occurrence of EWSR1-FLI1
motif_names$genes <- str_replace(motif_names$genes, pattern="EWSR1-FLI1", replacement="EWSR1")

#Make Motif names lower case to be compatible with gene annotation
motif_names$genes <-paste0(substr(motif_names$genes,1,1),tolower(substr(motif_names$genes,2,nchar(motif_names$genes)))) 

# use the findmarkers function to calculate gene expression and activity for every motif-gene combination
motif.ls = motif_names$motif

df.ls <- lapply(motif.ls, function(motif) {
  print(motif)
  gene <- motif_names[motif_names$motif == motif,]$genes
  print(gene)
  DefaultAssay(PTC) <- "chromvar"
  aver_chromvar <- FindAllMarkers(
    object = PTC,
    only.pos = FALSE,
    features=motif,
    test.use = 'LR',
    min.pct = 0,
    mean.fxn = rowMeans,
    logfc.threshold = 0 #find all cluster-specific motifs
  )
  # compute the average expression from the rna assay
  DefaultAssay(PTC) <- "RNA"
  PTC <- NormalizeData(PTC)
  aver_exp <- FindAllMarkers(PTC, assays="RNA",min.pct=0, logfc.threshold=0, features=gene, verbose=FALSE, only.pos = FALSE)
  
  # join the matrices and return null if gene expression is not detected for all cell types
  mat <- tryCatch(full_join(aver_chromvar, aver_exp, by = "cluster"),
                  error=function(e) NULL)
  if(!is.null(mat)) {
    df <- data.frame(chromvar=mat$avg_log2FC.x, rna=mat$avg_log2FC.y, celltype=mat$cluster,
                     motif=mat$gene.x, gene=mat$gene.y, chrom_pval=mat$p_val_adj.x, gene_pval=mat$p_val_adj.y)
  } else {
    return(NULL)
  }
})

# aggregate the motif-gene df
df <- bind_rows(df.ls)
df <- dplyr::arrange(df, gene, motif)
write.csv(df, file = "./Subclustering_Analysis/Correlation_chromVar_GEX.csv", row.names = F)

sessionInfo()