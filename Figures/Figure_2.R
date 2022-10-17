#Title: Plots for Figure 2

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)
library(WGCNA)

setwd()

#Read in the Seurat object from Kirita et al., PNAS, 2020, PMID: 32571916 which was kindly provided by the Humphreys Lab. 
load("./Kirita_et_al/IRI.all.RData")
IRI.all

#################################### Figure 2A - Heatmap of DEG in cell type-specific comparison of 12h post AKI versus control data #################################### 
### Identify differentially expressed genes in per cell type comparison of 12h post AKI versus control data
#Subset 12h post AKI and control data
Idents(IRI.all) <- "time"

Ctrl_12h <- subset(IRI.all, idents=c("Control", "12hours")) 
Idents(Ctrl_12h) <- "name"
celltypes  <- levels(Ctrl_12h)

#Check number of cells per cell type and time point
Ctrl_12h@meta.data$name_time <- paste0(Ctrl_12h@meta.data$name,sep="_",Ctrl_12h@meta.data$time)
Idents(Ctrl_12h) <- "name_time"
table(Idents(Ctrl_12h))

# Identify differentially expressed genes between 12h post AKI and control per cluster
celltypes <- celltypes[which(!celltypes=="MD")] #exclude MD from analysis, because Ctrl_12h contains < 50 cells per group
Idents(Ctrl_12h) <- "name"
Ctrl_12h <- subset(Ctrl_12h,idents=celltypes)

DEG_all_up <- data.frame("p_val"=c(),"avg_log2FC"=c(),"pct.1"=c(),"pct.2"=c(),"p_val_adj"=c(),"Gene"=c())
DEG_all_down <- data.frame("p_val"=c(),"avg_log2FC"=c(),"pct.1"=c(),"pct.2"=c(),"p_val_adj"=c(),"Gene"=c())

for(i in 1:length(celltypes)){
  SUB <- subset(Ctrl_12h,idents=celltypes[i])
  SUB <- NormalizeData((SUB), normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(SUB) <- "time"
  DEG_12hvsCtrl <- FindMarkers(SUB, ident.1 = "12hours", ident.2="Control", min.pct = 0.05, logfc.threshold = 0.25,only.pos = FALSE)
  DEG_12hvsCtrl$Gene <- row.names(DEG_12hvsCtrl)
  DEG_12hvsCtrl <- DEG_12hvsCtrl[DEG_12hvsCtrl$p_val_adj < 0.01,]
  write.csv(DEG_12hvsCtrl,paste0("./Kirita_et_al/2021_08_04_DEG_",paste(celltypes[i],collapse="_"),"_12hvsCtrl_logFC0.25_pct10.05_adjpval_0.01.csv"),row.names = F)
  DEG_12hvsCtrl_up <- DEG_12hvsCtrl[which(DEG_12hvsCtrl$avg_log2FC>=0.25),]
  DEG_12hvsCtrl_down <- DEG_12hvsCtrl[which(DEG_12hvsCtrl$avg_log2FC<=-0.25),]
  DEG_all_up <- rbind(DEG_all_up, DEG_12hvsCtrl_up)
  DEG_all_down <- rbind(DEG_all_down, DEG_12hvsCtrl_down)
}

write.csv(DEG_all_up,"./Kirita_et_al/IRI.all_12hpIRIvsCtrl_per_celltype_logFC0.25_pct0.05_padj0.01_UP.csv",row.names = T)
write.csv(DEG_all_down,"./Kirita_et_al/IRI.all_12hpIRIvsCtrl_per_celltype_logFC0.25_pct0.05_padj0.01_DOWN.csv",row.names = T)

# Get average expression per cluster
Idents(Ctrl_12h) <- "name_time"
Ctrl_12h <- NormalizeData(Ctrl_12h)
cluster.averages.up <- data.frame(AverageExpression(Ctrl_12h, features = unique(DEG_all_up$Gene)))
colnames(cluster.averages.up) <- substr(colnames(cluster.averages.up),5,20)
cluster.averages.down <- data.frame(AverageExpression(Ctrl_12h, features = unique(DEG_all_down$Gene)))
colnames(cluster.averages.down) <- substr(colnames(cluster.averages.down),5,20)
new_order <- c("Pod_Control","PEC_Control","PTS1_Control","PTS2_Control","PTS3_Control","NewPT1_Control","NewPT2_Control","DTL.ATL_Control","MTAL_Control","CTAL1_Control","CTAL2_Control",
               "DCT_Control","DCT.CNT_Control","CNT_Control","PC1_Control","PC2_Control","ICA_Control","ICB_Control","Uro_Control","EC1_Control","EC2_Control","Fib_Control","Per_Control","Mø_Control","Tcell_Control",
               "Pod_12hours","PEC_12hours","PTS1_12hours","PTS2_12hours","PTS3_12hours","NewPT1_12hours","NewPT2_12hours","DTL.ATL_12hours","MTAL_12hours","CTAL1_12hours","CTAL2_12hours",
               "DCT_12hours","DCT.CNT_12hours","CNT_12hours","PC1_12hours","PC2_12hours","ICA_12hours","ICB_12hours","Uro_12hours","EC1_12hours","EC2_12hours","Fib_12hours","Per_12hours","Mø_12hours","Tcell_12hours"
)
cluster.averages.up <- cluster.averages.up[,new_order]
write.csv(cluster.averages.up,"./Kirita_et_al/IRI.all_12hpIRIvsCtrl_per_celltype_cluster.averages.up.csv",row.names = T)
cluster.averages.down <- cluster.averages.down[,new_order]
write.csv(cluster.averages.down,"./Kirita_et_al/IRI.all_12hpIRIvsCtrl_per_celltype_cluster.averages.down.csv",row.names = T)

### WGCNA analysis on downregulated genes
#Transpose the data for further analysis
cluster.averages.down <- t(cluster.averages.down)
dim(unique(cluster.averages.down))
#Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cluster.averages.down, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#We choose the power 7

cluster.averages.down_net = blockwiseModules(cluster.averages.down, power = 7,
                                           TOMType = "unsigned", minModuleSize = 20,
                                           reassignThreshold = 0, mergeCutHeight = 0.25,
                                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                                           saveTOMs = TRUE,
                                           saveTOMFileBase = "IRI_12h_DOWN",
                                           verbose = 3)

table(cluster.averages.down_net$colors)

#Saving the data
cluster.averages.down_moduleLabels = cluster.averages.down_net$colors
cluster.averages.down_moduleColors = labels2colors(cluster.averages.down_net$colors)
cluster.averages.down_MEs = cluster.averages.down_net$MEs;
cluster.averages.down_geneTree = cluster.averages.down_net$dendrograms[[1]];
save(cluster.averages.down_MEs, cluster.averages.down_moduleLabels, cluster.averages.down_moduleColors, cluster.averages.down_geneTree,
     file = "./Kirita_et_al/IRI_12h_vs_Ctrl_DOWN-networkConstruction-auto.RData")

#Get gene names of all modules and plot genes and modules across clusters
modulenumber <- as.list(c(0:18))
cluster.averages.down_modules <- as.list(modulenumber)

for(i in 0:length(modulenumber)){
  cluster.averages.down_modules[[i+1]] <- names(cluster.averages.down_moduleLabels)[cluster.averages.down_moduleLabels==i]
}

cluster.averages.down_modules <- cluster.averages.down_modules[c(1:19)]

cluster.averages.down_modules_genes <- data.frame("Gene"=c(),"Module"=c())

for(i in 1:length(cluster.averages.down_modules)){
  Genes <- cluster.averages.down_modules[[i]]
  df <- data.frame("Gene"=Genes,"Module"=rep((i-1),length(Genes)))
  cluster.averages.down_modules_genes <- rbind(cluster.averages.down_modules_genes,df)
}

cluster.averages.down_modules_genes <- cluster.averages.down_modules_genes[-1,]
write.csv(cluster.averages.down_modules_genes ,"./Kirita_et_al/WGCNA_downregulated_modules_genes.csv",row.names = F)

#Reorder data according to modules
cluster.averages.down <- t(cluster.averages.down)
cluster.averages.down <- cluster.averages.down[as.character(cluster.averages.down_modules_genes$Gene),]
head(cluster.averages.down)

row_anno <- as.data.frame(cluster.averages.down_modules_genes[,"Module"])
row.names(row_anno) <- cluster.averages.down_modules_genes$Gene
row_anno$`cluster.averages.down_modules_genes[, "Module"]` <- as.character(row_anno$`cluster.averages.down_modules_genes[, "Module"]`)

annotdf <- data.frame(row.names = row.names(cluster.averages.down), Module = as.character(row_anno$`cluster.averages.down_modules_genes[, "Module"]`))

#Plot modules across clusters
pdf("Figure_2A1.pdf",height=8)
pheatmap(mat = t(cluster.averages.down), 
         cluster_cols = F,
         color = colorRampPalette(colors = c("blue", "white", "red"))(255),
         scale = "column", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotdf,
         #cellwidth = 15, # Make the cells wider
         show_colnames = F,
         show_rownames = T,
         cluster_rows = F)
dev.off()

### WGCNA analysis on upregulated genes
#Transpose the data for further analysis
cluster.averages.up <- t(cluster.averages.up)

#Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cluster.averages.up, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#We choose the power 7

cluster.averages.up_net = blockwiseModules(cluster.averages.up, power = 7,
                              TOMType = "unsigned", minModuleSize = 20,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "IRI_12h_UP",
                              verbose = 3)

table(cluster.averages.up_net$colors)

#Saving the data
cluster.averages.up_moduleLabels = cluster.averages.up_net$colors
cluster.averages.up_moduleColors = labels2colors(cluster.averages.up_net$colors)
cluster.averages.up_MEs = cluster.averages.up_net$MEs;
cluster.averages.up_geneTree = cluster.averages.up_net$dendrograms[[1]];
save(cluster.averages.up_MEs, cluster.averages.up_moduleLabels, cluster.averages.up_moduleColors, cluster.averages.up_geneTree,
     file = "./Kirita_et_al/IRI_12h_vs_Ctrl_UP-networkConstruction-auto.RData")

#Get gene names of all modules and plot genes and modules across clusters
modulenumber <- as.list(c(0:22))
cluster.averages.up_modules <- as.list(modulenumber)

for(i in 0:length(modulenumber)){
  cluster.averages.up_modules[[i+1]] <- names(cluster.averages.up_moduleLabels)[cluster.averages.up_moduleLabels==i]
}

cluster.averages.up_modules <- cluster.averages.up_modules[c(1:23)]

cluster.averages.up_modules_genes <- data.frame("Gene"="1","Module"=1)

for(i in 1:length(cluster.averages.up_modules)){
  Genes <- cluster.averages.up_modules[[i]]
  df <- data.frame("Gene"=Genes,"Module"=rep((i-1),length(Genes)))
  cluster.averages.up_modules_genes <- rbind(cluster.averages.up_modules_genes,df)
}

cluster.averages.up_modules_genes <- cluster.averages.up_modules_genes[-1,]
write.csv(cluster.averages.up_modules_genes ,"./Kirita_et_al/WGCNA_upregulated_modules_genes.csv",row.names = F)

#Reorder data according to modules
cluster.averages.up <- t(cluster.averages.up)
cluster.averages.up <- cluster.averages.up[as.character(cluster.averages.up_modules_genes$Gene),]
head(cluster.averages.up)

row_anno <- as.data.frame(cluster.averages.up_modules_genes[,"Module"])
row.names(row_anno) <- cluster.averages.up_modules_genes$Gene
row_anno$`cluster.averages.up_modules_genes[, "Module"]` <- as.character(row_anno$`cluster.averages.up_modules_genes[, "Module"]`)

annotdf <- data.frame(row.names = row.names(cluster.averages.up), Module = as.character(row_anno$`cluster.averages.up_modules_genes[, "Module"]`))

#Plot modules across clusters
pdf("Figure_2A2.pdf",height=8)
pheatmap(mat = t(cluster.averages.up), 
         cluster_cols = F,
         color = colorRampPalette(colors = c("blue", "white", "red"))(255),
         scale = "column", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotdf,
         #cellwidth = 15, # Make the cells wider
         show_colnames = F,
         show_rownames = T,
         cluster_rows = F)
dev.off()

#################################### Figure 2B - bar plot of gene set enrichment analysis of cell type-specific marker genes amongst genes that are downregulated at 12 hours post IRI #################################### 
# usage: gene set enrichment analysis of cell type-specific marker genes amongst genes that are downregulated at 12 hours post IRI 
# top 20 differentially expressed genes per cell type in control samples will be used as as cell-type-specific marker gene set

#Get list of all genes in the dataset as reference data set
myFeatures <- row.names(Ctrl_12h@assays$RNA@data)

###Get marker gene set (top 20 differentially expressed genes between cell types in control samples)
#subset Control samples only
Idents(Ctrl_12h) <- "time"
Ctrl <- subset(Ctrl_12h, idents="Control") 
Idents(Ctrl) <- "name"

Ctrl <- NormalizeData(Ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
Ctrl_DEG <- FindAllMarkers(Ctrl, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
Ctrl_DEG <- Ctrl_DEG[which(Ctrl_DEG$p_val_adj < 0.01),]
allmarkergenes <- Ctrl_DEG[which(Ctrl_DEG$cluster%in%celltypes),]
allmarkergenes <- as.data.frame(allmarkergenes %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))

#### Gene set enrichment analysis of marker gene sets in genes that are downregulated in 12h post IRI versus control data using Fisher's exact test
Fisher_DF <- data.frame("cluster"=c(),"p value"=c(),"A"=c(),"B"=c(),"C"=c(),"D"=c())

for(i in 1:length(celltypes)){
  querygenes <- read.csv(paste0("./Kirita_et_al/2021_08_04_DEG_",paste(celltypes[i],collapse="_"),"_12hvsCtrl_logFC0.25_pct10.05_adjpval_0.01.csv"))
  querygenes<- querygenes[which(querygenes$avg_log2FC<0),"Gene"]
  markergenes <- as.character(allmarkergenes[which(allmarkergenes$cluster ==celltypes[i]),"gene"])
  A <- querygenes[which(querygenes %in% markergenes)]
  B <- length(myFeatures[which(myFeatures %in% markergenes & !myFeatures %in% A)])
  A <- length(A)  
  
  C <- querygenes[which(!querygenes %in% markergenes)]
  D <- length(myFeatures[which(!myFeatures %in% markergenes & !myFeatures %in% C)])
  C <- length(C)
  
  print(sum(A,B,C,D))
  
  myCT <- matrix(c(A,B,C,D), ncol=2,nrow=2)
  test.result <- fisher.test(myCT,  y = NULL, hybrid = FALSE, #all default values
                             hybridPars = c(expect = 5, percent = 80, Emin = 1),
                             control = list(), or = 1, alternative = "two.sided",
                             conf.int = TRUE, conf.level = 0.95,
                             simulate.p.value = FALSE, B = 2000)
  newDF <- data.frame("p value"= test.result$p.value, "cluster"= celltypes[i],"A"=A,"B"=B,"C"=C,"D"=D)
  Fisher_DF <- rbind(Fisher_DF, newDF)
  rm(querygenes,A,B,C,D,newDF)
}

# adjust P values for multiple testing 
padj <- p.adjust(Fisher_DF$p.value, method = "BH", n = length(Fisher_DF$p.value))
Fisher_DF$p_val_adjust <- padj
Fisher_DF$minlog10padjust <- -log10(Fisher_DF$p_val_adjust)
Fisher_DF$Percent <- (Fisher_DF$A/20)*100 #Calculate the percentage of marker genes in the query data set

#Make a barplot with the percentage of marker genes in downregulated genes per cluster. Color the bar by padjust
Fisher_DF$cluster <- factor(Fisher_DF$cluster,levels=rev(c("Pod","PEC","PTS1","PTS2","PTS3","NewPT1","NewPT2","DTL-ATL","MTAL","CTAL1","CTAL2","DCT",
                                                           "DCT-CNT","CNT","PC1","PC2","ICA","ICB","Uro","EC1","EC2","Fib","Per","Mø","Tcell")))

pdf("Figure_2B.pdf")
ggplot(Fisher_DF, aes(Percent, cluster, fill = minlog10padjust)) + geom_bar(stat = "identity") + scale_fill_continuous("-log10(p.adjust)",low="blue", high="red")+
  xlab("Percent of marker genes [%]")+ylab("Cell type") + theme_classic()
dev.off()


