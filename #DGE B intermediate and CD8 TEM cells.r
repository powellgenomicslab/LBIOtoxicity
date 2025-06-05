#DGE B intermediate and CD8 TEM cells

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(jsonlite)
library(data.table)
library(devtools)
library(Azimuth)
library(ComplexHeatmap)
library(DESeq2)
library(pheatmap)
library(viridis)
library(EnhancedVolcano)
library(dittoSeq)
library(presto)
library(tictoc)

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
LBtox
colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"


#just baseline cells
Baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))
Idents(Baseline) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Idents(Baseline) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#baseline no LB142 and LB035 (female
Baseline.New <- subset(Baseline.no142, subset=Sex==c("Male"))

#baseline B intermediate cells subset
Baseline.Bint.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("B intermediate"))
Idents(Baseline.Bint.no142) <- "Population"
Baseline.Bcells

#baseline analysis regressing out effect of sex
DefaultAssay(Baseline.Bint.no142) <- "RNA"
DEGS_Bint_regress <- FindMarkers(Baseline.Bint.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
head(DEGS_Bint_regress)
DEGS_Bint_regress <- DEGS_Bint_regress %>%
            tibble::rownames_to_column(var="Gene")
write.csv(DEGS_Bint_regress, "B Intermediate Differentially expressed Genes Baseline cohort vs control.CSV") 

#Repeat volcano plot - regressing out effect of sex
plot1 <- EnhancedVolcano(DEGS_Bint_regress, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bint_regress$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot B intermediate Cells Baseline cohort vs control.png", height = 7.5, width = 7.5)

#repeat baseline analysis with males only
DefaultAssay(Baseline.New.Bcells) <- "RNA"
Baseline.New.Bcells <- subset(x=Baseline.New, subset=predicted.celltype.l2==c("B intermediate"))
Idents(Baseline.New.Bcells) <- "Population"
DEGS_Bint_Males <- FindMarkers(Baseline.New.Bcells, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
head(DEGS_Bint_Males)
DEGS_Bint_Males <- DEGS_Bint_Males %>%
            tibble::rownames_to_column(var="Gene")

write.csv(DEGS_Bint_Males, "B Intermediate Differentially expressed Genes Males Only.CSV")   
DEGS_Bint <- read.csv("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/B intermediate cell markers at BASELINE control vs cohort nil LB142.csv")
DEGS_Bint

#Male only plot
TopGenes <- rownames(DEGS_Bint_Males[order(DEGS_Bint_Males$p_val_adj, decreasing=FALSE)[1:25],])
DEGS_Bint_Males$TopGene <- ifelse(DEGS_Bint_Males$Gene %in% TopGenes, TRUE, FALSE)


#Enhanced Volcano Plot - youtube tutorial
plot2 <- EnhancedVolcano(DEGS_Bint_Males, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bint_Males$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot B intermediate Cells top 10 genes Malesv2.png", height = 7.5, width = 7.5)

#B intermediate cells DGE
DEGS_Bint <- FindMarkers(Baseline.Bcells.no142, ident.1= "Cohort", ident.2="Control",  vars.to.regress=c("Pool", "Sex"))
head(DEGS_Bint)
DEGS_Bint <- DEGS_Bint %>%
            tibble::rownames_to_column(var="Gene")

DEGS_Bint
write.csv(DEGS_Bint, "B Intermediate Differentially expressed Genes no LB142.CSV")     

getwd()

##select the top25 genes for plotting
TopGenes <- rownames(DEGS_Bint[order(DEGS_Bint$p_val_adj, decreasing=FALSE)[1:25],])
DEGS_Bint$TopGenes <- ifelse(DEGS_Bint$Gene %in% TopGenes, TRUE, FALSE)

DefaultAssay(Baseline.Bcells) <- "RNA"

#Re-analyse ALL cell types
Idents(LBtox$predicted.celltype.l2) <- "Cohort"
ResponseMarkers <- FindAllMarkers(mostpbmc, only.pos=TRUE)
head(ResponseMarkers)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5, min.cells.per.ident > 20)
dim(ResponseMarkersFilter)
plot11 <- DoHeatmap(subset(mostpbmc_celltypes$`Cytotoxic_CD8+_Tcell_S100B+`, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot11, filename = "Heatmap_DGE_CytotoxicCD8.png", height = 4.5, width = 7.5)

#create new function for all cell types
find_markers<- function(x,celltype){

cat("Analyzing group ", celltype, "\n")

#repeat code for all cell types
Idents(x) <- "Responder"
ResponseMarkers <- FindAllMarkers(mostpbmc, only.pos=TRUE)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5)
plot11 <- DoHeatmap(subset(x, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot11, filename = paste("Heatmap_DGE_cells_",celltype,".png"), height = 4.5, width = 7.5


#Heatmap
Heatmap.Data <- GetAssayData(Baseline.Bcells.no142, slot='scale.data', features=TopGenes)

#Volcano Plot - Lachlan's code
library(ggrepel)

plot1 <- ggplot(DEGS_Bint, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=TopGene)) +
    geom_text_repel(data=DEGS_Bint[1:10,], aes(x=avg_log2FC, y=-log10(p_val_adj), label=Gene)) +
    ggtitle("Volcano Plot: B intermediate cells") +
    xlab("log2 fold change") +
    ylab("-log10 FDR") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  scale_color_discrete(name = "DEG")

  ggsave(plot1,filename="Volcano Plot B intermediate Cells top 10 genes.png", height = 7.5, width = 7.5)

#Enhanced Volcano Plot - youtube tutorial
plot1 <- EnhancedVolcano(DEGS_Bint, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bint$Gene)
  ggsave(plot1,filename="Volcano Plot B intermediate Cells top 10 genes.png", height = 7.5, width = 7.5)

#Heatmap B intermediate Cells

#convert to dataframe -sanbiomics youtube
DEGS_Bint.df <- as.data.frame(DEGS_Bint)
mat <- counts(DEGS_Bint.df, normalized=TRUE)
t(apply(DEGS_Bint.df, 1, scale))

plot2 <- Heatmap(DEGS_Bint.df, cluster_rows=T, cluster_columns=T, column_labels=colnames(DEGS_bint.df), name="Z score" )


#pheatmap for B intermediate cell genes
plot3 <- pheatmap(df.bint, scale="row", clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="complete", color=colorRampPalette(c("blue", "white", "red"))(50))
ggsave(plot3, filename="Heatmap for B intermediate cell genes.png", height=20, width=10)

#baseline CD8 TEM cells
Baseline.CD8TEMcells.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD8 TEM"))
Idents(Baseline.CD8TEMcells.no142) <- "Population"

#CD8 TEM cells DGE
CD8TEM.markers.DGE.no142 <- FindMarkers(Baseline.CD8TEMcells.no142, ident.1= "Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
head(CD8TEM.markers.DGE.no142)
tail(CD8TEM.markers.DGE.no142)
CD8TEM.markers.DGE.no142 <- CD8TEM.markers.DGE.no142 %>%
            tibble::rownames_to_column(var="Gene")
head(CD8TEM.markers.DGE.no142)

write.csv(CD8TEM.markers.DGE.no142, "CD8 TEM Differentially expressed Genes Baseline cohort vs control.CSV")   

#Volcano Plot CD8 TEM
plot3 <- EnhancedVolcano(CD8TEM.markers.DGE.no142, x="avg_log2FC", y="p_val_adj", lab=CD8TEM.markers.DGE.no142$Gene, FCcutoff=0.1)
  ggsave(plot3,filename="Volcano Plot CD8 TEM cells DEGs Baseline.png", height = 7.5, width = 7.5)

#save CSV
write.csv(CD8TEM.markers.DGE.no142, "CD8 TEM cell markers at BASELINE control vs cohort nil LB142.csv")

#Heatmap CD8TEM cells
plot2 <- DoHeatmap(subset(CD8TEM.markers.DGE.no142, downsample = 100), features = features, size = 3)
