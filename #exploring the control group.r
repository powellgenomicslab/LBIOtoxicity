#exploring the control group

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

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
LBtox
colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Idents(LBtox.no142) <- "predicted.celltype.l2"
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Cohort.no142 <- subset(x=LBtox.no142, subset=Population==c("Cohort"))
Control.no142 <- subset(x=LBtox.no142, subset=Population==c("Control"))
Idents(Baseline.no142) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#Baseline Control Cells
Baseline <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Baseline.Control <- subset(x=Baseline, subset=Population==c("Control"))

#Control Cells
Control <- subset(x=LBtox.no142, subset=Population==c("Control"))

#specific cell subtypes

#NK cells
NK.Control <- subset(x=Control, subset=predicted.celltype.l2==c("NK"))
Idents(NK.Control) <- "Timepoint"
DefaultAssay(NK.Control) <- "RNA"
DEGS_NKControl <- FindMarkers(NK.Control, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_NKControl <- DEGS_NKControl %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_NKControl)

write.csv(DEGS_NKControl, "NK cells Differentially expressed Genes Control C1 vs C2.CSV")   

#volcano plot NK cells
plot1 <- EnhancedVolcano(DEGS_NKControl, x="avg_log2FC", y="p_val_adj", lab=DEGS_NKControl$Gene, FCcutoff=0.5)
  ggsave(plot1,filename="Volcano Plot NK cells control group C1 vs C2.png", height = 7.5, width = 7.5)

#CD14 mono
CD14.Control <- subset(x=Control, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(CD14.Control) <- "Timepoint"
DefaultAssay(CD14.Control) <- "RNA"
DEGS_CD14mono <- FindMarkers(CD14.Control, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD14mono <- DEGS_CD14mono %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD14mono)

write.csv(DEGS_CD14mono, "CD14 Mono cells Differentially expressed Genes Control C1 vs C2.CSV")   

#volcano plot CD14 Mono cells
plot2 <- EnhancedVolcano(DEGS_CD14mono, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD14mono$Gene, FCcutoff=0.5)
  ggsave(plot2,filename="Volcano Plot CD14 Mono cells control group C1 vs C2.png", height = 7.5, width = 7.5)

#CD16 mono
CD16.Control <- subset(x=Control, subset=predicted.celltype.l2==c("CD16 Mono"))
Idents(CD16.Control) <- "Timepoint"
DefaultAssay(CD16.Control) <- "RNA"
DEGS_CD16mono <- FindMarkers(CD16.Control, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD16mono <- DEGS_CD16mono %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD16mono)

write.csv(DEGS_CD16mono, "CD16 Mono cells Differentially expressed Genes Control C1 vs C2.CSV")   

#volcano plot CD16 Mono cells
plot3 <- EnhancedVolcano(DEGS_CD16mono, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD16mono$Gene, FCcutoff=0.5)
  ggsave(plot3,filename="Volcano Plot CD16 Mono cells control group C1 vs C2.png", height = 7.5, width = 7.5)

#CD4 TCM cells
CD4TCM.Control <- subset(x=Control, subset=predicted.celltype.l2==c("CD4 TCM"))
Idents(CD4TCM.Control) <- "Timepoint"
DefaultAssay(CD4TCM.Control) <- "RNA"
DEGS_CD4TCM <- FindMarkers(CD4TCM.Control, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4TCM <- DEGS_CD4TCM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4TCM)

write.csv(DEGS_CD4TCM, "CD4 TCM cells Differentially expressed Genes Control C1 vs C2.CSV")   

#volcano plot CD4 TCM cells
plot3 <- EnhancedVolcano(DEGS_CD4TCM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4TCM$Gene, FCcutoff=0.5)
  ggsave(plot3,filename="Volcano Plot CD4 TCM cells control group C1 vs C2.png", height = 7.5, width = 7.5)


#CD8 TEM cells
CD8TEM.Control <- subset(x=Control, subset=predicted.celltype.l2==c("CD8 TEM"))
Idents(CD8TEM.Control) <- "Timepoint"
DefaultAssay(CD8TEM.Control) <- "RNA"
DEGS_CD8TEM <- FindMarkers(CD8TEM.Control, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8TEM <- DEGS_CD8TEM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD8TEM)

write.csv(DEGS_CD8TEM, "CD8 TEM cells Differentially expressed Genes Control C1 vs C2.CSV")   

#volcano plot CD8 TEM cells
plot4 <- EnhancedVolcano(DEGS_CD8TEM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8TEM$Gene, FCcutoff=0.5)
  ggsave(plot4,filename="Volcano Plot CD8 TEM cells control group C1 vs C2.png", height = 7.5, width = 7.5)