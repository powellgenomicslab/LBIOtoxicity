#DGE between timepoints 
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
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Cohort.no142 <- subset(x=LBtox.no142, subset=Population==c("Cohort"))
Control.no142 <- subset(x=LBtox.no142, subset=Population==c("Control"))
Idents(Baseline.no142) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#Cohort baseline no 142
cohort.baseline.no142 <- subset (x=Baseline.no142, subset=Population==c("Cohort"))

#C2 no LB142
C2.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("C2"))
Idents(C2.no142) <- "predicted.celltype.l2"
table(C2.no142)

#control group no 142
control.baseline.no142 <- subset(x=Baseline.no142, subset=Population==c("Control"))
Idents(Control.no142) <- "Timepoint"


#Cohort C2
cohort.C2.no142 <- subset(x= C2.no142, subset=Population==c("Cohort"))
Idents(cohort.C2.no142) <- "predicted.celltype.l2"
table(Idents(cohort.C2.no142))

#Control C2
control.C2.no142 <- subset(x= C2.no142, subset=Population==c("Control"))
Idents(control.C2.no142) <- "predicted.celltype.l2"
table(Idents(control.C2.no142))

#subset out individual cell types and compare C1/C2
#CD4 TCM Cells
Cohort.CD4.TCM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD4 TCM"))
Idents(Cohort.CD4.TCM) <- "Timepoint"

DefaultAssay(Cohort.CD4.TCM) <- "RNA"
DEGS_CD4_TCM <- FindMarkers(Cohort.CD4.TCM, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4_TCM <- DEGS_CD4_TCM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4_TCM)

write.csv(DEGS_CD4_TCM, "CD4 TCM Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot
plot1 <- EnhancedVolcano(DEGS_CD4_TCM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4_TCM$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot CD4 TCM cells C1 vs C2.png", height = 7.5, width = 7.5)

#B intermediate cells
Cohort.Bint <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("B intermediate"))
Idents(Cohort.Bint) <- "Timepoint"

DefaultAssay(Cohort.Bint) <- "RNA"
DEGS_Bint <- FindMarkers(Cohort.Bint, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_Bint <- DEGS_Bint %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_Bint)

write.csv(DEGS_Bint, "B intermediate cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot B intermediate cells
plot1 <- EnhancedVolcano(DEGS_Bint, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bint$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot B intermediate cells C1 vs C2.png", height = 7.5, width = 7.5)

#CD8 TEM Cells
Cohort.CD8TEM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD8 TEM"))
Idents(Cohort.CD8TEM) <- "Timepoint"

DefaultAssay(Cohort.CD8TEM) <- "RNA"
DEGS_CD8TEM <- FindMarkers(Cohort.CD8TEM, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8TEM <- DEGS_CD8TEM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD8TEM)

write.csv(DEGS_CD8TEM, "CD8 TEM cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot B intermediate cells
plot1 <- EnhancedVolcano(DEGS_CD8TEM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8TEM$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot CD8 TEM cells C1 vs C2.png", height = 7.5, width = 7.5)

#NK Cells
Cohort.NK <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("NK"))
Idents(Cohort.NK) <- "Timepoint"
DefaultAssay(Cohort.NK) <- "RNA"
DEGS_NK <- FindMarkers(Cohort.NK, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_NK <- DEGS_NK %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_NK)

write.csv(DEGS_NK, "NK cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot NK cells
plot2 <- EnhancedVolcano(DEGS_NK, x="avg_log2FC", y="p_val_adj", lab=DEGS_NK$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot NK cells C1 vs C2.png", height = 7.5, width = 7.5)

#CD4 TEM
Cohort.CD4TEM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD4 TEM"))
Idents(Cohort.CD4TEM) <- "Timepoint"
DefaultAssay(Cohort.CD4TEM) <- "RNA"
DEGS_CD4TEM <- FindMarkers(Cohort.CD4TEM, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4TEM <- DEGS_CD4TEM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4TEM)

write.csv(DEGS_CD4TEM, "CD4 TEM cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot CD4 TEM cells
plot2 <- EnhancedVolcano(DEGS_CD4TEM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4TEM$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot CD4 TEM cells C1 vs C2.png", height = 7.5, width = 7.5)

#CD4 naive
Cohort.CD4naive <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD4 Naive"))
Idents(Cohort.CD4naive) <- "Timepoint"
DefaultAssay(Cohort.CD4naive) <- "RNA"
DEGS_CD4naive <- FindMarkers(Cohort.CD4naive, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4naive <- DEGS_CD4naive %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4naive)

write.csv(DEGS_CD4naive, "CD4 naive cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot CD4 naive cells
plot3 <- EnhancedVolcano(DEGS_CD4naive, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4naive$Gene, FCcutoff=0.1)
  ggsave(plot3,filename="Volcano Plot CD4 Naive cells C1 vs C2.png", height = 7.5, width = 7.5)

#B naive
Cohort.Bnaive <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("B naive"))
Idents(Cohort.Bnaive) <- "Timepoint"
DefaultAssay(Cohort.Bnaive) <- "RNA"
DEGS_Bnaive <- FindMarkers(Cohort.Bnaive, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_Bnaive <- DEGS_Bnaive %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_Bnaive)

write.csv(DEGS_Bnaive, "B naive cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot B naive cells
plot4 <- EnhancedVolcano(DEGS_Bnaive, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bnaive$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot B Naive cells C1 vs C2.png", height = 7.5, width = 7.5)

#CD14 monocytes
Cohort.CD14mono <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(Cohort.CD14mono) <- "Timepoint"
DefaultAssay(Cohort.CD14mono) <- "RNA"
DEGS_CD14Mono <- FindMarkers(Cohort.CD14mono, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD14Mono <- DEGS_CD14Mono %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD14Mono)

write.csv(DEGS_CD14Mono, "CD14 Monocyte cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot CD14 Mono cells
plot6 <- EnhancedVolcano(DEGS_CD14Mono, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD14Mono$Gene, FCcutoff=0.5)
  ggsave(plot6,filename="Volcano Plot CD14 Monocyte cells C1 vs C2.png", height = 7.5, width = 7.5)

  #CD16 monocytes
Cohort.CD16mono <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD16 Mono"))
Idents(Cohort.CD16mono) <- "Timepoint"
DefaultAssay(Cohort.CD16mono) <- "RNA"
DEGS_CD16Mono <- FindMarkers(Cohort.CD16mono, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD16Mono <- DEGS_CD16Mono %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD16Mono)

write.csv(DEGS_CD16Mono, "CD16 Monocyte cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot CD16 Mono cells
plot7 <- EnhancedVolcano(DEGS_CD16Mono, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD16Mono$Gene, FCcutoff=0.1)
  ggsave(plot7,filename="Volcano Plot CD16 Monocyte cells C1 vs C2.png", height = 7.5, width = 7.5)

#CD8 TCM
Cohort.CD8TCM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD8 TCM"))
Idents(Cohort.CD8TCM) <- "Timepoint"
DefaultAssay(Cohort.CD8TCM) <- "RNA"
DEGS_CD8TCM <- FindMarkers(Cohort.CD8TCM, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8TCM <- DEGS_CD8TCM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD8TCM)

write.csv(DEGS_CD8TCM, "CD8 TCM cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot CD8 TCM cells
plot8 <- EnhancedVolcano(DEGS_CD8TCM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8TCM$Gene, FCcutoff=0.1)
  ggsave(plot8,filename="Volcano Plot CD8 TCM cells C1 vs C2.png", height = 7.5, width = 7.5)

#CD8 Naive
Cohort.CD8naive <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD8 Naive"))
Idents(Cohort.CD8naive) <- "Timepoint"
DefaultAssay(Cohort.CD8naive) <- "RNA"
DEGS_CD8naive <- FindMarkers(Cohort.CD8naive, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8naive <- DEGS_CD8naive %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD8naive)

write.csv(DEGS_CD8naive, "CD8 Naive cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot CD8 Naive cells
plot9 <- EnhancedVolcano(DEGS_CD8naive, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8naive$Gene, FCcutoff=0.1)
  ggsave(plot9,filename="Volcano Plot CD8 Naive cells C1 vs C2.png", height = 7.5, width = 7.5)

#Treg
Cohort.Treg <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("Treg"))
Idents(Cohort.Treg) <- "Timepoint"
DefaultAssay(Cohort.Treg) <- "RNA"
DEGS_Treg <- FindMarkers(Cohort.Treg, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_Treg <- DEGS_Treg %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_Treg)

write.csv(DEGS_Treg, "T reg cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot Treg cells
plot1 <- EnhancedVolcano(DEGS_Treg, x="avg_log2FC", y="p_val_adj", lab=DEGS_Treg$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot T reg cells C1 vs C2.png", height = 7.5, width = 7.5)

#gDT
Cohort.gdT <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("gdT"))
Idents(Cohort.gdT) <- "Timepoint"
DefaultAssay(Cohort.gdT) <- "RNA"
DEGS_gdT <- FindMarkers(Cohort.gdT, ident.1="Baseline", ident.2="C2", vars.to.regress=c("Pool", "Sex"))
DEGS_gdT <- DEGS_gdT %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_gdT)

write.csv(DEGS_gdT, "gdT cells Differentially expressed Genes Cohort C1 vs C2.CSV")   

#volcano plot gdT cells
plot2 <- EnhancedVolcano(DEGS_gdT, x="avg_log2FC", y="p_val_adj", lab=DEGS_gdT$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot gdT cells C1 vs C2.png", height = 7.5, width = 7.5)