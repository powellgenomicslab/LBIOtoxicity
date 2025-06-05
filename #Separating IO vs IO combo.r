#Separating IO vs IO combo

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

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_Merged_updatedmetadata_postQC.RDS")
LBtox

#Reminder of treatment type
colnames(LBtox[[]])
unique(LBtox$Treatment)

Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#Subset out Baseline cells
Baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))
Idents(Baseline) <- "predicted.celltype.l2"

#just cohort cells (from baseline)
Cohort.Baseline <- subset(x=Baseline, subset=Population==c("Cohort"))
Control.Baseline <-subset(x=Baseline, subset=Population==c("Control"))

#PD1 cohort baseline
PD1.Cohort.Baseline <- subset(x=Cohort.Baseline.no142, subset =Treatment==c("PD1"))
Combo.Cohort.Baseline <- subset(x=Cohort.Baseline.no142, subset=Treatment==c("IO_Combo"))

#UMAP Cohort Baseline cells - split by treatment
p1 <- DimPlot(Cohort.Baseline.no142, reduction="umap", group.by = "predicted.celltype.l2", label.size = 3, repel=TRUE, label=TRUE)
p2 <- DimPlot(Cohort.Baseline.no142, reduction="umap", group.by = "Treatment", label.size = 3, repel=TRUE, label=TRUE)
plot2 <- p1 + p2
ggsave(plot2, filename = "UMAP for Baseline Cohort Cells - PD1 vs IO_Combo, no LB142.png", height = 7.5, width = 15)

#cell component breakdown
Idents(PD1.Cohort.Baseline) <- "predicted.celltype.l2"
table(Idents(PD1.Cohort.Baseline))
prop.table(table(Idents(PD1.Cohort.Baseline)))*100

Idents(Combo.Cohort.Baseline) <- "predicted.celltype.l2"
table(Idents(Combo.Cohort.Baseline))
prop.table(table(Idents(Combo.Cohort.Baseline)))*100

#subset out C2
C2.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("C2"))
Idents(Baseline.no142) <- "predicted.celltype.l2"

#just cohort cells (from baseline)
C2.no142 <- subset(x=C2.no142, subset=Population==c("Cohort"))

#PD1 cohort baseline
PD1.Cohort.C2 <- subset(x=C2.no142, subset =Treatment==c("PD1"))
Combo.Cohort.C2 <- subset(x=C2.no142, subset=Treatment==c("IO_Combo"))

#cell component breakdown
Idents(PD1.Cohort.C2) <- "predicted.celltype.l2"
table(Idents(PD1.Cohort.C2))
prop.table(table(Idents(PD1.Cohort.C2)))*100

Idents(Combo.Cohort.C2) <- "predicted.celltype.l2"
table(Idents(Combo.Cohort.C2))
prop.table(table(Idents(Combo.Cohort.C2)))*100

#median statistics on cell numbers
summary()