#Plotting out initial cell clusters

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

#all cell UMAP
plot4 <- DimPlot(LBtox.no142, reduction="umap", group.by= "predicted.celltype.l2", label=TRUE, repel=TRUE)+ NoLegend()
ggsave(plot4, filename="LBtox no 142 whole UMAP.png", width = 6.5, height = 6.5)

#Cohort baseline no 142
cohort.baseline.no142 <- subset (x=Baseline.no142, subset=Population==c("Cohort"))

C2.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("C2"))
Idents(C2.no142) <- "predicted.celltype.l2"
table(C2.no142)

#NK Cells
Cohort.NK <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("NK"))
Idents(Cohort.NK) <- "Timepoint"
p1 <- DimPlot(Cohort.NK, group.by = "Timepoint", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p1, filename="NK cells clustered by Timepoint Legend.png", width = 8.5, height = 6.5)

Control.NK <- subset(x=Control.no142, subset=predicted.celltype.l2==c("NK"))

Baseline.NK <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("NK"))
Idents(Baseline.NK) <- "Population"
p2 <- DimPlot(Baseline.NK, group.by = "Population", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p2, filename="NK cells at baseline clustered by population.png", width = 8.5, height = 6.5)

C2.NK <- subset(x=C2.no142, subset=predicted.celltype.l2==c("NK"))
Idents(C2.NK) <- "Population"
p3 <- DimPlot(C2.NK, group.by = "Population", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p3, filename="NK cells at C2 clustered by population.png", width = 8.5, height = 6.5)

#CD4 TCM cells
Cohort.CD4TCM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD4 TCM"))
Idents(Cohort.CD4TCM) <- "Timepoint"
p4 <- DimPlot(Cohort.CD4TCM, group.by = "Timepoint", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p4, filename="CD4 TCM cells clustered by Timepoint Legend.png", width = 8.5, height = 6.5)

Baseline.CD4TCM <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD4 TCM"))
Idents(Baseline.CD4TCM) <- "Population"
p5 <- DimPlot(Baseline.CD4TCM, group.by = "Population", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p5, filename="CD4 TCM cells at baseline clustered by population.png", width = 8.5, height = 6.5)

C2.CD4TCM <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD4 TCM"))
Idents(C2.CD4TCM) <- "Population"
p6 <- DimPlot(C2.CD4TCM, group.by = "Population", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p6, filename="CD4 TCM cells at C2 clustered by population.png", width = 8.5, height = 6.5)

#CD14 Monocytes
Cohort.CD14mono <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(Cohort.CD14mono) <- "Timepoint"
p7 <- DimPlot(Cohort.CD14mono, group.by = "Timepoint", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p7, filename="CD14 Mono cells clustered by Timepoint Legend.png", width = 8.5, height = 6.5)

Baseline.CD14mono <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(Baseline.CD14mono) <- "Population"
p8 <- DimPlot(Baseline.CD14mono, group.by = "Population", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p8, filename="CD14 Mono cells at baseline clustered by population.png", width = 8.5, height = 6.5)

C2.CD14mono <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(C2.CD14mono) <- "Population"
p9 <- DimPlot(C2.CD14mono, group.by = "Population", label = TRUE, label.size = 3, repel=TRUE)
ggsave(p9, filename="CD14 Mono cells at C2 clustered by population.png", width = 8.5, height = 6.5)
p1 <- FeaturePlot(C2.CD14mono, features=c("NEBL", "HLA-DQA2", "S100A8", "LGALS2", "ADK", "VCAN"), split.by="Population")
ggsave(p1, filename="CD14 Monocytes Feature Plot Mapping Genes.png", height=15, width=7.5)

#CD16 Monocytes
CD16mono <- subset(x=LBtox.no142, subset=predicted.celltype.l2==c("CD16 Mono"))
Idents(CD16mono) <- "Population"
p1 <- DimPlot(CD16mono, group.by="Population", label=TRUE, label.size = 3, repel=TRUE)
ggsave(p1, filename="CD16 monocytes clustered by cohort vs control.png", width=8.5, height=6.5)

Cohort.CD16mono <- subset(x=CD16mono, subset=Population==c("Cohort"))
Idents(Cohort.CD16mono) <- "Timepoint"
p2 <- DimPlot(Cohort.CD16mono, group.by="Timepoint", label=TRUE, label.size = 3, repel=TRUE)
ggsave(p2, filename="CD16 monocytes cohort cells clustered by timepoint.png", width=8.5, height=6.5)

Baseline.CD16mono <- subset(x=CD16mono, subset=Timepoint==c("Baseline"))
Idents(Baseline.CD16mono) <- "Population"
p3 <- DimPlot(Baseline.CD16mono, group.by="Population", label=TRUE, label.size = 3, repel=TRUE)
ggsave(p3, filename="CD16 monocytes baseline cells clustered by population.png", width=8.5, height=6.5)

C2.CD16mono <- subset(x=CD16mono, subset=Timepoint==c("C2"))
Idents(C2.CD16mono) <- "Population"
p4 <- DimPlot(C2.CD16mono, group.by="Population", label=TRUE, label.size = 3, repel=TRUE)
ggsave(p4, filename="CD16 monocytes C2 cells clustered by population.png", width=8.5, height=6.5)

#CD8 TEM cells
CD8TEM <- subset(x=LBtox.no142, subset=predicted.celltype.l2==c("CD8 TEM"))
Idents(CD8TEM) <- "Population"

Cohort.CD8TEM <- subset(x=CD8TEM, subset=Population==c("Cohort"))
Idents(Cohort.CD8TEM) <- "Timepoint"
p5 <- DimPlot(Cohort.CD8TEM, group.by="Timepoint", label=TRUE, label.size=3, repel=TRUE)
ggsave(p5, filename="CD8 TEM cohort cells clustered by timepoint.png", width=8.5, height=6.5)

Baseline.CD8TEM <- subset(x=CD8TEM, subset=Timepoint==c("Baseline"))
Idents(Baseline.CD8TEM) <- "Population"
p6 <- DimPlot(Baseline.CD8TEM, group.by="Population", label=TRUE, label.size=3, repel=TRUE)
ggsave(p6, filename="CD8 TEM baseline cells clustered by population.png", width=8.5, height=6.5)

C2.CD8TEM <- subset(x=CD8TEM, subset=Timepoint==c("C2"))
Idents(C2.CD8TEM) <- "Population"
p6 <- DimPlot(C2.CD8TEM, group.by="Population", label=TRUE, label.size=3, repel=TRUE)
ggsave(p6, filename="CD8 TEM C2 cells clustered by population.png", width=8.5, height=6.5)