#Making Ridge Plots for Expression of CTLA4 etc in specific cells

ibrary(Seurat)
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
install.packages("pheatmap")
library(pheatmap)
library(RColorBrewerlibrary)
library(viridis)

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_Merged_updatedmetadata_postQC.RDS")
LBtox
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#Fetch Data
Features_RP <- FetchData(LBtox, vars=c("CTLA4", "PDCD1"))

#subset population - baseline cohort cells
Baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))
Cohort.Baseline <- subset(x=Baseline, subset=Population==c("Cohort"))

#ridgeplot for Cohort Baseline cells and specific data

plot4 <- RidgePlot(Cohort.Baseline, features=Features_RP, ncol=2)
ggsave(plot4, filename="Ridgeplot Baseline Cells.png", height = 7.5, width = 7.5)

#subset baseline control cells
Control.Baseline <- subset(x=Baseline, subset=Population==c("Control"))

#Feature Plot for certain markers
plot1 <- FeaturePlot(LBtox, features=c("PDCD1", "JAK3", "CTLA4", "PTPN2", "FCGR2A", "ADAD1"))
ggsave(plot1, filename= "Feature Plot for Checkpoint Inhibitor Genes.png", height=10, width=7.5)