#Merge->QC->Analysis

#import packages
library(Seurat) 
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(sctransform)
library(jsonlite)
library(data.table)
library(Azimuth)

#set wd
#set working directory
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_Merged_NoQC_NoMetaData.RDS")
LBtox

#QC with new parameters as per Joseph
LBtox <- subset(LBtox, subset = nFeature_RNA>350 & percent.mt < 15)
LBtox

plot1 <- VlnPlot(LBtox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot1, filename="Violin Plot QC Merged Pools.png", height = 4.5, width = 7.5)

#save RDS
saveRDS(LBtox,file="/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_Merged_updatedmetadata_postQC.RDS" )

#PBMC workflow 
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_postQC.RDS")
LBtox <- SCTransform(LBtox, vars.to.regress="percent.mt", verbose=FALSE)

saveRDS(LBtox,file="/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_Merged_Metadata_postQC_postSCT.RDS" )

LBtox <- RunPCA(LBtox, verbose=FALSE)
LBtox <- RunUMAP(LBtox, dims=1:30, verbose=FALSE)
LBtox <- FindNeighbors(LBtox, dims= 1:30, verbose=FALSE)
LBtox <- FindClusters(LBtox, verbose=FALSE)
plot1 <- DimPlot(LBtox, reduction="umap", label=TRUE)
ggsave(plot1, filename="Unlabelled UMAP POSTQC.png", width = 8.5, height = 6.5)

#Azimuth Annotation
LBtox <- RunAzimuth(LBtox, reference = "pbmcref", )
LBtox
saveRDS(LBtox, file= "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_PostQC_Annotated.RDS")
getwd()

p1 <- DimPlot(LBtox, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel=TRUE) + NoLegend()
p2 <- DimPlot(LBtox, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
azimuth_classification <- p1 
plot3 <- azimuth_classification <- p1
ggsave(plot3, filename="LBTOX_combinedPOSTQC_azimuth_repel.png",width = 8.5, height = 6.5)

#Feature Plot - Caspase markers
plot4 <- FeaturePlot(LBtox, features=c("CASP3"))
ggsave(plot4, filename="Caspase Expression on General UMAP.png",width = 8.5, height = 6.5)