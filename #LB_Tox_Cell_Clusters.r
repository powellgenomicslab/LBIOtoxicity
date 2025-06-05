#LB_Tox_Cell_Clusters

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(jsonlite)
library(data.table)
library(devtools)
library(clustree)
library(Azimuth)

#install all packages
if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

#install scPred
devtools::install_github("immunogenomics/harmony", force=TRUE)
devtools::install_github("powellgenomicslab/scPred", force=TRUE)
library(scPred)

#install scTransform
devtools::install_github("satijalab/sctransform", ref = "develop")

#set working directory
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)
getwd()

#Read RDS
pbmc.tox <- LBtox
pbmc.tox <- readRDS(file = "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_metadata1.RDS")
pbmc.tox

#normalize data
pbmc.tox <- NormalizeData(pbmc.tox)

#scaling data
pbmc.tox <- ScaleData(pbmc.tox)

head(pbmc.tox)

#Azimuth Image
pbmc.tox.azi <- RunAzimuth(pbmc.tox, reference = "pbmcref")
pbmc.tox.azi
saveRDS(pbmc.tox.azi, file= "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/pbmctoxazimuth2.RDS")
getwd()

p1 <- DimPlot(pbmc.tox.azi, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmc.tox.azi, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
azimuth_classification <- p1 
plot1 <- azimuth_classification <- p1
ggsave(plot1, filename="LBTOX_combined_azimuth2.png",width = 8.5, height = 6.5)

setwd("/directflow/SCCGGroupShare/projects/jenli3/projects/LBTOX/azimuth_cell_class")

ggsave(plot1, filename="azimuth_classification_pbmctox.png",width = 10, height = 4.5)

#set idents to cell types
Idents(pbmc.tox.azi) <- 'predicted.celltype.l2'
levels(mostpbmc)

#Find Markers for Each Cluster
pbmc.markers <- FindAllMarkers(pbmc.tox.azi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

cluster1.markers <- FindMarkers(pbmc.tox.azi, ident.1 = "B intermediate", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

head(pbmc.markers)

#Feature Plot for Markers
plot1 <- FeaturePlot(pbmc.tox.azi, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A", "LEF1", "FOXP3"))
ggsave(plot1, filename="LBTOX_FeaturePlot_Azimuth.png",width = 10, height = 6)

getwd()

saveRDS(pbmc.markers, "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Markers.RDS")

#how many cells are in each cluster
table(Idents(pbmc.tox.azi))

#proportion of cells in each cluster?
prop.table(table(Idents(pbmc.tox.azi)))

#how does cluster membership vary by population?
table(Idents(pbmc.tox.azi), pbmc.tox.azi$Population)

#subset out a single indvidual and cluster
LB091 <- subset(x=pbmc.tox.azi, subset= MajoritySinglet_Individual_Assignment == c("LB091"))
LB091

#normalize
LB091 <- NormalizeData(LB091)

#scaling data
LB091 <- ScaleData(LB091)

#azimuth annotation
LB091.azi <- RunAzimuth(LB091, reference = "pbmcref")
LB091.azi

#UMAP clusters
p1 <- DimPlot(LB091.azi, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(LB091.azi, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
azimuth_classification <- p1 
plot1 <- azimuth_classification <- p1
ggsave(plot1, filename="LB091_combined_azimuth_classification_postscaling.png",width = 8.5, height = 6.5)

#set Idents
Idents(LB091.azi) <- 'predicted.celltype.l2'
levels(LB091.azi)

#summary of cells in LB091
table(Idents(LB091.azi))