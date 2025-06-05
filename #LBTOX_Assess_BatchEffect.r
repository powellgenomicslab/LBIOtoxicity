#LBTOX_Assess_BatchEffect

#install packages
library(Seurat) 
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(jsonlite)
library(data.table)
library(Azimuth)
library(sctransform)

#set wd
#set working directory
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#Read RDS
pbmc.tox1 <- readRDS(file = "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtoxCombined2.RDS")
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_PostQC_Annotated.RDS")

colnames(LBtox[[]])
unique(pbmc.tox$Pool)

#split dataset into three objects (pools)
tox.pool <-SplitObject(pbmc.tox1, split.by="Pool")
pool1 <- tox.pools[["LBTOX_Pool1"]]
pool2 <- tox.pools [["LBTox_Pool2"]]
pool3 <- tox.pools [["LBTox_Pool3"]]

tox.pool

# normalize and identify variable features for each dataset independently w scTransform

tox.pool <- lapply(X = tox.pool, FUN= function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = tox.pool)

#perform integration
immune.anchors <- FindIntegrationAnchors(object.list= tox.pool, anchor.features=features)

#create integrated data assay
tox.pools.combined <- IntegrateData(anchorset=immune.anchors)

#perform an integrated analysis
DefaultAssay(tox.pools.combined) <- "integrated"

tox.pools.combined

#run standard workflow - already normalized and scaled 
tox.pools.combined <- ScaleData(tox.pools.combined, verbose = FALSE)
tox.pools.combined <- RunPCA(tox.pools.combined, features=features, verbose =FALSE)
tox.pools.combined <- RunUMAP(tox.pools.combined, reduction = "pca", dims = 1:30)
tox.pools.combined <- FindNeighbors(tox.pools.combined, reduction = "pca", dims = 1:30)
tox.pools.combined <- FindClusters(tox.pools.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(Baseline, reduction = "umap", group.by = "Population")
p2 <- DimPlot(Baseline, reduction = "umap", label = TRUE, repel = TRUE)
plot5 <- p1 + p2
ggsave(plot5, filename = "UMAP for Cohort vs Control Baseline.png", height = 4.5, width = 9.5)

#assess for batch effects
plot5 <- DimPlot(LBtox, reduction = "umap", group.by = "Pool")
ggsave(plot5, filename = "UMAP for Pool Batch Effect.png", height = 4.5, width = 4.5)

plot6 <- DimPlot(LBtox, reduction = "umap", group.by = "MajoritySinglet_Individual_Assignment")
ggsave(plot6, filename = "UMAP for individuals.png", height = 4.5, width = 4.5)