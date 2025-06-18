# R script: Seurat merge, QC, clustering, annotation, and batch effect assessment

# ==============================
# Load packages and dependencies
# ==============================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratDisk)
library(magrittr)
library(sctransform)
library(Azimuth)
library(jsonlite)
library(data.table)
library(devtools)
library(clustree)

# Install additional packages if needed
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

devtools::install_github("immunogenomics/harmony", force = TRUE)
devtools::install_github("powellgenomicslab/scPred", force = TRUE)
devtools::install_github("satijalab/sctransform", ref = "develop")
library(scPred)

# =======================
# Set working directory
# =======================
setwd("/directflow/SCCGGroupShare/projects/jenli3/LBTOX")

# =====================
# Basic QC and filtering
# =====================
LBtox <- readRDS("SeuratObjects/LBtox_Merged_NoQC_NoMetaData.RDS")
LBtox <- subset(LBtox, subset = nFeature_RNA > 350 & percent.mt < 15)
saveRDS(LBtox, "SeuratObjects/LBtox_Merged_updatedmetadata_postQC.RDS")

# Visual QC
vln_qc <- VlnPlot(LBtox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(vln_qc, filename = "Violin Plot QC Merged Pools.png", height = 4.5, width = 7.5)

# =====================
# SCTransform and PCA/UMAP
# =====================
LBtox <- readRDS("SeuratObjects/LBTox_Merged_Metadata_postQC.RDS")
LBtox <- SCTransform(LBtox, vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(LBtox, "SeuratObjects/LBtox_Merged_Metadata_postQC_postSCT.RDS")

LBtox <- RunPCA(LBtox, verbose = FALSE)
LBtox <- RunUMAP(LBtox, dims = 1:30, verbose = FALSE)
LBtox <- FindNeighbors(LBtox, dims = 1:30, verbose = FALSE)
LBtox <- FindClusters(LBtox, verbose = FALSE)

umap_plot <- DimPlot(LBtox, reduction = "umap", label = TRUE)
ggsave(umap_plot, filename = "Unlabelled UMAP POSTQC.png", width = 8.5, height = 6.5)

# =======================
# Azimuth Annotation
# =======================
LBtox <- RunAzimuth(LBtox, reference = "pbmcref")
saveRDS(LBtox, "SeuratObjects/LBtox_PostQC_Annotated.RDS")

p1 <- DimPlot(LBtox, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 <- DimPlot(LBtox, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
ggsave(p1, filename = "LBTOX_combinedPOSTQC_azimuth_repel.png", width = 8.5, height = 6.5)

cas_plot <- FeaturePlot(LBtox, features = "CASP3")
ggsave(cas_plot, filename = "Caspase Expression on General UMAP.png", width = 8.5, height = 6.5)

saveRDS(LBtox, file = "SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")

# =======================
# Batch Effect Assessment
# =======================
pbmc.tox1 <- readRDS("SeuratObjects/LBtoxCombined2.RDS")

# Split by pool
tox.pool <- SplitObject(pbmc.tox1, split.by = "Pool")
tox.pool <- lapply(tox.pool, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# Integration
features <- SelectIntegrationFeatures(object.list = tox.pool)
anchors <- FindIntegrationAnchors(object.list = tox.pool, anchor.features = features)
tox.integrated <- IntegrateData(anchorset = anchors)

# Dimensionality reduction
DefaultAssay(tox.integrated) <- "integrated"
tox.integrated <- ScaleData(tox.integrated, verbose = FALSE)
tox.integrated <- RunPCA(tox.integrated, features = features, verbose = FALSE)
tox.integrated <- RunUMAP(tox.integrated, reduction = "pca", dims = 1:30)
tox.integrated <- FindNeighbors(tox.integrated, reduction = "pca", dims = 1:30)
tox.integrated <- FindClusters(tox.integrated, resolution = 0.5)

# Batch visualization
plot_pool <- DimPlot(LBtox, reduction = "umap", group.by = "Pool")
ggsave(plot_pool, filename = "UMAP for Pool Batch Effect.png", height = 4.5, width = 4.5)

plot_individuals <- DimPlot(LBtox, reduction = "umap", group.by = "MajoritySinglet_Individual_Assignment")
ggsave(plot_individuals, filename = "UMAP for individuals.png", height = 4.5, width = 4.5)
