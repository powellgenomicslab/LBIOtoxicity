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

#Read 10x data

LBtox1 <- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/221125/IRAE_LB10_pool1/GE/summary/filtered_feature_bc_matrix")
LBtox2 <- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/221125/IRAE_LB10_pool2/GE/summary/filtered_feature_bc_matrix")
LBtox3 <- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/221125/IRAE_LB10_pool3/GE/summary/filtered_feature_bc_matrix")

#Read RDS
pbmc.tox <- readRDS(file = "/directflow/SCCGGroupShare/projects/jenli3/pbmctoxwholepostSCT.rds")
pbmc.tox

#make into seurat objects

pbmc.tox1 <- CreateSeuratObject(counts=LBtox1, project ="Toxicity", min.cells = 3, min.features=200)
pbmc.tox1
pbmc.tox2 <- CreateSeuratObject(counts=LBtox2, project ="Toxicity", min.cells = 3, min.features=200)
pbmc.tox2
pbmc.tox3 <- CreateSeuratObject(counts=LBtox3, project= "Toxicity", min.cells = 3, min.features=200)
pbmc.tox3

#integrate seurat objects
pbmc.tox <- merge(pbmc.tox1, y=c(pbmc.tox2, pbmc.tox3), add.cell.ids = c("Pool1", "Pool2", "Pool3"), project="CombinedToxPBMC")
pbmc.tox

# QC individual pools
pbmc.tox1[["percent.mt"]] <- PercentageFeatureSet(pbmc.tox1, pattern = "^MT-")
pbmc.tox2[["percent.mt"]] <- PercentageFeatureSet(pbmc.tox2, pattern = "^MT-")
pbmc.tox3[["percent.mt"]] <- PercentageFeatureSet(pbmc.tox3, pattern = "^MT-")

#QC whole pool
pbmc.tox[["percent.mt"]] <- PercentageFeatureSet(pbmc.tox, pattern = "^MT-")

#violin plots
plot1 <- VlnPlot(pbmc.tox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot1, filename = "TOX_VLN_wholepool.png", height = 4.5, width = 7.5)

plot2 <- VlnPlot(pbmc.tox2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(plot2, filename = "TOX2_VLN1.png", height = 4.5, width = 7.5)

plot3 <- VlnPlot(pbmc.tox3, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(plot3, filename = "TOX3_VLN1.png", height = 4.5, width = 7.5)


#filtering cells
pbmc.tox1 <- subset(pbmc.tox1, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & percent.mt < 12)
pbmc.tox1

pbmc.tox2 <- subset(pbmc.tox2, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & percent.mt < 12)
pbmc.tox2

pbmc.tox3 <- subset(pbmc.tox3, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & percent.mt < 12)
pbmc.tox3

pbmc.tox <- subset(pbmc.tox, subset= nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt <12)
pbmc.tox

#repeat VLN plots after filtering
plot4 <- VlnPlot(pbmc.tox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot4, filename = "LBtox_VLN_postfilter_WHOLEPOOL.png", height = 4.5, width = 7.5)

plot5 <- VlnPlot(pbmc.tox2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot5, filename = "LBtox2_VLN_postfilter.png", height = 4.5, width = 7.5)

plot6 <- VlnPlot(pbmc.tox3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot6, filename = "LBtox3_VLN_postfilter.png", height = 4.5, width = 7.5)

plot7 <- VlnPlot(pbmc.tox3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot6, filename = "LBtox3_VLN_postfilter.png", height = 4.5, width = 7.5)

#save Seurat Object
saveRDS(pbmc.tox1, file="/directflow/SCCGGroupShare/projects/jenli3/pbmctox1.rds")
saveRDS(pbmc.tox2, file="/directflow/SCCGGroupShare/projects/jenli3/pbmctox2.rds")
saveRDS(pbmc.tox3, file="/directflow/SCCGGroupShare/projects/jenli3/pbmctox3.rds")
saveRDS(pbmc.tox, file="/directflow/SCCGGroupShare/projects/jenli3/pbmctoxwhole.rds")

#scaling data
pbmc.tox1 <- ScaleData(pbmc.tox1)
pbmc.tox2 <- ScaleData(pbmc.tox2)
pbmc.tox3 <- ScaleData(pbmc.tox3)
pbmc.tox <- ScaleData(pbmc.tox)

#normalising pool 1 with scTransform
pbmc.tox <-SCTransform(pbmc.tox, vars.to.regress="percent.mt", verbose=FALSE)
pbmc.tox <- RunPCA(pbmc.tox, verbose = FALSE)
pbmc.tox <- RunUMAP(pbmc.tox, dims = 1:30, verbose = FALSE)

head(pbmc.tox)

#clustree to determine resolution
pbmc.tox[['TSNE']] <- CreateDimReducObject(embeddings = pbmc.tox,
                                         key = "tSNE_", assay="RNA")


#initial clustering and UMAP pool 1
pbmc.tox <- RunPCA(pbmc.tox, features = VariableFeatures(object = pbmc.tox))
pbmc.tox <- FindNeighbors(pbmc.tox, dims = 1:10)
pbmc.tox <- FindClusters(pbmc.tox, resolution = 1.5)
pbmc.tox <- RunUMAP(pbmc.tox, dims = 1:20)
plot8 <- DimPlot(pbmc.tox, reduction = "umap")
ggsave(plot8, filename = "LBtox_resolution1.5_UMAP_wholepool.png", height = 4.5, width = 7.5)

#azimuth tutorial
reference <- LoadH5Seurat("/directflow/SCCGGroupShare/projects/jenli3/mostpbmc/pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE, raster=FALSE) + NoLegend()

#find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc.tox,
  normalization.method = "SCT",
  reference.reduction = "spca",
)

#save RDS post processing 
saveRDS(pbmc.tox, file="/directflow/SCCGGroupShare/projects/jenli3/pbmctoxwholepostSCT.rds")

#transfer cell type labels and protein data
pbmc.tox <- MapQuery(anchorset = anchors,
  query = pbmc.tox,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap")

#visualise annotated cells pool 1
p1 = DimPlot(pbmc.tox, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc.tox, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
plot8 <- p1 + p2
ggsave(plot8, filename = "LBtox_annotated_UMAP_WholePool.png", height = 4.5, width = 7.5)

#Find Markers

pbmc.markers <- FindAllMarkers(pbmc.tox, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


cluster0.markers <- FindMarkers(pbmc.tox, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

plot9 <- FeaturePlot(pbmc.tox, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
    ggsave(plot9, filename = "FeaturesperCluster.png", height = 4.5, width = 7.5)

#normalising pool 2 with scTransform
pbmc.tox2 <-SCTransform(pbmc.tox2, vars.to.regress="percent.mt", verbose=FALSE)
pbmc.tox2 <- RunPCA(pbmc.tox2, verbose = FALSE)
pbmc.tox2 <- RunUMAP(pbmc.tox2, dims = 1:30, verbose = FALSE)

#initial clustering and UMAP pool 2
pbmc.tox2 <- RunPCA(pbmc.tox2, features = VariableFeatures(object = pbmc.tox2))
pbmc.tox2 <- FindNeighbors(pbmc.tox2, dims = 1:10)
pbmc.tox2 <- FindClusters(pbmc.tox2, resolution = 0.5)
pbmc.tox2 <- RunUMAP(pbmc.tox2, dims = 1:20)
plot9 <- DimPlot(pbmc.tox2, reduction = "umap")
ggsave(plot9, filename = "LBtox2_initial_UMAP.png", height = 4.5, width = 7.5)

#find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc.tox2,
  normalization.method = "SCT",
  reference.reduction = "spca",
)

#transfer cell type labels and protein data
pbmc.tox2 <- MapQuery(anchorset = anchors,
  query = pbmc.tox2,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap")

#visualise annotated cells
p1 = DimPlot(pbmc.tox2, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc.tox2, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
plot10 <- p1 + p2
ggsave(plot10, filename = "LBtox2_annotated_UMAP.png", height = 4.5, width = 7.5)

##Pool 3 clusters

#normalising pool 3 with scTransform
pbmc.tox3 <-SCTransform(pbmc.tox3, vars.to.regress="percent.mt", verbose=FALSE)
pbmc.tox3 <- RunPCA(pbmc.tox3, verbose = FALSE)
pbmc.tox3 <- RunUMAP(pbmc.tox3, dims = 1:30, verbose = FALSE)

#initial clustering and UMAP pool 3
pbmc.tox3 <- RunPCA(pbmc.tox3, features = VariableFeatures(object = pbmc.tox3))
pbmc.tox3 <- FindNeighbors(pbmc.tox3, dims = 1:10)
pbmc.tox3 <- FindClusters(pbmc.tox3, resolution = 0.5)
pbmc.tox3 <- RunUMAP(pbmc.tox3, dims = 1:20)
plot13 <- DimPlot(pbmc.tox3, reduction = "umap")
ggsave(plot13, filename = "LBtox3_initial_UMAP.png", height = 4.5, width = 7.5)

#find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc.tox3,
  normalization.method = "SCT",
  reference.reduction = "spca",
)

#transfer cell type labels and protein data
pbmc.tox3 <- MapQuery(anchorset = anchors,
  query = pbmc.tox3,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap")

#visualise annotated cells
p1 = DimPlot(pbmc.tox3, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc.tox3, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
plot14 <- p1 + p2
ggsave(plot14, filename = "LBtox3_annotated_UMAP.png", height = 4.5, width = 7.5)

#scPred annotation
reference <- scPred::pbmc_1
query <- scPred::pbmc_2

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(reference, dims = 1:30)

  DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

#training classifiers with scPred
reference <- getFeatureSpace(reference, "cell_type")
library(getFeatureSpace)
install.packages("getFeatureSpace")

#umap scpred reference plot
plot5 <- DimPlot(scpredref, group.by = "cell_type", label = TRUE, repel = TRUE)
ggsave(plot5, filename = "scPredref_UMAP.png", height = 4.5, width = 7.5)

scpredref <- getFeatureSpace(scpredref, "cell_type")
scpredref <- trainModel(scpredref)


#plot scPred predictions
plot6 <- DimPlot(mostpbmc1, group.by = "scpred_prediction", reduction = "scpred")
ggsave(plot6, filename = "scPred_predictions1.png", height = 4.5, width = 7.5)

#UMAP with scPred aligned data
mostpbmc1 <- RunUMAP(mostpbmc1, reduction = "scpred", dims = 1:30)
plot7 <- DimPlot(mostpbmc1, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
ggsave(plot7, filename = "scPred_UMAP_annotated1.png", height = 4.5, width = 7.5)
