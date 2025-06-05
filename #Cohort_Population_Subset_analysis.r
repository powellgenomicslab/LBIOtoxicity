#Cohort_Population_Subset_analysis
#subset out a single indvidual and cluster
Baseline <- subset(x=LBtox, subset= Timepoint == c("Baseline"))
Baseline

unique(LBtox$predicted.celltype.l2)

CD4Baseline <- subset(x=Baseline, subset = predicted.celltype.l2 == c("CD4 TEM"))


#normalize
Baseline <- NormalizeData(Baseline)

#scaling data
Baseline <- ScaleData(Baseline)

#azimuth annotation
Baseline.azi <- RunAzimuth(Baseline, reference = "pbmcref")
Baseline.azi

#UMAP clusters
p1 <- DimPlot(cohort.azi, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(cohort.azi, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
azimuth_classification <- p1 
plot2 <- azimuth_classification <- p1
ggsave(plot2, filename="Cohort_Azimuth_Classification.png",width = 8.5, height = 6.5)

#set Idents
Idents(CD4Baseline) <- 'Population'
levels(Baseline)

#find Markers for clusters in cohort population
cohort.azi.markers <- FindAllMarkers(cohort.azi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cohort.azi.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

cluster1.markers <- FindMarkers(pbmc.tox.azi, ident.1 = "B intermediate", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#find markers for baseline CD4 TEM cells between cohort and controls
CD4.baseline.markers <- FindAllMarkers(CD4Baseline, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
baseline.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

CD4.baseline.markers


#summary of cells in LB091
table(Idents(cohort.azi))

#DGE between Treg cells and CD8 TEM cells
reg.tem.de.markers <- FindMarkers(cohort.azi, ident.1="CD8 TEM", ident.2 = "Treg")

#view results
head(reg.tem.de.markers)

#split by cohort/control - UMAP
Read object
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/pbmc_toxicity_completemetadata.RDS")
cohort.pool <-SplitObject(LBtox, split.by="Population")

cohort.pool

#normalize
cohort.pool <- lapply(X = cohort.pool, FUN= function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

features <- SelectIntegrationFeatures(object.list = cohort.pool)

#perform integration
immune.anchors <- FindIntegrationAnchors(object.list= cohort.pool, anchor.features=features)

#create integrated data assay
cohort.pool <- IntegrateData(anchorset=immune.anchors)

#perform an integrated analysis
DefaultAssay(cohort.pool) <- "integrated"

#run standard workflow - already normalized and scaled 
cohort.pool <- ScaleData(cohort.pool, verbose = FALSE)
cohort.pool <- RunPCA(cohort.pool, features=features, verbose =FALSE)
cohort.pool <- RunUMAP(cohort.pool, reduction = "pca", dims = 1:30)
cohort.pool <- FindNeighbors(cohort.pool, reduction = "pca", dims = 1:30)
cohort.pool <- FindClusters(cohort.pool, resolution = 0.5)

# Visualization
p1 <- DimPlot(cohort.pool, reduction = "umap", group.by = "Population")
p2 <- DimPlot(cohort.pool, reduction = "umap", label = TRUE, repel = TRUE)
plot3 <- p1 + p2
ggsave(plot3, filename = "UMAP for Cohort vs Control.png", height = 4.5, width = 7.5)

#annotation
pbmc.tox.azi <- RunAzimuth(cohort.pool, reference = "pbmcref")

p1 <- DimPlot(pbmc.tox.azi, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmc.tox.azi, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
azimuth_classification <- p1 
plot1 <- azimuth_classification <- p1
ggsave(plot1, filename="cohort.pool.azi.png",width = 8.5, height = 6.5)