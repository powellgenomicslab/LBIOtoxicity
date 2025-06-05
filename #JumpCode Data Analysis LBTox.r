#JumpCode Data Analysis


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
getwd()

#Read in raw data
LBtox3 <- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/jenli3//R_and_D_JumpCode_results/analyses/IRAE_LB10_pool3_JumpCode/GE_COUNT/summary/filtered_feature_bc_matrix/")

#create Seurat Objects for individual pools
LBtox3 <- CreateSeuratObject(counts = LBtox3, min.cells=3, min.features=200)

#read in demultiplexing and doublet detection results
demuxafy.tox.1 <- read.table("/directflow/SCCGGroupShare/projects/himaro/demuxafy/toxicity_pbmc/output/combined_results/IRAE_LB10_pool1/combined_results_w_combined_assignments.tsv", sep = "\t", header=TRUE)
rownames(demuxafy.tox.1) <- demuxafy.tox.1$Barcode

demuxafy.tox.2 <- read.table("/directflow/SCCGGroupShare/projects/himaro/demuxafy/toxicity_pbmc/output/combined_results/IRAE_LB10_pool2/combined_results_w_combined_assignments.tsv", sep = "\t", header=TRUE)
rownames(demuxafy.tox.2) <- demuxafy.tox.2$Barcode

demuxafy.tox.3 <- read.table("/directflow/SCCGGroupShare/projects/himaro/demuxafy/toxicity_pbmc/output/combined_results/IRAE_LB10_pool3/combined_results_w_combined_assignments.tsv", sep = "\t", header=TRUE)
rownames(demuxafy.tox.3) <- demuxafy.tox.3$Barcode

#Add demuxafy data to each pool object
LB.tox.1 <- AddMetaData (LB.tox.1, demuxafy.tox.1)
LB.tox.2 <- AddMetaData (LB.tox.2, demuxafy.tox.2)
LBtox3 <- AddMetaData (LBtox3, demuxafy.tox.3)

LBtox3

#saveRDS
saveRDS(LBtox3, filename="/directflow/SCCGGroupShare/projects/jenli3/LBTOX/JumpCode Analysis/LBtox3_Jumpcode.RDS")

#QC
#mitochondrial counts
LBtox3[["percent.mt"]] <- PercentageFeatureSet(LBtox3, pattern = "^MT-")

#violin plots
plot1 <- VlnPlot(LBtox3, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(plot1, filename = "LBtox3_VLN.png", height = 4.5, width = 7.5)

#metrics for each count
summary(LBtox3$nFeature_RNA)
summary(LBtox3$nCount_RNA)
summary(LBtox3$percent.mt)

#cluster pre-filtering
LBtox3 <- SCTransform(LBtox3, vars.to.regress="percent.mt", verbose=FALSE)
LBtox3
LBtox3 <- RunPCA(LBtox3, verbose=FALSE)
LBtox3 <- RunUMAP(LBtox3, dims=1:30, verbose=FALSE)
LBtox3 <- FindNeighbors(LBtox3, dims= 1:30, verbose=FALSE)
LBtox3 <- FindClusters(LBtox3, verbose=FALSE)
plot2 <- DimPlot(LBtox3, reduction="umap", label=TRUE)
ggsave(plot2, filename="Unlabelled JumpCode UMAP pre QC.png", width = 8.5, height = 6.5)

#save RDS
saveRDS(LBtox3,file="/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/JumpCode_Data_preQC_postSCT.RDS" )

#Azimuth Annotation
LBtox3.preQC <- RunAzimuth(LBtox3, reference = "pbmcref")
LBtox3.preQC
saveRDS(LBtox3.preQC, file= "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/JumpCode Analysis/LBtox_PreQC_Annotated.RDS")
getwd()

p1 <- DimPlot(LBtox3.preQC, group.by = "predicted.celltype.l2", label = TRUE, repel=TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(LBtox3.preQC, group.by = "prediction.score.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
azimuth_classification <- p1 
plot3 <- azimuth_classification <- p1
ggsave(plot3, filename="JumpCode_PREQC_annotated.png",width = 8.5, height = 6.5)

#readRDS

JumpCode <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/JumpCode Analysis/LBtox_PreQC_Annotated.RDS")

#filtering
#QC with new parameters as per Joseph
LBtox3 <- subset(JumpCode, subset = nFeature_RNA>350 & percent.mt < 15)
LBtox3
Idents(object=LBtox3)
levels(x=LBtox3)

plot4 <- VlnPlot(LBtox3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL, ncol = 3)
ggsave(plot4, filename="Violin Plot post QC JumpCode Data.png", height = 4.5, width = 7.5)

#save RDS
saveRDS(LBtox,file="/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/Jumpcode_Data_postQC_postSCT.RDS" )

#cluster post QC
plot5 <- DimPlot(LBtox3, group.by = "predicted.celltype.l2", label = TRUE, repel=TRUE, label.size = 3) + NoLegend()
ggsave(plot5, filename="JumpCode_postQC_annotated.png",width = 8.5, height = 6.5)

#metrics post QC
summary(LBtox3$nFeature_RNA)
summary(LBtox3$nCount_RNA)
summary(LBtox3$percent.mt)