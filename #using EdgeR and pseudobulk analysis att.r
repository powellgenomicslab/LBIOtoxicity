#using EdgeR and pseudobulk analysis attempt

install.packages("scater")
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)
library(limma)

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_Merged_updatedmetadata_postQC.RDS")
LBtox
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#Extract counts
counts <- GetAssayData(object=LBtox,slot="counts", assay="RNA")
counts

#create single-cell experiment object
sce <- SingleCellExperiment(assays=list(counts=counts),
                            colData=metadata)

                            