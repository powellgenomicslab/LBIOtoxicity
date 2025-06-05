#Complex Heatmap DGE B cells

#DGE B intermediate and CD8 TEM cells

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
library(DropletUtils)


my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
LBtox
colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#Read CSV
DEGS_Bint <- read.csv("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/B intermediate cell markers at BASELINE control vs cohort nil LB142.csv")
DEGS_Bint
DEGS_Bint <- DEGS_Bint %>%
            tibble::rownames_to_column(var="Gene")

#select the top25 genes for plotting
TopGenes <- rownames(DEGS_Bint[order(DEGS_Bint$p_val_adj, decreasing=FALSE)[1:25],])
DEGS_Bint$TopGenes <- ifelse(DEGS_Bint$Gene %in% TopGenes, TRUE, FALSE)

#select top 15
DEGS_Bint %>% group_by("Population") %>% top_n(n=15, wt=avg_log2FC) -> top15
toplot_gene <- top15$X

#scale the mat
mat <- LBtox[["RNA"]]@data[toplot_gene,] %>% as.matrix()
mat <- t(scale(t(mat)))
head(mat)


Population <- LBtox@meta.data$Population
head(LBtox@meta.data)

library(circlize)
col

#making the heatmap
plot1 <- Heatmap(mat)
ggsave(plot1, filename="B intermediate cells heatmap.pdf", height=7.5, width=4.5)

