#DGE CD8 TEM cells

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
library(tictoc)

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
LBtox
colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Idents(Baseline) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#baseline no LB142 and LB035 (female
Baseline.New <- subset(Baseline.no142, subset=Sex==c("Male"))

#baseline CD8 TEM cells subset
Baseline.CD8.no142 <- subset(x=Baseline.New, subset=predicted.celltype.l2==c("CD8 TEM"))
Baseline.CD8 <- ScaleData(object = Baseline.Bcells, verbose=FALSE)
Idents(Baseline.CD8.no142) <- "Population"

#DEGS Analysis
DefaultAssay(Baseline.CD8.no142) <- "RNA"
DEGS_CD8_Males <- FindMarkers(Baseline.CD8.no142, ident.1="Cohort", ident.2="Control", vars.to.regress="Pool")
head(DEGS_CD8_Males)
DEGS_CD8_Males <- DEGS_CD8_Males %>%
            tibble::rownames_to_column(var="Gene")

write.csv(DEGS_CD8_Males, "CD8 TEM Differentially expressed Genes Males Only.CSV")   

#Male only plot
TopGenes <- rownames(DEGS_Bint_Males[order(DEGS_Bint_Males$p_val_adj, decreasing=FALSE)[1:25],])
DEGS_Bint_Males$TopGene <- ifelse(DEGS_Bint_Males$Gene %in% TopGenes, TRUE, FALSE)

#Enhanced Volcano Plot - youtube tutorial
plot2 <- EnhancedVolcano(DEGS_CD8_Males, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8_Males$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot CD8 TEM Cells top 10 genes Males.png", height = 7.5, width = 7.5)


#Idents for baseline cells
Idents(Baseline.New) <- "predicted.celltype.l2"
table(Idents(Baseline.New))

unique(Baseline.New$predicted.celltype.l2)
