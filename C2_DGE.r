#Cycle 2 cohort vs control DEG analysis
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
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)

#C2 no LB142
C2.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("C2"))
C2.cohort.no142 <- subset(x=C2.no142, subset=Population==c("Cohort"))
C2.control.no142 <- subset(x=C2.no142, subset=Population==c("Control"))

#Set Idents
Idents(C2.cohort.no142) <- "predicted.celltype.l2"
Idents(C2.control.no142) <- "predicted.celltype.l2"
Idents(C2.no142) <- "predicted.celltype.l2"
Idents(Baseline.no142) <- "predicted.celltype.l2"
Idents(Cohort.no142) <- "predicted.celltype.l2"

#identify numbers of cells in each category
table(Idents(C2.cohort.no142))
table(Idents(C2.control.no142))
table(Idents(C2.no142), C2.no142$Population)
table(Idents(Cohort.no142), Cohort.no142$Timepoint)

#Cell proportions
prop.table(table(Idents(Cohort.Baseline)))*100

#CD14 Monocytes
C2.CD14.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(C2.CD14.no142) <- "Population"

#DEGS Analysis CD14 Mono
DefaultAssay(C2.CD14.no142) <- "RNA"
DEGS_CD14Mono_C2 <- FindMarkers(C2.CD14.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
head(DEGS_CD14Mono_C2)
DEGS_CD14Mono_C2 <- DEGS_CD14Mono_C2 %>%
            tibble::rownames_to_column(var="Gene")

write.csv(DEGS_CD14Mono_C2, "CD14 Mono Differentially expressed Genes C2.CSV")   

#Volcano Plot CD14 Mono
plot1 <- EnhancedVolcano(DEGS_CD14Mono_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD14Mono_C2$Gene, FCcutoff=0.5)
  ggsave(plot1,filename="Volcano Plot CD14 Monocytes at C2.png", height = 7.5, width = 10)

  #CD4 TCM 
C2.CD4.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD4 TCM"))
Idents(C2.CD4.no142) <- "Population"

#DEGS Analysis CD4 TCM
DefaultAssay(C2.CD4.no142) <- "RNA"
DEGS_CD4_TCM_C2 <- FindMarkers(C2.CD4.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4_TCM_C2 <- DEGS_CD4_TCM_C2 %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4_TCM_C2)

write.csv(DEGS_CD4_TCM_C2, "CD4 TCM Differentially expressed Genes at C2 cohort vs control.CSV")   

#Volcano Plot CD4 TCM
plot2 <- EnhancedVolcano(DEGS_CD4_TCM_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4_TCM_C2$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot CD4 TCM C2.png", height = 7.5, width = 7.5)

#B naive (949 cells)
C2.Bnaive.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("B naive"))
Idents(C2.Bnaive.no142) <- "Population"

  #DEGS Analysis B naive

DefaultAssay(C2.Bnaive.no142) <- "RNA"
DEGS_Bnaive_C2 <- FindMarkers(C2.Bnaive.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_Bnaive_C2 <- DEGS_Bnaive_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_Bnaive_C2)

write.csv(DEGS_Bnaive_C2, "B naive cells Differentially expressed Genes at C2.CSV")   

#Volcano Plot B naive
plot3 <- EnhancedVolcano(DEGS_Bnaive_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bnaive_C2$Gene, FCcutoff=0.1)
  ggsave(plot3,filename="Volcano Plot B naive DEGs at C2.png", height = 7.5, width = 7.5)

#NK 
C2.NK.no142 <- subset(x= C2.no142, subset=predicted.celltype.l2==c("NK"))
Idents(C2.NK.no142) <- "Population"

#DEGS Analysis NK cells
DefaultAssay(C2.NK.no142) <- "RNA"
DEGS_NK_C2 <- FindMarkers(C2.NK.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_NK_C2 <- DEGS_NK_C2 %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_NK_C2)

write.csv(DEGS_NK_C2, "NK cells Differentially expressed Genes.CSV")   

#Volcano Plot NK cells
plot4 <- EnhancedVolcano(DEGS_NK_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_NK_C2$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot NK cells DEGs at C2.png", height = 7.5, width = 7.5)


#CD16 Mono (847 cells)
C2.CD16.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD16 Mono"))
Idents(C2.CD16.no142) <- "Population"

#DEGS Analysis CD16 mono cells
DefaultAssay(C2.CD16.no142) <- "RNA"
DEGS_CD16_C2 <- FindMarkers(C2.CD16.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD16_C2 <- DEGS_CD16_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD16_C2)

write.csv(DEGS_CD16_C2, "CD16 Monocyte cells Differentially expressed Genes at C2 cohort vs control.CSV")   

#Volcano Plot CD16 cells
plot1 <- EnhancedVolcano(DEGS_CD16_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD16_C2$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot CD16 mono cells DEGs at C2 cohort vs control.png", height = 7.5, width = 7.5)

#CD4 naive (415 cells)
C2.CD4naive.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD4 Naive"))
Idents(C2.CD4naive.no142) <- "Population"

#DEGS Analysis CD4 naive cells
DefaultAssay(C2.CD4naive.no142) <- "RNA"
DEGS_CD4n_C2 <- FindMarkers(C2.CD4naive.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4n_C2 <- DEGS_CD4n_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD4n_C2)

write.csv(DEGS_CD4n_C2, "CD4 Naive cells Differentially expressed Genes.CSV")   

#Volcano Plot CD4 naive cells
plot6 <- EnhancedVolcano(DEGS_CD4n_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4n_C2$Gene, FCcutoff=0.1)
  ggsave(plot6,filename="Volcano Plot CD4 naive cells DEGs C2.png", height = 7.5, width = 7.5)

#CD8 naive (415 cells)
C2.CD8naive.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD8 Naive"))
Idents(C2.CD8naive.no142) <- "Population"

#DEGS Analysis CD8 naive cells
DefaultAssay(C2.CD8naive.no142) <- "RNA"
DEGS_CD8n_C2 <- FindMarkers(C2.CD8naive.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8n_C2 <- DEGS_CD8n_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD8n_C2)

write.csv(DEGS_CD8n_C2, "CD8 Naive cells Differentially expressed Genes C2.CSV")   

#Volcano Plot CD8 naive cells
plot7 <- EnhancedVolcano(DEGS_CD8n_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8n_C2$Gene, FCcutoff=0.1)
  ggsave(plot7,filename="Volcano Plot CD8 naive cells DEGs C2.png", height = 7.5, width = 7.5)

#CD8 TEM
C2.CD8TEM.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD8 TEM"))
Idents(C2.CD8TEM.no142) <- "Population"

#DEGS Analysis CD8 TEM
DefaultAssay(C2.CD8TEM.no142) <- "RNA"
DEGS_CD8T_C2 <- FindMarkers(C2.CD8TEM.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8T_C2 <- DEGS_CD8T_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD8T_C2)

write.csv(DEGS_CD8T_C2, "CD8 TEM cells Differentially expressed Genes C2.CSV")   

#Volcano Plot CD8 TEM cells
plot8 <- EnhancedVolcano(DEGS_CD8T_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8T_C2$Gene, FCcutoff=0.1)
  ggsave(plot8,filename="Volcano Plot CD8 TEM cells DEGs C2.png", height = 7.5, width = 7.5)

  #CD4 TEM
C2.CD4TEM.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD4 TEM"))
Idents(C2.CD4TEM.no142) <- "Population"

#DEGS Analysis CD4 TEM
DefaultAssay(C2.CD4TEM.no142) <- "RNA"
DEGS_CD4T_C2 <- FindMarkers(C2.CD4TEM.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4T_C2 <- DEGS_CD4T_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD4T_C2)

write.csv(DEGS_CD4T_C2, "CD4 TEM cells Differentially expressed Genes C2.CSV")   

#Volcano Plot CD8 TEM cells
plot9 <- EnhancedVolcano(DEGS_CD4T_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4T_C2$Gene, FCcutoff=0.1)
  ggsave(plot9,filename="Volcano Plot CD4 TEM cells DEGs C2.png", height = 7.5, width = 7.5)

#gdT
C2.gdT.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("gdT"))
Idents(C2.gdT.no142) <- "Population"

#DEGS Analysis gdT 
DefaultAssay(C2.gdT.no142) <- "RNA"
DEGS_gdT_C2 <- FindMarkers(C2.gdT.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_gdT_C2 <- DEGS_gdT_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_gdT_C2)

write.csv(DEGS_gdT_C2 , "gdT cells Differentially expressed Genes C2.CSV")   

#Volcano Plot gdT cells
plot9 <- EnhancedVolcano(DEGS_gdT_C2 , x="avg_log2FC", y="p_val_adj", lab=DEGS_gdT_C2$Gene, FCcutoff=0.1)
  ggsave(plot9,filename="Volcano Plot gdT cells DEGs C2.png", height = 7.5, width = 7.5)

#B intermediate
C2.Bint.no142 <- subset(x= C2.no142, subset=predicted.celltype.l2==c("B intermediate"))
Idents(C2.Bint.no142) <- "Population"

#DEGS Analysis B intermediate cells
DefaultAssay(C2.Bint.no142) <- "RNA"
DEGS_Bint_C2 <- FindMarkers(C2.Bint.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_Bint_C2 <- DEGS_Bint_C2 %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_Bint_C2)

write.csv(DEGS_Bint_C2, "B Intermediate cells Differentially expressed Genes C2 cohort vs control.CSV")   

#Volcano Plot B intermediate cells
plot1 <- EnhancedVolcano(DEGS_Bint_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bint_C2$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot B intermediate cells DEGs at C2.png", height = 7.5, width = 7.5)

  #CD14 Mono 
C2.CD14.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(C2.CD14.no142) <- "Population"

#DEGS Analysis CD14 mono cells
DefaultAssay(C2.CD14.no142) <- "RNA"
DEGS_CD14_C2 <- FindMarkers(C2.CD14.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD14_C2 <- DEGS_CD14_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD14_C2)

write.csv(DEGS_CD14_C2, "CD14 Monocyte cells Differentially expressed Genes at C2.CSV")   

#Volcano Plot CD16 cells
plot5 <- EnhancedVolcano(DEGS_CD14_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD14_C2$Gene, FCcutoff=0.1)
  ggsave(plot5,filename="Volcano Plot CD14 mono cells DEGs at C2.png", height = 7.5, width = 7.5)

#CD8 TCM
C2.CD8TCM.no142 <- subset(x=C2.no142, subset=predicted.celltype.l2==c("CD8 TCM"))
Idents(C2.CD8TCM.no142) <- "Population"

#DEGS Analysis CD8 TCM mono cells
DefaultAssay(C2.CD8TCM.no142) <- "RNA"
DEGS_CD8TCM_C2 <- FindMarkers(C2.CD8TCM.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8TCM_C2 <- DEGS_CD8TCM_C2 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD8TCM_C2)

write.csv(DEGS_CD8TCM_C2, "CD8 TCM cells Differentially expressed Genes at C2 cohort vs control.CSV")   

#Volcano Plot CD8 TCM cells
plot4 <- EnhancedVolcano(DEGS_CD8TCM_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8TCM_C2$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot CD8 TCM cells DEGs at C2.png", height = 7.5, width = 7.5)

#T reg cells
C2.Treg.no142 <- subset(x= C2.no142, subset=predicted.celltype.l2==c("Treg"))
Idents(C2.Treg.no142) <- "Population"

#DEGS Analysis T reg cells
DefaultAssay(C2.Treg.no142) <- "RNA"
DEGS_Treg_C2 <- FindMarkers(C2.Treg.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_Treg_C2 <- DEGS_Treg_C2 %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_Treg_C2)

write.csv(DEGS_Treg_C2, "T reg cells Differentially expressed Genes C2 cohort vs control.CSV")   

#Volcano Plot T reg cells
plot2 <- EnhancedVolcano(DEGS_Treg_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_Treg_C2$Gene, FCcutoff=0.1)
  ggsave(plot2,filename="Volcano Plot T reg cells DEGs at C2.png", height = 7.5, width = 7.5)

#gdT cells
C2.gdT.no142 <- subset(x= C2.no142, subset=predicted.celltype.l2==c("gdT"))
Idents(C2.gdT.no142) <- "Population"

#DEGS Analysis gdT cells
DefaultAssay(C2.gdT.no142) <- "RNA"
DEGS_gdT_C2 <- FindMarkers(C2.gdT.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_gdT_C2 <- DEGS_gdT_C2 %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_gdT_C2)

write.csv(DEGS_gdT_C2, "gdT cells Differentially expressed Genes C2 cohort vs control.CSV")   

#Volcano Plot gdT cells
plot3 <- EnhancedVolcano(DEGS_gdT_C2, x="avg_log2FC", y="p_val_adj", lab=DEGS_gdT_C2$Gene, FCcutoff=0.1)
  ggsave(plot3,filename="Volcano Plot gdT cells DEGs at C2.png", height = 7.5, width = 7.5)