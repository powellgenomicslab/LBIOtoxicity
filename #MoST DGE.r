#MoST DGE

#Install packages
library(Seurat) 
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(jsonlite)
library(data.table)

unlink("/home/jenli3/R/x86_64-conda-linux-gnu-library/4.2/vctrs", recursive = TRUE)
install.packages("vctrs")

### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

#export graphics in VS Code
install.packages("httpgd")
devtools::install_github("nx10/httpgd")  
install.packages("languageserver")


#install Monocle3
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c("monocle"))
library(monocle)

#set working directory
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/mostpbmc"
setwd(my_dir)
getwd()

#import MoST PBMC RDS
mostpbmc <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/mostpbmc_azimuth.RDS")
mostpbmc
mostpbmc <- SetIdent(mostpbmc, value=pool)
Idents(mostpbmc) <- "orig.idents"


#visualising metadata frame 
head(mostpbmc@meta.data)
colnames(mostpbmc[[]])
table <- mostpbmc@meta.data %>% as.data.table
table
head(mostpbmc$ClinicalOutcome)
unique(mostpbmc$ClinicalOutcome)
unique(mostpbmc$CellType)
unique(mostpbmc$Timepoint)
unique(mostpbmc$Cancer)

#subset a population - eg tumour type
Colorectal <- subset(x = mostpbmc, subset = Cancer== c("Colorectal"))
HNSCC <- subset (x= mostpbmc, subset = Cancer ==c("H&N"))

#select a cell population
CD4 <- subset( x=mostpbmc, subset= CellType ==c("CD4+_cytotoxic_T_lymphocyte"))

#Set idents from a value in object metadata
Idents(CD4) <- 'ClinicalOutcome'
levels(CD4)

#Find Markers for CD4 cells between PR and PD
CD4.Response.Markers <- FindMarkers (CD4, ident.1="PD", ident.2="PR", test.usw=DE)

#view results
head(CD4.Response.Markers)

#set Idents to CELL TYPES
Idents(mostpbmc) <- 'predicted.celltype.l3'
levels(mostpbmc)

#find differentially expressed markers between CD4_TEM1 and CD8_TEM1
TEM.markers <-FindMarkers(mostpbmc, ident.1="CD4 TEM_1", ident.2 = "CD8 TEM_1")
head(TEM.markers)


#creating new variables for response vs non response
mostpbmc$Responder <- ifelse(mostpbmc$ClinicalOutcome %in% c("CR", "PR"), "Response", "Nonresponse")
mostpbmc_celltypes <- SplitObject(mostpbmc, "CellType")
mostpbmc_tumours <- SplitObject(mostpbmc, "Cancer")


##find markers response vs non response
Idents(mostpbmc) <- "Responder"
ResponseMarkers <- FindAllMarkers(mostpbmc, only.pos=TRUE)
head(ResponseMarkers)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5)
plot10 <- DoHeatmap(subset(mostpbmc, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot10, filename = "Heatmap_DGE_response2.png", height = 4.5, width = 7.5)

##Cytotoxic CD8 T cells response vs no response
Idents(mostpbmc_celltypes$`Cytotoxic_CD8+_Tcell_S100B+`) <- "Responder"
ResponseMarkers <- FindAllMarkers(mostpbmc, only.pos=TRUE)
head(ResponseMarkers)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5)
dim(ResponseMarkersFilter)
plot11 <- DoHeatmap(subset(mostpbmc_celltypes$`Cytotoxic_CD8+_Tcell_S100B+`, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot11, filename = "Heatmap_DGE_CytotoxicCD8.png", height = 4.5, width = 7.5)


#repeat code for all cell types
Idents(x) <- "Responder"
ResponseMarkers <- FindAllMarkers(mostpbmc, only.pos=TRUE)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5)
plot11 <- DoHeatmap(subset(x, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot11, filename = paste("Heatmap_DGE_cells_",celltype,".png"), height = 4.5, width = 7.5)

mapply(find_markers, mostpbmc_celltypes, names(mostpbmc_celltypes))

##Practising subsetting for DGE analysis
#Subset CD4 Cells
CD4 <- subset( x=mostpbmc, subset= CellType ==c("CD4+_cytotoxic_T_lymphocyte"))
unique(CD4$Timepoint)
CD4baseline <- subset (x=CD4, subset = Timepoint == c("Baseline"))
head(CD4baseline)
CD4wk4 <- subset (x=CD4, subset=Timepoint == c("Week4"))


#make heatmap
plot1 <- DoHeatmap(mostpbmc, features=NULL, group.by="ident", size=3)
ggsave(plot1, filename = "Heatmap_DGE_AzimuthCells.png", height = 4.5, width = 7.5)
Idents(CD4) <- subset(x=mostpbmc, subset= predicted.celltypel3 ==c("CD4 TEM_1"))
getwd()

#Subset timepoints
baseline <- subset (x=mostpbmc, subset=Timepoint == c("Baseline"))
head(baseline)
unique(baseline$CellType)
Wk8 <- subset (x=mostpbmc, subset=Timepoint ==c("Week8"))

#Subset Cell Types
CD4 <- subset(x=mostpbmc, subset=CellType ==c("CD4+_cytotoxic_T_lymphocyte"))
CD4
unique(CD4$Cancer)

#Subset Cancer Types
Colorectal <- subset(x = mostpbmc, subset = Cancer== c("Colorectal"))
head(Colorectal)
unique(Colorectal$CellType)
HNSCC <- subset (x= mostpbmc, subset = Cancer ==c("H&N"))
table(HNSCC$Cancer, HNSCC$CellType)
unique(HNSCC$CellType)

ColorectalCD4 <- subset(x=Colorectal, subset=CellType==c("CD4+_cytotoxic_T_lymphocyte"))
HNSCCCD4 <- subset (x=HNSCC, subset=CellType ==c("CD4+_cytotoxic_T_lymphocyte"))
HNSCCCD4

#stash identity classes


#Set Idents
Idents(object=mostpbmc) <- "ColorectalCD4"
CRC.markers <- FindMarkers(mostpbmc, ident.1= "ColorectalCD4")
CRC.markers

Idents(object=mostpbmc) <- "HNSCCCD4"
levels(x=mostpbmc)
HNSCC.markers <- FindMarkers(mostpbmc, ident.1="HNSCCCD4")
cells.1 <- WhichCells(object=mostpbmc, ident= "HNSCCCD4")
length(x= cells.1)
head(HNSCCCD4)



#Find markers


#performing DGE between two subsets - CD4 NK cells in HNSCC vs CD4NK cells in CRC


##subset out TWO cancer subtypes from object
CRCHN <- subset(x=mostpbmc, Cancer %in% c("Colorectal","HNSCC"))

#Identify what cells are present in that new object
table(CRCHN$Cancer, CRCHN$CellType)
Idents(CD4) <- 'Cancer'
head(CD4@meta.data)

#Find Markers of specific cell groups
CD4markers <- FindMarkers(CD4, ident.1="Colorectal", ident.2="HNSCC")

##
Idents(mostpbmc_celltypes$`Cytotoxic_CD8+_Tcell_S100B+`) <- "Responder"
ResponseMarkers <- FindAllMarkers(mostpbmc, only.pos=TRUE)
head(ResponseMarkers)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5)
dim(ResponseMarkersFilter)
plot11 <- DoHeatmap(subset(mostpbmc_celltypes$`Cytotoxic_CD8+_Tcell_S100B+`, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot11, filename = "Heatmap_DGE_CytotoxicCD8.png", height = 4.5, width = 7.5)
