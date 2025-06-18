#Add SNP Metadata

#import packages
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

#set Wd
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)
getwd()

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_CompleteMetadata_V3.RDS")

#Add metadata for rs11171739 genotype
#CC for 
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB019"] <- "CC"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "CC"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "CC"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB098"] <- "CC"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "CC"

#CT for:
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO181219A-2"] <- "CT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO050220A-5"] <- "CT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "CT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "CT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "CT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB105"] <- "CT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "CT"

#TT for:
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB062"] <- "TT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB042"] <- "TT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "TT"
LBtox@meta.data$rs11171739[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB121"] <- "TT"

unique(LBtox$rs11171739)


#Add metadata for rs4899554 genotype
#CC genotype
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB019"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO181219A-2"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO050220A-5"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB062"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB098"] <- "CC"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB105"] <- "CC"

#CT genotype
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "CT"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB042"] <- "CT"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB121"] <- "CT"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "CT"
LBtox@meta.data$rs4899554[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "CT"

#Add rs11666543 genotype
#GG genotype
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO181219A-2"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB062"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB042"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB098"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB105"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB121"] <- "GG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "GG"

#AG genotype
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB019"] <- "AG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "AG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "AG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "AG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "AG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "AG"
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "AG"

#AA genotype
LBtox@meta.data$rs11666543[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO050220A-5"] <- "AA"


saveRDS(LBtox, "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_SNP_Metadata.RDS")
