#Finding Markers for Each Cluster

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

#import seurat object
pbmc.tox.azi <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/pbmctoxazi_scale_normalize.RDS")

#visualising metadata frame 
head(pbmc.tox.azi@meta.data)
colnames(LBtox[[]])
unique(LBtox$MajoritySinglet_Individual_Assignment)
unique(LBtox$Pool)
table <- LBtox.no142@meta.data %>% as.data.table
table
colnames(LBtox.no142[[]])

#add column - Population (Cohort vs Control)
#Cohort = LB002 (LBIO181219A-2), LB005 (LBIO050220A-5), LB035, LB029, LB091, LB098, LB105, LB142
#Control = LB019, LB062, LB128, LB020, LB013, LB121, LB042, LB207
LBtox$Population <- ifelse(LBtox$MajoritySinglet_Individual_Assignment %in% c("LBIO181219A-2","LBIO050220A-5","LB035","LB029", "LB091",  "LB098" , "LB105","LB142", "LB142-C1", "LB128"), "Cohort", "Control")
unique(LBtox$Population)
index(LBtox$Population)

#Add column - treatment (PD1 vs IO_Combo)
#PD1 =  LB002, LB019, LB062, LB098, LB121, LB105, LB042
#IO_Combo = LB035, LB128, LB029, LB020, LB091, LB014, LB142, LB207
LBtox$Treatment <- ifelse (LBtox$MajoritySinglet_Individual_Assignment %in% c("LBIO181219A-2", "LBIO050220A-5", "LB019","LB062", "LB098", "LB121" ,"LB105", "LB042"), "PD1", "IO_Combo")  

#add data column
#timepoint - baseline, C2, C3
#baseline = ALL THE REST
#C2 = "LB035" in "LBTox_pool2", LB091 in "LBTox_pool3", "LB105" in "LBTox_pool2"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB035" & LBtox$Pool =="LBTox_pool2"] <- "C2"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB091" & LBtox$Pool =="LBTox_pool3"] <- "C2"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB105" & LBtox$Pool =="LBTox_pool2"] <- "C2"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB128" & LBtox$Pool =="LBTox_pool1"] <- "C2"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB014" & LBtox$Pool =="LBTox_pool3"] <- "C2"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB042" & LBtox$Pool =="LBTox_pool3"] <- "C2"

#C3 = "LB091" in "LBTox_pool2", "LB105" in "LBTox_pool1"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB091" & LBtox$Pool =="LBTox_pool2"] <- "C3"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB105" & LBtox$Pool =="LBTox_pool1"] <- "C3"
LBtox@meta.data$Timepoint[LBtox$MajoritySinglet_Individual_Assignment == "LB014" & LBtox$Pool =="LBTox_pool1"] <- "C3"


#baseline = the rest
LBtox$Timepoint[is.na(LBtox$Timepoint)] <- "Baseline"
unique(LBtox$Timepoint)
unique(LBtox$predicted.celltype.l2)

LBtox@meta.data$Timepoint

#add data column #Cancer - RCC, melanoma, NSCLC, cSCC
#RCC = "LB091","LB029", "LB020", "LB014", "LB142-C1"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "RCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "RCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "RCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "RCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "RCC"

#melanoma = "LB035", "LB128", "LB207-C1"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "Melanoma"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "Melanoma"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "Melanoma"

#NSCLC = "LBIO050220A-5", "LBIO181219A-2", "LB019", "LB062"  
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LBIO050220A-5"] <- "NSCLC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LBIO181219A-2"] <- "NSCLC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LB019"] <- "NSCLC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LB062"] <- "NSCLC"

#SCC = "LB098" , "LB121", "LB105" , "LB042" 
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LB098"] <- "SCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LB121"] <- "SCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LB105"] <- "SCC"
LBtox@meta.data$Cancer[LBtox@meta.data$MajoritySinglet_Individual_Assignment =="LB042"] <- "SCC"

unique(LBtox$Cancer)

#Add meta data for sex
LBtox@meta.data$Sex[LBtox$MajoritySinglet_Individual_Assignment == "LB035" & LBtox$Pool =="LBTox_pool2"] <- "Female"
LBtox@meta.data$Sex[LBtox$MajoritySinglet_Individual_Assignment == "LB035" & LBtox$Pool =="LBTox_pool3"] <- "Female"
unique(LBtox$Sex)
LBtox$Sex[is.na(LBtox$Sex)] <- "Male"

#Add metadata for rs1131017 genotype
#CC for 
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB019"] <- "CC"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "CC"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "CC"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB098"] <- "CC"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "CC"

#CG rs1131017 genotype
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO050220A-5"] <- "CG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO181219A-2"] <- "CG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "CG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "CG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB105"] <- "CG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "CG"

#GG rs1131017 genotype
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB062"] <- "GG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "GG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "GG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB121"] <- "GG"
LBtox@meta.data$rs1131017[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB042"] <- "GG"

##Change CC to 0
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB019"] <- "0"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "0"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "0"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB098"] <- "0"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "0"

#CG rs1131017 genotype
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO050220A-5"] <- "1"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO181219A-2"] <- "1"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "1"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "1"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB105"] <- "1"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "1"

#GG rs1131017 genotype
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB062"] <- "2"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "2"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "2"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB121"] <- "2"
LBtox@meta.data$SNP[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB042"] <- "2"

#set Idents for Cell type
unique(pbmc.tox.azi$predicted.celltype.l2)
Idents(pbmc.tox.azi) <- 'predicted.celltype.l2'
levels(pbmc.tox.azi)

#set Idents for Population
Idents(pbmc.tox.azi) <- 'Population'

LBtox

#Save Seurat Object with All Meta Data
saveRDS(LBtox, "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
saveRDS(LBtox, "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_CompleteMetadata_V2.RDS")
saveRDS(LBtox, "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_postQC_postLogNormalize.RDS")
