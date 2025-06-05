#ribosomal gene analysis

#load packages
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
install.packages("SoupX")
library(SoupX)
library(DropletUtils)
library(DoubletFinder)
library(knitr)

#set working directory
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

#import raw data to look into ambient ribosomal RNA

#identify ribosomal genes 
grep("^RP[LS]", rownames(LBtox.no142@assays$RNA@counts), value=TRUE)
LBtox.no142[["percent.rb"]] <- PercentageFeatureSet(LBtox, pattern="^RP[LS]")
head(LBtox.no142$percent.rb)
Idents(LBtox) <- "Operator"
summary(LBtox.no142$percent.rb)

#violin plot according to pools
plot4 <- VlnPlot(LBtox, features = c("percent.rb", "percent.mt"))
ggsave(plot4, filename="Ribosomal and Mitochondrial gene distribution by Operator whole LBtox.png", height = 4.5, width = 7.5)

#distribution of ribosomal genes percentages
as_tibble(LBtox.no142[[]], rownames="Cell.Barcode") -> qc.metrics
head(qc.metrics)

plot2 <- qc.metrics %>%
    ggplot(aes(percent.rb)) + 
    geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
    ggtitle("Distribution of Ribosomal Percentages") +
    geom_vline(xintercept=10)


ggsave(plot2, filename="Distribution ribosomal percentages.png", height= 4.5, width=7.5)

#Feature Plot
plot1 <- FeaturePlot(LBtox.no142, features=c("RPS26","RPS24", "RPS27","RPS9", "RPS4Y1", "RPS4X"), split.by="Population")
ggsave(plot1, filename="Feature Plot of Ribosomal Genes Cohort vs Control.png", height=15, width=7.5)

#smaller plots with individual genes
plot2 <- FeaturePlot(LBtox.no142, features=c("RPS26"), split.by="Population")
ggsave(plot2, filename="Feature Plot of RPS26 Cohort vs Control.png", height=4.5, width=7.5)

plot3 <- FeaturePlot(LBtox.no142, features=c("RPS24"), split.by="Population")
ggsave(plot3, filename="Feature Plot of RPS24 Cohort vs Control.png", height=4.5, width=7.5)

plot4 <- FeaturePlot(LBtox.no142, features=c("RPS4Y1"), split.by="Population")
ggsave(plot4, filename="Feature Plot of RPS4Y1 Cohort vs Control.png", height=4.5, width=7.5)

#subset according to PBMC processing method (JL, KP, TK)
#Operator - JL, KP, TK
#JL = LB002_C
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LBIO181219A-2" & LBtox$Pool =="LBTox_pool1"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LBIO050220A-5" & LBtox$Pool =="LBTox_pool2"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB019" & LBtox$Pool =="LBTox_pool1"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB035" & LBtox$Pool =="LBTox_pool3"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB035" & LBtox$Pool =="LBTox_pool2"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB020" & LBtox$Pool =="LBTox_pool2"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB014" & LBtox$Pool =="LBTox_pool1"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB014" & LBtox$Pool =="LBTox_pool2"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB014" & LBtox$Pool =="LBTox_pool3"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB042" & LBtox$Pool =="LBTox_pool2"] <- "JL"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB062" & LBtox$Pool =="LBTox_pool3"] <- "JL"

#KP operator
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB128" & LBtox$Pool =="LBTox_pool2"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB128" & LBtox$Pool =="LBTox_pool1"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB128" & LBtox$Pool =="LBTox_pool2"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB029" & LBtox$Pool =="LBTox_pool3"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB091" & LBtox$Pool =="LBTox_pool1"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB091" & LBtox$Pool =="LBTox_pool2"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB091" & LBtox$Pool =="LBTox_pool3"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB098" & LBtox$Pool =="LBTox_pool1"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB121" & LBtox$Pool =="LBTox_pool3"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB105" & LBtox$Pool =="LBTox_pool2"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB105" & LBtox$Pool =="LBTox_pool1"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB042" & LBtox$Pool =="LBTox_pool3"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB142-C1" & LBtox$Pool =="LBTox_pool1"] <- "KP"
LBtox@meta.data$Operator[LBtox$MajoritySinglet_Individual_Assignment == "LB207-C1" & LBtox$Pool =="LBTox_pool3"] <- "PK"


unique(LBtox$Operator)
head(LBtox$Operator)