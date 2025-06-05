
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
library(speckle)
devtools::install_github("rpolicastro/scProportionTest")

#set wd
#set working directory
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#Read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
LBtox

#Summary of cell clusters and numbers - baseline
colnames(LBtox[[]])

#Cohort_Population_Subset_analysis
#subset out a single indvidual and cluster
Baseline <- subset(x=LBtox, subset= Timepoint == c("Baseline"))
Baseline
Cohort.Baseline <- subset(x=Baseline, subset=Population==c("Cohort"))
Cohort.Baseline
Control.Baseline <- subset(x=Baseline, subset=Population==c("Control"))
Control.Baseline

unique(LBtox$predicted.celltype.l2)
unique(LBtox$MajoritySinglet_Individual_Assignment)

saveRDS(Baseline@meta.data, '/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/Baseline_Metadata.RDS')

#Cell numbers
Idents(Baseline) <- "predicted.celltype.l2"
Idents(Cohort.Baseline) <- "predicted.celltype.l2"
Idents(Control.Baseline) <- "predicted.celltype.l2"
table(Idents(Cohort.Baseline))
table(Idents(Control.Baseline))
table(Idents(Baseline.no142))

#Cell proportions
prop.table(table(Idents(Cohort.Baseline)))*100
prop.table(table(Idents(Control.Baseline)))*100

res <- t.test(Control.Baseline$predicted.celltype.l2=="CD8 TEM", Cohort.Baseline$predicted.celltype.l2=="CD8 TEM")
res

res2 <- t.test(Control.Baseline$predicted.celltype.l2=="B intermediate", Cohort.Baseline$predicted.celltype.l2=="B intermediate")
res2

#Subset Out Cohort Patients
Cohort <- subset(x=LBtox, subset= Population == c("Cohort"))
Cohort
Cohort.Baseline <- subset(x=Cohort, subset=Timepoint==c("Baseline"))
Cohort.C2 <- subset(x=Cohort, subset=Timepoint==c("C2"))
Cohort.C3 <- subset(x=Cohort, subset=Timepoint==c("C3"))

#baseline vs C2 cell counts and proportions
LBtox
LBtoxC2 <- subset(x=LBtox, subset=Timepoint==c("C2"))
LBtoxC2.cohort <- subset(x=LBtoxC2, subset=Population==c("Cohort"))
LBtoxC2.cohort
LBtoxC2.control <- subset(x=LBtoxC2, subset=Population==c("Control"))
LBtoxC2.control
Idents(LBtoxC2.cohort) <- "predicted.celltype.l2"
Idents(LBtoxC2.control) <- "predicted.celltype.l2"
table(Idents(LBtoxC2.cohort))
prop.table(table(Idents(LBtoxC2.cohort)))*100
table(Idents(LBtoxC2.control))
prop.table(table(Idents(LBtoxC2.control)))*100

#repeating above analysis without LB142
LBtox.no142

#baseline - control vs cohort NO LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Baseline.no142 <- subset(x=LBtox.no142, subset= Timepoint == c("Baseline"))
Baseline.no142
Cohort.Baseline.no142 <- subset(x=Baseline.no142, subset=Population==c("Cohort"))
Cohort.Baseline.no142
Control.Baseline.no142 <- subset(x=Baseline.no142, subset=Population==c("Control"))
Control.Baseline.no142

#baseline control vs cohort NO LB142 tables
Idents(Cohort.Baseline.no142) <- "predicted.celltype.l2"
Idents(Control.Baseline.no142) <- "predicted.celltype.l2"
table(Idents(Cohort.Baseline.no142))
prop.table(table(Idents(Cohort.Baseline.no142)))*100
table(Idents(Control.Baseline.no142))
prop.table(table(Idents(Control.Baseline.no142)))*100

#C2 control vs cohort NO LB142
Btox
LBtoxC2.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("C2"))
LBtoxC2.cohort.no142 <- subset(x=LBtoxC2.no142, subset=Population==c("Cohort"))
LBtoxC2.cohort
LBtoxC2.control.no142 <- subset(x=LBtoxC2.no142, subset=Population==c("Control"))
LBtoxC2.control.no142
Idents(LBtoxC2.cohort.no142) <- "predicted.celltype.l2"
Idents(LBtoxC2.control.no142) <- "predicted.celltype.l2"
table(Idents(LBtoxC2.cohort.no142))
prop.table(table(Idents(LBtoxC2.cohort.no142)))*100
table(Idents(LBtoxC2.control.no142))
prop.table(table(Idents(LBtoxC2.control.no142)))*100

#C2 IN GENRAL
C2.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("C2"))
Idents(C2.no142) <- "predicted.celltype.l2"
table(Idents(C2.no142))

#UMAP

p1 <- DimPlot(LBtoxC2.cohort, reduction="umap", group.by = "predicted.celltype.l2", label.size = 3, repel=TRUE, label=TRUE)
p2 <- DimPlot(LBtoxC2.control, reduction="umap", group.by = "predicted.celltype.l2", label.size = 3, repel=TRUE, label=TRUE)
plot1 <- p1 + p2
ggsave(plot1,filename="UMAP_Cohort vs Control_C2.png", width = 15.5, height = 10.5 )

#visualising cohort patients only
Cohort <- subset(LBtox, subset=Population==c("Cohort"))

#Stacked Barchart
Cohort.metadata <- Cohort@meta.data %>%
    group_by(Timepoint, predicted.celltype.l2) %>%
    summarise(counts = n()) %>%
    group_by(Timepoint) %>%
    mutate(percentage = counts/sum(counts)*100)
​

# Stacked bar chart
library(tidyverse)
install.packages("tidyverse")

p1 <- Cohort.metadata %>%
    ggplot(aes(fill=predicted.celltype.l2, y=percentage, x=Timepoint)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab("Timepoint") +
    ylab("Percentage") +
    theme_bw() +
    ggtitle("Cell type percentages: cohort group over timepoints") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_discrete("Cell type")
​
ggsave(filename = "cohort_celltype_percentage_stackedbarchart.pdf", plot = p1, device = "pdf", width = 20, height = 20, units = "cm")
​


#Stacked Barchart
metadata.summary <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/Baseline_Metadata.RDS") %>%
    group_by(Timepoint, predicted.celltype.l2) %>%
    summarise(counts = n()) %>%
    group_by(Timepoint) %>%
    mutate(percentage = counts/sum(counts)*100)
​
​
​
# Plot stacked barchart ---------------------------------------------------
​
# Stacked bar chart
library(tidyverse)
install.packages("tidyverse")

Cohort.metadata <- Cohort@meta.data %>%
    group_by(Timepoint, predicted.celltype.l2) %>%
    summarise(counts = n()) %>%
    group_by(Timepoint) %>%
    mutate(percentage = counts/sum(counts)*100)

metadata.no142 <- Baseline.no142@meta.data %>%
    group_by(Population, predicted.celltype.l2) %>%
    summarise(counts=n()) %>%
    group_by(Population) %>%
    mutate(percentage=counts/sum(counts)*100)


p1 <- metadata.no142 %>%
    ggplot(aes(fill=predicted.celltype.l2, y=percentage, x=Population)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab("Population") +
    ylab("Percentage") +
    theme_bw() +
    ggtitle("Cell type percentages: cohort vs control") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_discrete("Cell type")
​
ggsave(filename = "Baseline Cohort vs Control Barchart no 142.pdf",
       plot = p1,
       device = "pdf", width = 20, height = 20, units = "cm")
​

--

# Define cell types that represent at least %5 in one sample
common.cell.types <- c("B intermediate","B naive","CD14 Mono","CD16 Mono",
                       "CD4 TCM", "CD8 TEM", "NK")
​
# Read in and transform data
metadata.summary <- readRDS("~/Downloads/Baseline_Metadata.RDS") %>%
    mutate(predicted.celltype.l2 = case_when(predicted.celltype.l2 %in% common.cell.types ~ predicted.celltype.l2,
                                             TRUE ~ "other")) %>%
    group_by(Population, predicted.celltype.l2) %>%
    summarise(counts = n()) %>%
    group_by(Population) %>%
    mutate(percentage = counts/sum(counts)*100)
​
​
​
# Plot stacked barchart ---------------------------------------------------
​
# Stacked bar chart
p1 <- metadata.summary %>%
    ggplot(aes(fill=predicted.celltype.l2, y=percentage, x=Population)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab("Population") +
    ylab("Percentage") +
    theme_bw() +
    ggtitle("Cell type percentages: cohort vs control") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_discrete("Cell type")
​
ggsave(filename = "~/Downloads/stacked-barchart.pdf",
       plot = p1,
       device = "pdf", width = 20, height = 20, units = "cm")