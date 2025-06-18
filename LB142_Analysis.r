#exploring LB142-C1
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

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBtox_PostQC_Annotated.RDS")
unique(LBtox$MajoritySinglet_Individual_Assignment)

#subset out LB142-C1
LB142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"))
LB142

#cluster LB142
plot2 <- DimPlot(LB142, reduction="umap", group.by = "predicted.celltype.l2", label.size = 3, repel=TRUE, label=TRUE)
ggsave(plot2, filename="LB142_UMAP_annotated.png", width = 8.5, height = 6.5)

#cell numbers and proportions
#set Idents
Idents(LB142) <- 'predicted.celltype.l2'
levels(LB142)

#summary of cells in LB142
table(Idents(LB142))
prop.table(table(Idents(LB142)))

#Find Markers for B intermediate cells
B_int.markers <- FindMarkers(LB142, ident.1="B intermediate")
head(B_int.markers)
write.csv(B_int.markers, "B intermediate cell markers.csv")

#caspase markers for LB142
plot5 <- FeaturePlot(LB142, features=c("CASP3"))
ggsave(plot5, filename="Caspase Expression on LB142 UMAP.png",width = 8.5, height = 6.5)

#EXCLUDE LB142-C1
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)

#cluster in cell type
plot3 <- DimPlot(LBtox.no142, reduction="umap", group.by = "predicted.celltype.l2", label.size = 3, repel=TRUE, label=TRUE)
ggsave(plot3, filename="UMAP_noLB142_annotated.png", width = 8.5, height = 6.5)

#cluster in cohort vs control
plot4 <- DimPlot(LBtox.no142, reduction="umap", group.by = "Population", label.size = 3, repel=TRUE, label=TRUE)
ggsave(plot4, filename="UMAP_noLB142_cohortvscontrol.png", width = 8.5, height = 6.5)

#cell numbers and proportions
Idents(LBtox.no142) <- 'predicted.celltype.l2'
table(Idents(LBtox.no142))
prop.table(table(Idents(LBtox.no142)))

#control vs cohort without LB142
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Cohort.Baseline.no142 <- subset(x=Baseline.no142, subset=Population==c("Cohort"))
Cohort.Baseline.no142
Control.Baseline.no142
Control.Baseline.no142 <- subset(x=Baseline.no142, subset=Population==c("Control"))

#numbers and proportions
Idents(Baseline.no142) <- "predicted.celltype.l2"
Idents(Cohort.Baseline.no142) <- "predicted.celltype.l2"
Idents(Control.Baseline.no142) <- "predicted.celltype.l2"
table(Idents(Cohort.Baseline.no142))
table(Idents(Control.Baseline.no142))

#Cell proportions
prop.table(table(Idents(Cohort.Baseline.no142)))
prop.table(table(Idents(Control.Baseline.no142)))


#barplot
LBtox.no142.metadata <- LBtox.no142@meta.data %>%
    group_by(Timepoint, predicted.celltype.l2) %>%
    summarise(counts = n()) %>%
    group_by(Timepoint) %>%
    mutate(percentage = counts/sum(counts)*100)

p1 <- LBtox.no142@meta.data %>%
    ggplot(aes(fill="predicted.celltype.l2", y=percentage, x=Population)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab("Population") +
    ylab("Percentage") +
    theme_bw() +
    ggtitle("Cell type percentages: cohort vs control no LB142") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_discrete("Cell type")
​
ggsave(filename = "Cell proportions NO LB142.pdf",
       plot = p1,
       device = "pdf", width = 20, height = 20, units = "cm")