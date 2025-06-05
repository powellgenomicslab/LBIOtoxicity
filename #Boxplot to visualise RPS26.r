#Boxplot to visualise RPS26 

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(jsonlite)
library(data.table)
library(tidyverse)
library(devtools)
library(Azimuth)
library(ComplexHeatmap)
library(pheatmap)
library(viridis)
library(EnhancedVolcano)
library(dittoSeq)
library(presto)
library(ggbeeswarm)
library(ggpubr)
library(ggsignif)


my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_CompleteMetadata.RDS")
colnames(LBtox[[]])
Idents(LBtox.no142) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"
table <- LBtox.no142@meta.data %>% as.data.table
table
head(LBtox.no142[[]])
rownames(LBtox.no142)
table <- LBtox.no142@meta.data %>% as.data.table
table
colnames(LBtox.no142)
dim(LBtox.no142)
genes <- rownames(LBtox.no142)
genes

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Idents(LBtox.no142) <- "predicted.celltype.l2"
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Cohort.no142 <- subset(x=LBtox.no142, subset=Population==c("Cohort"))
Control.no142 <- subset(x=LBtox.no142, subset=Population==c("Control"))
Idents(Baseline.no142) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#RPS26 in NK cells
plot1 <- dittoBoxPlot(Cohort.NK, "RPS26", group.by="Timepoint")
ggsave(plot1, filename="NK cohort cells RPS26 boxplot.png", height=6.5, width=6.5)

plot2 <- dittoBoxPlot(Control.NK, "RPS26", group.by="Timepoint")
ggsave(plot2, filename="NK control cells RPS26 boxplot.png", height=6.5, width=6.5)

#significance test with NK cells
cohort.NK <- GetAssayData(subset(Co))

#RPS26 in CD14 Monocytes
Cohort.CD14Mono <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Control.CD14Mono <- subset(x=Control.no142, subset=predicted.celltype.l2==c("CD14 Mono"))

plot3 <- dittoBoxPlot(Cohort.CD14Mono, "RPS26", group.by="Timepoint")
ggsave(plot3, filename="CD14 Monocytes cohort cells RPS26 boxplot.png", height=6.5, width=6.5)

plot4 <- dittoBoxPlot(Control.CD14Mono, "RPS26", group.by="Timepoint")
ggsave(plot4, filename="CD14 Monocytes Control cells RPS26 boxplot.png", height=6.5, width=6.5)

#CD4 TCM
Cohort.CD4TCM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD4 TCM"))
Control.CD4TCM <- subset(x=Control.no142, subset=predicted.celltype.l2==c("CD4 TCM"))

plot5 <- dittoBoxPlot(Cohort.CD4TCM, "RPS26", group.by="Timepoint")
ggsave(plot5, filename="CD4 TCM cohort cells RPS26 boxplot.png", height=6.5, width=6.5)

plot6 <- dittoBoxPlot(Control.CD4TCM, "RPS26", group.by="Timepoint")
ggsave(plot6, filename="CD4 TCM Control cells RPS26 boxplot.png", height=6.5, width=6.5)

#CD8 TEM cells
Cohort.CD8TEM <- subset(x=Cohort.no142, subset=predicted.celltype.l2==c("CD8 TEM"))
Control.CD8TEM <- subset(x=Control.no142, subset=predicted.celltype.l2==c("CD8 TEM"))

plot7 <- dittoBoxPlot(Cohort.CD8TEM, "RPS26", group.by="Timepoint")
ggsave(plot7, filename="CD8 TEM cohort cells RPS26 boxplot.png", height=6.5, width=6.5)

plot8 <- dittoBoxPlot(Control.CD8TEM, "RPS26", group.by="Timepoint")
ggsave(plot8, filename="CD8 TEM Control cells RPS26 boxplot.png", height=6.5, width=6.5)


#Create new metadata with RPS26
library(dplyr)
index = which(rownames(LBtox@assays$SCT@counts)=="RPS26")
meta = LBtox@meta.data
head(meta)
meta$RPS26 = LBtox@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Timepoint, Population) %>% reframe(count=n(),exp=(RPS26))
res = as.data.frame(res)
res

#try with ggplot2
#make metadata table
LBtox.metadatasummary <- LBtox.no142@meta.data %>% 
group_by(MajoritySinglet_Individual_Assignment, predicted.celltype.l2, Population, Timepoint) %>%
summarise(counts=n()) %>%
group_by(MajoritySinglet_Individual_Assignment, Population, Timepoint) %>%
mutate(Percentage=counts/sum(counts)*100)
head(LBtox.metadatasummary)

RPS26.metadata <- res %>%
group_by(predicted.celltype.l2, Population, Timepoint) %>%
summarise(features=n(), exp) %>%
group_by( Population, Timepoint)
head(RPS26.metadata)

#create CSV with RPS26 data
write.csv(RPS26.metadata, file="Metadata_RPS26.csv")

#boxplot with dots ggplot
Population.palette <- c("Cohort"="#E69F00", "Control" = "#56B4E9" )


#chatGPT version
RPS26.metadata <- read.csv("Metadata_RPS26.csv") 
RPS.plot.list <- lapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS26.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = Timepoint,
                y = exp,
                color = as.factor(Population),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(Population)),
            width = 0.45,
            position = position_dodge(0.8),
            outlier.shape = NA,
            coef = 0
        ) +
stat_summary(
            fun = median,
            geom = "line",
            aes(group = Timepoint, color = as.factor(Population)),
            position = position_dodge(0.8),
            size = 0.8
        ) +  
        geom_quasirandom(dodge.width = 0.8,
                         width = 0.04,
                         alpha = 0.75) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = Population.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS26 expression") +
        ylim(c(0, 1.1 * max(RPS26.metadata$exp)))

    return(p1)
})

ggsave(filename="Plotting RPS26 expression according to cell type 4.pdf", plot= wrap_plots(RPS.plot.list, ncol=2, guides="collect"),
device="pdf", height=40, width=20, units="cm")

# Create a directory to save the individual plots
dir.create("individual_plots_updatedRPS26", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS.plot.list)) {
    filename <- paste0("individual_plots_updatedRPS26/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#calculate p values
result.RPS26 <- sapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    x <- filter(RPS26.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[LBtox$Population=="Cohort"]
    group2 <- x$exp[LBtox$Population=="Control"]
    result.RPS26 <- wilcox.test(group1, group2)
    return(result.RPS26$p.value)
})

result.RPS26


#violin plot expression
features <- c("RPS26", "RPS4X", "RPS4Y1")
plot1 <- VlnPlot(LBtox.no142, features = features, slot=counts)
ggsave(plot1, filename="RPS Gene Expression per cell type 4.png", height=10, width=25)


#Just baseline TOX vs CONTROL RPS26 expression
#chatGPT version
RPS.plot.list <- lapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS26.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = Timepoint$Baseline,
                y = exp,
                color = as.factor(Population),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(Population)),
            width = 0.45,
            position = position_dodge(0.8),
            outlier.shape = NA,
            coef = 0
        )  +
        geom_quasirandom(dodge.width = 0.8,
                         width = 0.04,
                         alpha = 0.75) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = Population.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS26 expression") +
        ylim(c(0, 1.1 * max(RPS26.metadata$exp)))

        # Apply wilcoxon signed-rank test between tox/control pairs
        sigFunc <- function(x) {
            if (x < 0.001) {"***"}
            else if (x < 0.01) {"**"}
            else if (x < 0.05) {"*"}
            else {NA}
        }

        p1 <- p1 + geom_signif(
            comparisons = list(c("Cohort", "Control")),
            test = "wilcox.test",
            map_signif_level = sigFunc,
            colour = "black",
            y_position = 30
        )
})


# Create a directory to save the individual plots
dir.create("individual_plots_RPS26Baseline", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS.plot.list)) {
    filename <- paste0("individual_plots_RPS26Baseline/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#Plot out RPS26 expression at BASELINE only between toxicity and control
baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))

#make RPS26 a variable
library(dplyr)
index = which(rownames(baseline@assays$SCT@counts)=="RPS26")
meta = baseline@meta.data
head(meta)
meta$RPS26 = baseline@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Population) %>% reframe(count=n(),exp=(RPS26))
res = as.data.frame(res)
res

#make metadata table
baseline.metadatasummary <- baseline@meta.data %>% 
group_by(MajoritySinglet_Individual_Assignment, predicted.celltype.l2, Population) %>%
summarise(counts=n()) %>%
group_by(MajoritySinglet_Individual_Assignment, Population) %>%
mutate(Percentage=counts/sum(counts)*100)
head(baseline.metadatasummary)

RPS26.baseline.metadata <- res %>%
group_by(predicted.celltype.l2, Population) %>%
summarise(features=n(), exp) %>%
group_by( Population)
head(RPS26.baseline.metadata)

#ggplot RPS26 at baseline
RPS.baseline.plot <- lapply(unique(RPS26.baseline.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS26.baseline.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = Population,
                y = exp,
                color = as.factor(Population),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(Population)),
            width = 0.45,
            position = position_dodge(0.8),
            outlier.shape = NA,
            coef = 0
        )  +
        geom_quasirandom(dodge.width = 0.8,
                         width = 0.04,
                         alpha = 0.75) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = Population.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS26 expression") +
        ylim(c(0, 1.1 * max(RPS26.baseline.metadata$exp)))

        # Apply wilcoxon signed-rank test between tox/control pairs
        sigFunc <- function(x) {
            if (x < 0.001) {"***"}
            else if (x < 0.01) {"**"}
            else if (x < 0.05) {"*"}
            else {NA}
        }

        p1 <- p1 + geom_signif(
            comparisons = list(c("Cohort", "Control")),
            test = "wilcox.test",
            map_signif_level = sigFunc,
            colour = "black",
            y_position = 30
        )
})

# Create a directory to save the individual plots
dir.create("individual_plots_RPS26Baseline", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS.baseline.plot)) {
    filename <- paste0("individual_plots_RPS26Baseline/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS.baseline.plot[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#calculate p value
result.baseline <- sapply(unique(RPS26.baseline.metadata$predicted.celltype.l2), function(i){
    x <- filter(RPS26.baseline.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[x$Population=="Cohort"]
    group2 <- x$exp[x$Population=="Control"]
    result.baseline <- wilcox.test(group1, group2)
    return(result.baseline$p.value)
})

result.baseline

