#JUN analysis
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
library(DESeq2)
library(pheatmap)
library(viridis)
library(EnhancedVolcano)
library(dittoSeq)
library(presto)
library(ggsci)
library(ggbeeswarm)
library(ggsignif)

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_CompleteMetadata.RDS")
baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))

#make metadata table for JUN
index = which(rownames(LBtox@assays$SCT@counts)=="JUN")
meta = LBtox@meta.data
head(meta)
meta$JUN = LBtox@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Timepoint, Population) %>% reframe(count=n(),exp=(JUN))
res = as.data.frame(res)
res

#make metadata table for JUN across timepoints
index = which(rownames(LBtox@assays$SCT@counts)=="JUN")
meta = LBtox@meta.data
head(meta)
meta$JUN = LBtox@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Population, Timepoint) %>% reframe(count=n(),exp=(JUN))
res = as.data.frame(res)
res

#Read it again to work on it later
res2 <- read.csv("Metadata_JUN.csv") 

#metadata table for JUN at baseline
index = which(rownames(baseline@assays$SCT@counts)=="JUN")
meta = baseline@meta.data
head(meta)
meta$JUN = baseline@assays$SCT@counts[index, ]
res1 = meta %>% group_by(predicted.celltype.l2, Population) %>% reframe(count=n(),exp=(JUN))
res1 = as.data.frame(res)
res1

JUN.baseline.metadata <- res1 %>%
group_by(predicted.celltype.l2, Population) %>%
summarise(features=n(), exp) %>%
group_by(Population)
head(JUN.baseline.metadata)

#ggplot JUN at baseline
JUN.baseline.plot <- lapply(unique(JUN.baseline.metadata$predicted.celltype.l2), function(i){
    p1 <- JUN.baseline.metadata %>%
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
        ylab("JUN expression") +
        ylim(c(0, 1 * max(JUN.baseline.metadata$exp)))

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
dir.create("individual_plots_JUNbaseline", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(JUN.baseline.plot)) {
    filename <- paste0("individual_plots_JUNbaseline/plot_", i, ".pdf")
    ggsave(filename = filename, plot = JUN.baseline.plot[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}



#calculate p value
result.baseline <- sapply(unique(JUN.baseline.metadata$predicted.celltype.l2), function(i){
    x <- filter(JUN.baseline.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[x$Population=="Cohort"]
    group2 <- x$exp[x$Population=="Control"]
    result.baseline <- wilcox.test(group1, group2)
    return(result.baseline$p.value)
})

result.baseline