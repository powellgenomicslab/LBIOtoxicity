#FOS and JUN analysis

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

#make metadata table for FOS
library(dplyr)
index = which(rownames(LBtox@assays$SCT@counts)=="FOS")
meta = LBtox@meta.data
head(meta)
meta$FOS = LBtox@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Timepoint, Population) %>% reframe(count=n(),exp=(FOS))
res = as.data.frame(res)
res

#make metadata table for FOS at baseline
index = which(rownames(baseline@assays$SCT@counts)=="FOS")
meta = baseline@meta.data
head(meta)
meta$FOS = baseline@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Population) %>% reframe(count=n(),exp=(FOS))
res = as.data.frame(res)
res

FOS.baseline.metadata <- res %>%
group_by(predicted.celltype.l2, Population) %>%
summarise(features=n(), exp) %>%
group_by(Population)
head(FOS.baseline.metadata)

#save res as a metadata frame with FOS
write.csv(res, file="Metadata_FOS.csv")

#Read it again to work on it later
res <- read.csv("Metadata_FOS.csv") 


fos.metadata <- res %>%
group_by(predicted.celltype.l2, Population, Timepoint) %>%
summarise(features=n(), exp) %>%
group_by(Population, Timepoint)
head(fos.metadata)

#boxplot with dots ggplot
Population.palette <- c("Cohort"="#E69F00", "Control" = "#56B4E9" )


#chatGPT version
FOS.plot.list <- lapply(unique(fos.metadata$predicted.celltype.l2), function(i){
    p1 <- fos.metadata %>%
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
        geom_quasirandom(dodge.width = 0.8,
                         width = 0.04,
                         alpha = 0.75) +
       geom_line(
            aes(x = Timepoint,
                y = median(exp),
                color = as.factor(Population)
            ),
            position = position_dodge(0.8),
            linewidth = 1.5
        ) +                    
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = Population.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("FOS expression") +
        ylim(c(0, 50))+
        scale_x_discrete(labels = c("Baseline", "C2", "C3")) 

sigFunc <- function(p) {
  if (p < 0.001) {"***"}
  else if (p < 0.01) {"**"}
  else if (p < 0.05) {"*"}
  else {""}
}

         p1 <- p1 + geom_signif(
    comparisons = list(c("Cohort", "Control")),
    test = "wilcox.test",
    map_signif_level = sigFunc,
    colour = "black",
    y_position = 30
  )

  return(p1)
})

# Create a directory to save the individual plots
dir.create("individual_plots_FOS", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(FOS.plot.list)) {
    filename <- paste0("individual_plots_FOS/plot_", i, ".pdf")
    ggsave(filename = filename, plot = FOS.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#calculate p values
result.FOS <- sapply(unique(fos.metadata$predicted.celltype.l2), function(i){
    x <- filter(fos.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[LBtox$Population=="Cohort"]
    group2 <- x$exp[LBtox$Population=="Control"]
    result.FOS <- wilcox.test(group1, group2)
    return(result.FOS$p.value)
})

result.FOS


#make metadata table for JUN
library(dplyr)
index = which(rownames(LBtox@assays$SCT@counts)=="JUN")
meta = LBtox@meta.data
head(meta)
meta$JUN = LBtox@assays$SCT@counts[index, ]
res2 = meta %>% group_by(predicted.celltype.l2, Timepoint, Population) %>% reframe(count=n(),exp=(JUN))
res2 = as.data.frame(res2)
res2

#save res as a metadata frame with FOS
write.csv(res2, file="Metadata_JUN.csv")

#Read it again to work on it later
res2 <- read.csv("Metadata_JUN.csv") 

#ggplot FOS at baseline
FOS.baseline.plot <- lapply(unique(FOS.baseline.metadata$predicted.celltype.l2), function(i){
    p1 <- FOS.baseline.metadata %>%
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
        ylab("FOS expression") +
        ylim(c(0, 50))

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
dir.create("individual_plots_FOSBaseline", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(FOS.baseline.plot)) {
    filename <- paste0("individual_plots_FOSBaseline/plot_", i, ".pdf")
    ggsave(filename = filename, plot = FOS.baseline.plot[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#calculate p value
result.baseline <- sapply(unique(FOS.baseline.metadata$predicted.celltype.l2), function(i){
    x <- filter(FOS.baseline.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[x$Population=="Cohort"]
    group2 <- x$exp[x$Population=="Control"]
    result.baseline <- wilcox.test(group1, group2)
    return(result.baseline$p.value)
})

result.baseline