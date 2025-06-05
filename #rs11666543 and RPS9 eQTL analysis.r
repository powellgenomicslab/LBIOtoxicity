#rs11666543 and RPS9 eQTL analysis
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

LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_SNP_Metadata.RDS"
table <- LBtox@meta.data %>% as.data.table
table
LBtox
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment")

#Chi square analysis
tbl <- table(unique(LBtox@meta.data[,c("Population","rs11666543", "MajoritySinglet_Individual_Assignment")])[,c("Population","rs11666543")])
tbl
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

chi_square_result<- chisq.test(tbl)
chi_square_result

alpha = 0.05
if (chi_square_result$p.value < alpha) {
  if (chi_square_result$statistic > 0) {
    cat("Reject the null hypothesis. There is a significant association in the expected direction.\n")
  } else {
    cat("Reject the null hypothesis. There is a significant association in the opposite direction.\n")
  }
} else {
  cat("Fail to reject the null hypothesis. There is no significant association.\n")
}

LBtox.baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))

#create RPS9
index = which(rownames(LBtox.baseline@assays$SCT@counts)=="RPS9")
meta = LBtox.baseline@meta.data
head(meta)
meta$RPS9 = LBtox.baseline@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Population, rs11666543) %>% reframe(count=n(),exp=(RPS9))
res = as.data.frame(res)
res

RPS9.metadata <- res %>%
group_by(predicted.celltype.l2, Population, rs11666543) %>%
summarise(features=n(), exp) %>%
group_by( Population, rs11666543)
head(RPS9.metadata)

rs11666543.palette <- c("AA"="#E69F00", "AG" = "#56B4E9", "GG" = "#999999" )

#chatGPT version
RPS9.plot.list <- lapply(unique(RPS9.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS9.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = rs11666543,
                y = exp,
                color = as.factor(rs11666543),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(rs11666543)),
            width = 0.45,
            position = position_dodge(0.8),
            outlier.shape = NA,
            coef = 0
        ) +
        geom_quasirandom(dodge.width = 0.8,
                         width = 0.04,
                         alpha = 0.75) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = rs11666543.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS9 expression by rs11666543 genotype") +
        ylim(c(0, 1 * max(RPS9.metadata$exp)))

})

# Create a directory to save the individual plots
dir.create("individual_plots_eQTL_RPS9", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS9.plot.list)) {
    filename <- paste0("individual_plots_eQTL_RPS9/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS9.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}


