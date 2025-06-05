#rs11171739 and RPS26 eQTL analysis
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

LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_SNP_Metadata.RDS")
table <- LBtox@meta.data %>% as.data.table
table
LBtox
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

#Chi square analysis
tbl2 <- table(unique(LBtox@meta.data[,c("Population","rs11171739", "MajoritySinglet_Individual_Assignment")])[,c("Population","rs11171739")])
tbl2
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

chi_square_result<- chisq.test(tbl2)
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

#create RPS26
library(dplyr)
index = which(rownames(LBtox.baseline@assays$SCT@counts)=="RPS26")
meta = LBtox.baseline@meta.data
head(meta)
meta$RPS26 = LBtox.baseline@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Timepoint, Population, rs11171739) %>% reframe(count=n(),exp=(RPS26))
res = as.data.frame(res)
res


RPS26.metadata <- res %>%
group_by(predicted.celltype.l2, Population, rs11171739) %>%
summarise(features=n(), exp) %>%
group_by( Population, rs11171739)
head(RPS26.metadata)

rs11171739.palette <- c("CC"="#E69F00", "CT" = "#56B4E9", "TT" = "#999999" )

#chatGPT version
RPS26.plot.list <- lapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS26.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = rs11171739,
                y = exp,
                color = as.factor(rs11171739),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(rs11171739)),
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
        scale_colour_manual(values = rs11171739.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS26 expression by rs11171739") +
        ylim(c(0, 1 * max(RPS26.metadata$exp)))

        sigFunc <- function(x) {
  if (x < 0.001) {"***"}
  else if (x < 0.01) {"**"}
  else if (x < 0.05) {"*"}
  else {""}
}

         p1 <- p1 + geom_signif(
    comparisons = list(c("CC", "CT", "TT")),
    test = "wilcox.test",
    map_signif_level = sigFunc,
    colour = "black",
    y_position = 30
  )

  return(p1)
})


# Create a directory to save the individual plots
dir.create("individual_plots_eQTL_rs11171739", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS26.plot.list)) {
    filename <- paste0("individual_plots_eQTL_rs11171739/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS26.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#calculate p values
result.RPS26 <- sapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    x <- filter(RPS26.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[LBtox$rs11171739=="CC"]
    group2 <- x$exp[LBtox$rs11171739=="TT"]
    result.RPS26 <- wilcox.test(group1, group2)
    return(result.RPS26$p.value)
})

result.RPS26

