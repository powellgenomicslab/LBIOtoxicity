#rs6032664 analysis
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
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_SNP_Metadata.RDS")
table <- LBtox@meta.data %>% as.data.table
table
LBtox
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

#Add metadata for rs6032664 genotype
#TT for 
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO050220A-5"] <- "TT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LBIO181219A-2"] <- "TT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB062"] <- "TT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB035"] <- "TT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB029"] <- "TT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB098"] <- "TT"

#AT rs6032664 genotype
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB019"] <- "AT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB020"] <- "AT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB091"] <- "AT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB014"] <- "AT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB121"] <- "AT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB105"] <- "AT"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB042"] <- "AT"

#AA rs6032664 genotype
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB207-C1"] <- "AA"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB142-C1"] <- "AA"
LBtox@meta.data$rs6032664[LBtox@meta.data$MajoritySinglet_Individual_Assignment == "LB128"] <- "AA"

unique(LBtox$rs6032664)

#save new metadata with new SNP 
saveRDS(LBtox, "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_CompleteMetadata_V3.RDS")

#Chi square analysis
table2 <- table(unique(LBtox@meta.data[,c("Population","rs6032664", "MajoritySinglet_Individual_Assignment")])[,c("Population","rs6032664")])
table2
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

chi_square_result<- chisq.test(table2)
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


#barplot for rs6032664 distribution
plot4 <- ggplot(LBtox@meta.data) + aes(x=Population, fill=rs6032664) + geom_bar(position="dodge")
ggsave(plot4, filename="rs1131017 distribution according to toxicity group separate bars 4.png", height=7.5, width=7.5)

LBtox.baseline <- subset(x=LBtox, subset=Timepoint==c("Baseline"))

#Genes of interest
index = which(rownames(LBtox.baseline@assays$SCT@counts)=="RPS9")
meta = LBtox.baseline@meta.data
head(meta)
meta$RPS9 = LBtox.baseline@assays$SCT@counts[index, ]
res = meta %>% group_by(predicted.celltype.l2, Timepoint, Population, rs6032664) %>% reframe(count=n(),exp=(RPS9))
res = as.data.frame(res)
res

RPS9.metadata <- res %>%
group_by(predicted.celltype.l2, Population, rs6032664) %>%
summarise(features=n(), exp) %>%
group_by( Population, rs6032664)
head(RPS9.metadata)

rs6032664.palette <- c("AA"="#E69F00", "AT" = "#56B4E9", "TT" = "#999999" )

#chatGPT version
RPS9.plot.list <- lapply(unique(RPS9.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS9.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = rs6032664,
                y = exp,
                color = as.factor(rs6032664),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(rs6032664)),
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
        scale_colour_manual(values = rs6032664.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS9 expression by rs6032664") +
        ylim(c(0, 1 * max(RPS9.metadata$exp)))

        sigFunc <- function(x) {
  if (x < 0.001) {"***"}
  else if (x < 0.01) {"**"}
  else if (x < 0.05) {"*"}
  else {""}
}

         p1 <- p1 + geom_signif(
    comparisons = list(c("TT", "AT", "AA")),
    test = "wilcox.test",
    map_signif_level = sigFunc,
    colour = "black",
    y_position = 30
  )

  return(p1)
})


# Create a directory to save the individual plots
dir.create("individual_plots_eQTL_rs6032664", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS9.plot.list)) {
    filename <- paste0("individual_plots_eQTL_rs6032664/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS9.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#calculate p values
result.RPS9 <- sapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    x <- filter(RPS26.metadata, predicted.celltype.l2==i)
    group1 <- x$exp[LBtox$rs11171739=="CC"]
    group2 <- x$exp[LBtox$rs11171739=="TT"]
    result.RPS26 <- wilcox.test(group1, group2)
    return(result.RPS26$p.value)
})

result.RPS26


