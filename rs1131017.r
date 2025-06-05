
#Rs1131017 analysis

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
table <- LBtox@meta.data %>% as.data.table
table
LBtox
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

table2 <- table(unique(LBtox@meta.data[,c("Population","rs1131017", "MajoritySinglet_Individual_Assignment")])[,c("Population","rs1131017")])
table2
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

#barplot for rs1131017 distribution
plot4 <- ggplot(LBtox@meta.data) + aes(x=Population, fill=rs1131017) + geom_bar(position="dodge")
ggsave(plot4, filename="rs1131017 distribution according to toxicity group separate bars 4.png", height=7.5, width=7.5)

tbl = table(unique(LBtox@meta.data[,c("Population", "MajoritySinglet_Individual_Assignment")]))
tbl

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

#RPS26 eQTL analysis ggplot
#try with ggplot2
#make metadata table
c

#boxplot with dots ggplot
Population.palette <- c("Cohort"="#E69F00", "Control" = "#56B4E9" )


#chatGPT version
RPS26.plot.list <- lapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS26.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = rs1131017,
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
        facet_wrap(~rs1131017, scale = "free") +
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
        ylab("RPS26 expression by rs1131017 genotype") +
        ylim(c(0, 1.1 * max(RPS26.metadata$exp)))

sigFunc <- function(x) {
  if (x < 0.001) {"***"}
  else if (x < 0.01) {"**"}
  else if (x < 0.05) {"*"}
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

# Function for significance annotation


ggsave(filename="Plotting RPS26 expression according to cell type and rs1131017 genotype.pdf", plot= wrap_plots(RPS.plot.list, ncol=2, guides="collect"),
device="pdf", height=40, width=20, units="cm")

# Create a directory to save the individual plots
dir.create("individual_plots_eQTL_2", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS.plot.list)) {
    filename <- paste0("individual_plots_eQTL_2/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS26.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
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

##RPS26 eQTL analysis NOT separated by population
rs1131017.palette <- c("CC"="#E69F00", "CG" = "#56B4E9", "GG" = "#999999" )

#chatGPT version
RPS26.plot.list <- lapply(unique(RPS26.metadata$predicted.celltype.l2), function(i){
    p1 <- RPS26.metadata %>%
        filter(predicted.celltype.l2 == i) %>%
        ggplot(
            mapping = aes(
                x = rs1131017,
                y = exp,
                color = as.factor(rs1131017),
            )
        ) +
        geom_boxplot(
            aes(color = as.factor(rs1131017)),
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
        scale_colour_manual(values = rs1131017.palette, name = "") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab(i) +
        theme(axis.title.x = element_text(size = 7)) +
        ylab("RPS26 expression by rs1131017 genotype") +
        ylim(c(0, 1.1 * max(RPS26.metadata$exp)))

})


# Create a directory to save the individual plots
dir.create("individual_plots_eQTL_2", showWarnings = FALSE)

# Save each plot separately
for (i in seq_along(RPS26.plot.list)) {
    filename <- paste0("individual_plots_eQTL_2/plot_", i, ".pdf")
    ggsave(filename = filename, plot = RPS26.plot.list[[i]], device = "pdf", height = 8, width = 10, units = "cm")
}

#regression analysis of RPS26 genotype and populations per cell type
Idents(LBtox) <- "predicted.celltype.l2"
levels(LBtox$predicted.celltype.l2)
class(LBtox$predicted.celltype.l2)
Idents(LBtox)

LBtox$predicted.celltype.l2 <- as.factor(LBtox$predicted.celltype.l2)

# Verify the levels
if(is.factor(LBtox$predicted.celltype.l2)) {LBtox$predicted.celltype.l2 <- as.character(LBtox$predicted.celltype.l2)}
LBtox$predicted.celltype.l2 <- as.factor(LBtox$predicted.celltype.l2)

# Set the cell identity
Idents(LBtox) <- "predicted.celltype.l2"
levels(LBtox$predicted.celltype.l2)
class(LBtox$predicted.celltype.l2)

# Subset for "CD14 Mono"
CD14_Mono_subset <- subset(LBtox, idents = "CD14 Mono")

 Set cell identities
Idents(LBtox) <- "predicted.celltype.l2"
str(LBtox)

remotes::install_version("Matrix", version = "1.6.1")

# Define the cell types you want to subset
selected_cell_types <- c("CD14 Mono", "CD4 Proliferating", "B intermediate", "CD4 TCM", "CD4 TEM", "CD8 TCM", "CD8 TEM", "NK", "Treg", "B memory", "CD16 Mono")

CD14.mono <- subset(LBtox,  subset=predicted.celltype.l2=="CD14 Mono")
CD14_Mono_subset <- subset(LBtox, subset = predicted.celltype.l2 == "CD14 Mono")
CD14_Mono_subset <- LBtox[LBtox$predicted.celltype.l2 == "CD14 Mono", ]

CD14.mono

dim(LBtox)

#Angli code regression
RPS26.metadata <- res %>%
group_by(predicted.celltype.l2, Population, SNP) %>%
summarise(features=n(), exp) %>%
group_by( Population, SNP)
head(RPS26.metadata)


summary(lm(exp~SNP * Population, LBtox))

install.packages("Matrix", version = "1.6-1", repos = "http://cran.us.r-project.org")