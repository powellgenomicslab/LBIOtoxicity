# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(magrittr)
library(scProportionTest)
library(speckle)
library(ggbeeswarm)
library(ggsignif)
library(tidyverse)

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/jenli3/LBTOX")

# Read data
LBtox <- readRDS("SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")

# Remove LB142
LBtox.no142 <- subset(LBtox, subset = MajoritySinglet_Individual_Assignment != "LB142-C1")

# Create metadata summary
create_metadata_summary <- function(seurat_obj) {
  seurat_obj@meta.data %>%
    group_by(MajoritySinglet_Individual_Assignment, predicted.celltype.l2, Population, Timepoint) %>%
    summarise(counts = n(), .groups = "drop") %>%
    group_by(MajoritySinglet_Individual_Assignment, Population, Timepoint) %>%
    mutate(Percentage = counts / sum(counts) * 100)
}

Baseline.no142 <- subset(LBtox.no142, subset = Timepoint == "Baseline")
LBtox.metadatasummary <- create_metadata_summary(Baseline.no142)

# Permutation test
prop_test <- sc_utils(LBtox)
prop_test <- permutation_test(prop_test, cluster_identity = "predicted.celltype.l2", sample_1 = "Cohort", sample_2 = "Control", sample_identity = "Population")
ggsave(permutation_plot(prop_test))

# Define plot function
plot_cell_proportions <- function(metadata, group_col, colour_palette, comparisons, y_pos = 30) {
  lapply(unique(metadata$predicted.celltype.l2), function(cell_type) {
    p <- metadata %>%
      filter(predicted.celltype.l2 == cell_type) %>%
      ggplot(aes_string(x = group_col, y = "Percentage", color = group_col)) +
      geom_boxplot(width = 0.45, position = position_dodge(0.8), outlier.shape = NA, coef = 0) +
      geom_quasirandom(dodge.width = 0.8, width = 0.04, alpha = 0.75) +
      theme_classic() +
      theme(legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      xlab(cell_type) + ylab("Percentage of cells (%)") +
      scale_colour_manual(values = colour_palette) +
      geom_signif(comparisons = list(comparisons), test = "wilcox.test", map_signif_level = TRUE, colour = "black", y_position = y_pos)

    if (cell_type != "B intermediate") {
      p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
    }
    return(p)
  })
}

# Plot Cohort vs Control at Baseline
Population.palette <- c("Cohort" = "#F58800", "Control" = "#2c7e8c")
plots <- plot_cell_proportions(LBtox.metadatasummary, "Population", Population.palette, c("Cohort", "Control"))
ggsave("Cell_Proportion_Significance_Baseline.pdf", wrap_plots(plots, ncol = 4, guides = "collect"), width = 28, height = 45, units = "cm")

# Wilcoxon test for each cell type
pvals <- sapply(unique(LBtox.metadatasummary$predicted.celltype.l2), function(ct) {
  data <- filter(LBtox.metadatasummary, predicted.celltype.l2 == ct)
  wilcox.test(data$Percentage[data$Population == "Cohort"], data$Percentage[data$Population == "Control"])$p.value
})

# Output p-values
print(pvals)

# Stacked bar chart for baseline cohort vs control
metadata.no142 <- Baseline.no142@meta.data %>%
  group_by(Population, predicted.celltype.l2) %>%
  summarise(counts = n(), .groups = "drop") %>%
  group_by(Population) %>%
  mutate(percentage = counts / sum(counts) * 100)

p_bar <- ggplot(metadata.no142, aes(fill = predicted.celltype.l2, y = percentage, x = Population)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("Population") +
  ylab("Percentage") +
  ggtitle("Cell type percentages: Cohort vs Control") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_discrete("Cell type")

ggsave("Baseline_Cohort_vs_Control_Barchart_no142.pdf", p_bar, width = 20, height = 20, units = "cm")