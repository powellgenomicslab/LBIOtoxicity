# Ribosomal Gene QC Analysis: ruling out technical artifact

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(patchwork)

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/jenli3/LBTOX")

# Load data
LBtox <- readRDS("SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
DefaultAssay(LBtox) <- "RNA"
Idents(LBtox) <- "predicted.celltype.l2"

# Exclude problematic sample LB142
LBtox.no142 <- subset(LBtox, subset = MajoritySinglet_Individual_Assignment != "LB142-C1")

# Assign operator metadata manually (JL, KP, PK)
# Create a reference table and apply via a join or match
operator_map <- data.frame(
  Sample = c("LBIO181219A-2", "LBIO050220A-5", "LB019", "LB035", "LB020", "LB014", "LB042", "LB062",
             "LB128", "LB029", "LB091", "LB098", "LB121", "LB105", "LB142-C1", "LB207-C1"),
  Pool = c("LBTox_pool1", "LBTox_pool2", "LBTox_pool1", "LBTox_pool3", "LBTox_pool2", 
           "LBTox_pool1", "LBTox_pool2", "LBTox_pool3",
           "LBTox_pool2", "LBTox_pool3", "LBTox_pool1", "LBTox_pool1", "LBTox_pool3",
           "LBTox_pool2", "LBTox_pool1", "LBTox_pool3"),
  Operator = c(rep("JL", 8), rep("KP", 6), "KP", "PK")
)

for (i in 1:nrow(operator_map)) {
  match_rows <- which(LBtox@meta.data$MajoritySinglet_Individual_Assignment == operator_map$Sample[i] & 
                      LBtox@meta.data$Pool == operator_map$Pool[i])
  LBtox@meta.data$Operator[match_rows] <- operator_map$Operator[i]
}

# Calculate ribosomal content
LBtox.no142[["percent.rb"]] <- PercentageFeatureSet(LBtox.no142, pattern = "^RP[LS]")

# Summarize
summary(LBtox.no142$percent.rb)

# Violin plot: Ribosomal and mitochondrial gene content by Operator
Idents(LBtox.no142) <- "Operator"
vln_plot <- VlnPlot(LBtox.no142, features = c("percent.rb", "percent.mt"), pt.size = 0)
ggsave(vln_plot, filename = "QC_Ribosomal_Mitochondrial_ByOperator.png", height = 4.5, width = 7.5)

# Histogram of ribosomal percentage
qc.metrics <- as_tibble(LBtox.no142[[]], rownames = "Cell.Barcode")
hist_plot <- ggplot(qc.metrics, aes(x = percent.rb)) + 
  geom_histogram(binwidth = 0.5, fill = "gold", color = "black") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  ggtitle("Distribution of Ribosomal Gene Percentages") +
  xlab("Percent Ribosomal Genes") + ylab("Cell Count")
ggsave(hist_plot, filename = "Histogram_Ribosomal_Percentages.png", height = 4.5, width = 7.5)

# FeaturePlots of ribosomal genes split by population
ribosomal_genes <- c("RPS26", "RPS24", "RPS27", "RPS9", "RPS4Y1", "RPS4X")
fp_all <- FeaturePlot(LBtox.no142, features = ribosomal_genes, split.by = "Population")
ggsave(fp_all, filename = "FeaturePlot_RibosomalGenes_SplitByPopulation.png", height = 15, width = 7.5)

# Individual FeaturePlots
for (gene in c("RPS26", "RPS24", "RPS4Y1")) {
  p <- FeaturePlot(LBtox.no142, features = gene, split.by = "Population")
  ggsave(p, filename = paste0("FeaturePlot_", gene, "_SplitByPopulation.png"), height = 4.5, width = 7.5)
}
