#DEG Analysis for LBTox Baseline no LB142 Github Copilot Code

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
library(ComplexHeatmap)
library(DESeq2)
library(pheatmap)
library(viridis)
library(EnhancedVolcano)
library(dittoSeq)
library(presto)
library(tictoc)
library(tibble)
library(RColorBrewer)

my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_SNP_Metadata_Scaled.RDS")
LBtox
colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#extract metadata table as CSV - just "Barcode" and "MajoritySinglet_Individual_Assignment"
metadata <- LBtox@meta.data

if ("Barcode" %in% colnames(metadata)) {
  selected_metadata <- metadata[, c("Barcode", "MajoritySinglet_Individual_Assignment")]
} else {
  # If "Barcode" is not a column, it might be the row names
  selected_metadata <- metadata[, "MajoritySinglet_Individual_Assignment", drop = FALSE]
  selected_metadata$Barcode <- rownames(metadata)
  selected_metadata <- selected_metadata[, c("Barcode", "MajoritySinglet_Individual_Assignment")]
}
print(selected_metadata)

#save metadata as CSV
write.csv(selected_metadata, file = "/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/Barcode_Individual_Metadata.csv", row.names = FALSE)

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Idents(Baseline.no142) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"



##please perform FindMarkers for all cell types between Cohort and Control (as above) and create heat maps for each cell type with more than 50 cells
# Subset cell types with more than 100 cells
cell_types <- names(table(Idents(Baseline.no142)))
cell_types <- cell_types[table(Idents(Baseline.no142)) > 100]

# Loop through each cell type
for (cell_type in cell_types) {
    # Subset the cell type
    subset_cells <- subset(Baseline.no142, subset = predicted.celltype.l2 == cell_type)
    
    # Perform FindMarkers between Cohort and Control
    markers <- FindMarkers(subset_cells, ident.1 = "Cohort", ident.2 = "Control")
    
    # Create volcano plot
    volcano_plot <- EnhancedVolcano(markers, lab = rownames(markers), x = "log2FoldChange", y = "p_val_adj")
    volcano_plot_file <- paste0("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/Plots/", cell_type, "_volcano_plot.png")
    ggsave(volcano_plot_file, volcano_plot)
    
    # Create heatmap
    heatmap <- DoHeatmap(subset_cells, features = rownames(markers)[1:50])
    heatmap_file <- paste0("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/Plots/", cell_type, "_heatmap.png")
    ggsave(heatmap_file, heatmap)
}