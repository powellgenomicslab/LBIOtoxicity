library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(EnhancedVolcano)
library(patchwork)
library(RColorBrewer)

# Setup and data loading
setwd("/directflow/SCCGGroupShare/projects/jenli3/LBTOX")
LBtox <- readRDS("SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
DefaultAssay(LBtox) <- "RNA"
Idents(LBtox) <- "predicted.celltype.l2"

# Remove LB142 and subset baseline
LBtox.no142 <- subset(LBtox, subset = MajoritySinglet_Individual_Assignment != "LB142-C1")
Baseline.no142 <- subset(LBtox.no142, subset = Timepoint == "Baseline")
Idents(Baseline.no142) <- "predicted.celltype.l2"

# Create output folders
dir.create("Volcano_Plots", showWarnings = FALSE)
dir.create("Heatmap_Plots", showWarnings = FALSE)
dir.create("DGE_Results", showWarnings = FALSE)

# Define minimum cells threshold
min_cells <- 100

# Get cell types meeting threshold
celltypes <- Baseline.no142@meta.data %>% 
  group_by(predicted.celltype.l2) %>% 
  summarise(n = n()) %>% 
  filter(n > min_cells) %>% 
  pull(predicted.celltype.l2)

# Loop over cell types
for (ct in celltypes) {
  cat("Analyzing cell type:", ct, "\n")
  
  # Subset Seurat object
  subset_obj <- subset(Baseline.no142, subset = predicted.celltype.l2 == ct)
  Idents(subset_obj) <- "Population"
  DefaultAssay(subset_obj) <- "RNA"
  
  # Run DEG analysis
  degs <- FindMarkers(subset_obj, ident.1 = "Cohort", ident.2 = "Control", vars.to.regress = c("Pool", "Sex")) %>%
    rownames_to_column("Gene")
  
  # Save raw results
  write.csv(degs, paste0("DGE_Results/", gsub(" ", "_", ct), "_DEGs.csv"), row.names = FALSE)
  
  # Volcano plot
  volcano <- EnhancedVolcano(degs,
                             x = "avg_log2FC",
                             y = "p_val_adj",
                             lab = degs$Gene,
                             FCcutoff = 0.5,
                             pCutoff = 0.05,
                             title = paste("Volcano Plot:", ct))
  
  ggsave(volcano, filename = paste0("Volcano_Plots/", gsub(" ", "_", ct), "_Volcano.png"), width = 7.5, height = 7.5)
  
  # Heatmap of top genes
  top_genes <- degs %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
    arrange(p_val_adj) %>%
    slice(1:25) %>%
    pull(Gene)
  
  if (length(top_genes) > 1) {
    heat <- DoHeatmap(subset_obj, features = top_genes) +
      scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "RdBu"))
    
    ggsave(heat, filename = paste0("Heatmap_Plots/", gsub(" ", "_", ct), "_Heatmap.png"), width = 7.5, height = 7.5)
  }
}

# Additional: Baseline vs C2/C3 comparison within Cohort and Control
for (group in c("Cohort", "Control")) {
  group_obj <- subset(LBtox.no142, subset = Population == group & Timepoint %in% c("Baseline", "C2", "C3"))
  Idents(group_obj) <- "Timepoint"
  
  for (ct in celltypes) {
    message("Comparing timepoints for ", ct, " in ", group)
    ct_obj <- subset(group_obj, subset = predicted.celltype.l2 == ct)
    DefaultAssay(ct_obj) <- "RNA"
    
    # Timepoint comparison: Baseline vs C2
    degs_time <- FindMarkers(ct_obj, ident.1 = "Baseline", ident.2 = "C2", vars.to.regress = c("Pool", "Sex")) %>%
      rownames_to_column("Gene")
    
    outname <- paste0("DGE_Results/", gsub(" ", "_", ct), "_", group, "_Baseline_vs_C2.csv")
    write.csv(degs_time, outname, row.names = FALSE)
    
    plot_time <- EnhancedVolcano(degs_time,
                                 x = "avg_log2FC",
                                 y = "p_val_adj",
                                 lab = degs_time$Gene,
                                 FCcutoff = 0.5,
                                 pCutoff = 0.05,
                                 title = paste0(group, " Baseline vs C2: ", ct))
    
    ggsave(plot_time, filename = paste0("Volcano_Plots/", group, "_", gsub(" ", "_", ct), "_Baseline_vs_C2.png"), width = 7.5, height = 7.5)
  }
}