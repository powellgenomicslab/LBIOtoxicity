#Repeat DGE for ALL other cell types

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
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_SNP_Metadata.RDS")
LBtox
colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Idents(Baseline.no142) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#baseline no LB142 and LB035 (female
Baseline.New <- subset(Baseline.no142, subset=Sex==c("Male"))

unique(LBtox$predicted.celltype.l2)

#Find All Markers
#repeat code for all cell types
Idents(LBtox.no142$predicted.celltype.l2) <- "Cohort"
ResponseMarkers <- FindAllMarkers(Baseline.no142, only.pos=TRUE, assay="RNA", min.cells.group=20 )
head(ResponseMarkers)
ResponseMarkersFilter <- filter(ResponseMarkers, p_val_adj < 0.05 & avg_log2FC >0.5)
dim(ResponseMarkersFilter)
plot11 <- DoHeatmap(subset(mostpbmc_celltypes$`Cytotoxic_CD8+_Tcell_S100B+`, downsample=100), features=ResponseMarkersFilter$gene, size=3)
ggsave(plot11, filename = "Heatmap_DGE_CytotoxicCD8.png", height = 4.5, width = 7.5)
mapply(find_markers, mostpbmc_celltypes, names(mostpbmc_celltypes))

#CD14 Mono (5505 cells)
Baseline.CD14.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD14 Mono"))
Idents(Baseline.CD14.no142) <- "Population"
head(Baseline.CD14.no142)

#scale data
Baseline.CD14.no142 <- ScaleData(Baseline.CD14.no142, verbose=FALSE)
colnames(Baseline.CD14.no142[[]])

#ensure rownames of gene expression matrix are gene names in uppercase

# Extract the gene expression matrix
gene_expression <- as.data.frame(as.matrix(Baseline.CD14.no142@assays$RNA@scale.data))
(as.matrix(Baseline.CD14.no142@assays$RNA@scale.data))
scaled_data <- GetAssayData(Baseline.CD14.no142, slot = "scale.data")
head(scaled_data)


# Extract the gene names from the row names of the gene expression matrix
gene_names <- rownames(gene_expression)

# Assign the gene names to the row names of the Seurat object
rownames(Baseline.CD14.no142@assays$RNA@scale.data) <- gene_names

#remove quotes
# Assuming your Seurat object is named Baseline.CD14.no142

# Convert row names to a character vector and remove quotes
row_names <- gsub("\"", "", as.character(rownames(Baseline.CD14.no142@assays$RNA@scale.data)))

# Assign modified row names back to the Seurat object
rownames(Baseline.CD14.no142@assays$RNA@scale.data) <- row_names
rownames(Baseline.CD14.no142@assays$RNA@scale.data)

#DEGS Analysis CD14 Mono
DefaultAssay(Baseline.CD14.no142) <- "RNA"
CD14MonoDEG <- FindMarkers(Baseline.CD14.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
head(CD14MonoDEG)


#convert row names to new column named Gene
CD14MonoDEG <- CD14MonoDEG %>% rownames_to_column(var = "Gene")
CD14MonoDEG$Gene <- paste0("\"", CD14MonoDEG$Gene, "\"")
head(CD14MonoDEG)

#extract top genes
top_genes <- CD14MonoDEG %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  pull(Gene)

  head(Baseline.CD14.no142)
head(top_genes)

top_genes <- toupper(top_genes)
rownames(Baseline.CD14.no142@assays$RNA@scale.data) <- toupper(rownames(Baseline.CD14.no142@assays$RNA@scale.data))
existing_genes <- intersect(top_genes, rownames(Baseline.CD14.no142@assays$RNA@data))
head(existing_genes)

#Create HeatMap
plot1<- DoHeatmap(Baseline.CD14.no142, features = existing_genes) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "RdBu"))

  ggsave(plot1, filename = "Heatmap CD14 Mono DEGs.png", height = 7.5, width = 7.5)

  #repeated with hierarchical clustering
 

# clustered heatmap
scaled_data <- Baseline.CD14.no142@assays$RNA@scale.data[existing_genes, ]
  
  # Generate the heatmap with hierarchical clustering
 png("heatmap CD14 Mono.png", width = 800, height = 800)
  
  plot2 <- pheatmap(scaled_data,
    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    scale = "row" )

    ggsave(plot2, filename = "Heatmap CD14 Mono DEGs.png", height = 7.5, width = 7.5)
  
  
  #can you please create individual heat maps on a loop for differential gene expression for each cell type in predicted.celltype.l2 ? 


write.csv(DEGS_CD14Mono, "CD14 Mono Differentially expressed Genes.CSV")    

#Volcano Plot CD14 Mono
plot1 <- EnhancedVolcano(DEGS_CD14Mono, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD14Mono$Gene, FCcutoff=0.5)
  ggsave(plot1,filename="Volcano Plot CD14 Monocytes3.png", height = 7.5, width = 10)


####please perform FindMarkers for all cell types between Cohort and Control (as above) and create heat maps for each cell type with more than 100 cells

# Loop through all cell subsets
cell_subsets <- c("CD14 Mono", "CD8 TEM", "NK", "CD16 Mono", "CD4 Naive", "Treg", "CD4 CTL", "CD4 TEM")
plots <- list()

for (subset in cell_subsets) {
  # Subset the Seurat object
  subset_obj <- subset(x = Baseline.no142, subset = predicted.celltype.l2 == subset)
  Idents(subset_obj) <- "Population"
  
  # Perform FindMarkers
  DefaultAssay(subset_obj) <- "RNA"
  DEGS_subset <- FindMarkers(subset_obj, ident.1 = "Cohort", ident.2 = "Control", vars.to.regress = c("Pool", "Sex"))
  DEGS_subset <- DEGS_subset %>% tibble::rownames_to_column(var = "Gene")
  
  # Save the results
  file_name <- paste0(subset, "_Differentially_expressed_Genes.CSV")
  write.csv(DEGS_subset, file_name)
  
  # Create Volcano Plot
  plot <- EnhancedVolcano(DEGS_subset, x = "avg_log2FC", y = "p_val_adj", lab = DEGS_subset$Gene, FCcutoff = 0.1)  # Subset the Seurat object based on cell type
  subset_seurat <- subset(x = Baseline.no142, subset = predicted.celltype.l2 == subset)
  Idents(subset_seurat) <- "Population"
  
  # Perform differential expression analysis
  DefaultAssay(subset_seurat) <- "RNA"
  DEGS_subset <- FindMarkers(subset_seurat, ident.1 = "Cohort", ident.2 = "Control", vars.to.regress = c("Pool", "Sex"))
  
  # Filter significant genes
  DEGS_subset <- filter(DEGS_subset, p_val_adj < 0.05 & avg_log2FC > 0.5)
  
  # Get the top genes
  top_genes <- rownames(DEGS_subset)
  
  # Create heatmap
  plot <- DoHeatmap(subset_seurat, features = top_genes)
  
  # Save the plot
  filename <- paste0("Heatmap_", gsub(" ", "", subset), ".png")
  ggsave(plot, filename = filename, height = 4.5, width = 7.5)
  
  # Store the plot in the list
  plots[[subset]] <- plot
}
# Create a directory to save the volcano plots
dir.create("Volcano_Plots")

# Save each volcano plot in the directory
ggsave(plot3, filename = "Volcano_Plots/DEGS_Bnaive.png", height = 4.5, width = 7.5)
ggsave(plot4, filename = "Volcano_Plots/DEGS_NK.png", height = 4.5, width = 7.5)
ggsave(plot5, filename = "Volcano_Plots/DEGS_CD16.png", height = 4.5, width = 7.5)
ggsave(plot6, filename = "Volcano_Plots/DEGS_CD4.png", height = 4.5, width = 7.5)
ggsave(plot7, filename = "Volcano_Plots/DEGS_Treg.png", height = 4.5, width = 7.5)
ggsave(plot8, filename = "Volcano_Plots/DEGS_CD4CTL.png", height = 4.5, width = 7.5)
ggsave(plot9, filename = "Volcano_Plots/DEGS_CD4TEM.png", height = 4.5, width = 7.5)


#create heatmaps for each cell subset DEG
# Loop through all cell subsets
cell_subsets <- c("CD14 Mono", "CD8 TEM", "NK", "CD16 Mono", "CD4 Naive", "Treg", "CD4 CTL", "CD4 TEM")
plots <- list()
# Create a directory to save the heatmaps
dir.create("Heatmap_Plots")

# Loop through all cell subsets
for (subset in cell_subsets) {
  # Subset the Seurat object for the specific cell subset
  subset_obj <- subset(x = Baseline.no142, subset = predicted.celltype.l2 == subset)
  
  # Perform DEG analysis for the cell subset
  DEG <- FindMarkers(subset_obj, ident.1 = "Cohort", ident.2 = "Control", vars.to.regress = c("Pool", "Sex"))
  
  # Filter DEGs based on significance
  DEG_filtered <- filter(DEG, p_val_adj < 0.05 & avg_log2FC > 0.5)
  
  # Select the top genes
  top_genes <- rownames(DEG_filtered)[1:100]
  
  # Create the heatmap
  heatmap_plot <- DoHeatmap(subset_obj, features = top_genes)
  
  # Save the heatmap plot in the directory
  ggsave(heatmap_plot, filename = paste0("Heatmap_Plots/Heatmap_", subset, ".png"), height = 4.5, width = 7.5)
  
  # Create the volcano plot
  volcano_plot <- EnhancedVolcano(DEG, x = "avg_log2FC", y = "p_val_adj", lab = DEG$Gene, FCcutoff = 0.5)
  
  # Save the volcano plot in the directory
  ggsave(volcano_plot, filename = paste0("Volcano_Plots/Volcano_", subset, ".png"), height = 4.5, width = 7.5)
}

# Access the plots using plots$subset_name

  *******

#Volcano Plot B naive
plot3 <- EnhancedVolcano(DEGS_Bnaive, x="avg_log2FC", y="p_val_adj", lab=DEGS_Bnaive$Gene, FCcutoff=0.1)
  ggsave(plot3,filename="Volcano Plot B naive DEGs.png", height = 7.5, width = 7.5)

#NK (2435 cells)
Baseline.NK.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("NK"))
Idents(Baseline.NK.no142) <- "Population"

#DEGS Analysis NK cells
DefaultAssay(Baseline.NK.no142) <- "RNA"
DEGS_NK <- FindMarkers(Baseline.NK.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_NK <- DEGS_NK %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_NK)

write.csv(DEGS_NK, "NK cells Differentially expressed Genes.CSV")   

#Volcano Plot NK cells
plot4 <- EnhancedVolcano(DEGS_NK, x="avg_log2FC", y="p_val_adj", lab=DEGS_NK$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot NK cells DEGs.png", height = 7.5, width = 7.5)


#CD16 Mono (847 cells)
Baseline.CD16.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD16 Mono"))
Idents(Baseline.CD16.no142) <- "Population"

#DEGS Analysis NK cells
DefaultAssay(Baseline.CD16.no142) <- "RNA"
DEGS_CD16 <- FindMarkers(Baseline.CD16.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD16 <- DEGS_CD16 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD16)

write.csv(DEGS_CD16, "CD16 Monocyte cells Differentially expressed Genes.CSV")   

#Volcano Plot CD16 cells
plot5 <- EnhancedVolcano(DEGS_CD16, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD16$Gene, FCcutoff=0.1)
  ggsave(plot5,filename="Volcano Plot CD16 mono cells DEGs.png", height = 7.5, width = 7.5)

#CD4 naive (415 cells)
Baseline.CD4.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD4 Naive"))
Idents(Baseline.CD4.no142) <- "Population"

#DEGS Analysis NK cells
DefaultAssay(Baseline.CD4.no142) <- "RNA"
DEGS_CD4 <- FindMarkers(Baseline.CD4.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4 <- DEGS_CD4 %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CD4)

write.csv(DEGS_CD4, "CD4 Naive cells Differentially expressed Genes.CSV")   

#Volcano Plot CD16 cells
plot5 <- EnhancedVolcano(DEGS_CD4, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4$Gene, FCcutoff=0.1)
  ggsave(plot5,filename="Volcano Plot CD4 mono cells DEGs.png", height = 7.5, width = 7.5)

#Treg (300 cells)
Baseline.Treg.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("Treg"))
Idents(Baseline.Treg.no142) <- "Population"

#DEGS Analysis T reg
DefaultAssay(Baseline.Treg.no142) <- "RNA"
DEGS_Treg <- FindMarkers(Baseline.Treg.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_Treg <- DEGS_Treg %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_Treg)

write.csv(DEGS_Treg, "Treg cells Differentially expressed Genes.CSV")   

#Volcano Plot T reg cells
plot5 <- EnhancedVolcano(DEGS_Treg, x="avg_log2FC", y="p_val_adj", lab=DEGS_Treg$Gene, FCcutoff=0.1)
  ggsave(plot5,filename="Volcano Plot Treg cells DEGs.png", height = 7.5, width = 7.5)

#CD4 CTL (322 cells)
Baseline.CTL.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD4 CTL"))
Idents(Baseline.CTL.no142) <- "Population"

#DEGS Analysis T reg
DefaultAssay(Baseline.CTL.no142) <- "RNA"
DEGS_CTL <- FindMarkers(Baseline.CTL.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CTL <- DEGS_CTL %>%
            tibble::rownames_to_column(var="Gene")

head(DEGS_CTL)

write.csv(DEGS_CTL, "CD4 CTL cells Differentially expressed Genes.CSV")   

TopGenes <- rownames(DEGS_CTL[order(DEGS_CTL$p_val_adj, decreasing=FALSE)[1:25],])
TopGenes

#annotating and ordering cells for dittoheatmap
plot1 <- dittoHeatmap(LBtox, TopGenes, annot.by=c("Population"))

#Volcano Plot CD4 CTL
plot5 <- EnhancedVolcano(DEGS_CTL, x="avg_log2FC", y="p_val_adj", lab=DEGS_CTL$Gene, title="CD4 CTL gene expression", pCutoff=10e-2)
  ggsave(plot5,filename="Volcano Plot CD4 CTL cells DEGs 4.png", height = 7.5, width = 7.5)

#CD4 TEM cells
Baseline.CD4TEM.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD4 TEM"))
Idents(Baseline.CD4TEM.no142) <- "Population"
DefaultAssay(Baseline.CD4TEM.no142) <- "RNA"
DEGS_CD4TEM<- FindMarkers(Baseline.CD4TEM.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4TEM <- DEGS_CD4TEM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4TEM)

write.csv(DEGS_CD4TEM, "CD4 TEM Differentially expressed Genes Baseline cohort vs control.CSV")   

#Volcano Plot CD4 TEM
plot4 <- EnhancedVolcano(DEGS_CD4TEM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4TEM$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot CD4 TEM cells DEGs Baseline.png", height = 7.5, width = 7.5)

#CD4 naive cells
Baseline.CD4naive.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD4 Naive"))
Idents(Baseline.CD4naive.no142) <- "Population"
DefaultAssay(Baseline.CD4naive.no142) <- "RNA"
DEGS_CD4naive<- FindMarkers(Baseline.CD4naive.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD4naive <- DEGS_CD4naive %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD4naive)

write.csv(DEGS_CD4naive, "CD4 Naive Differentially expressed Genes Baseline cohort vs control.CSV")   

#Volcano Plot CD4 naive cells
plot5 <- EnhancedVolcano(DEGS_CD4naive, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD4naive$Gene, FCcutoff=0.1)
  ggsave(plot5,filename="Volcano Plot CD4 Naive cells DEGs Baseline.png", height = 7.5, width = 7.5)


#CD8 TCM (217 cells)
Baseline.CD8TCM.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD8 TCM"))
Idents(Baseline.CD8TCM.no142) <- "Population"
DefaultAssay(Baseline.CD8TCM.no142) <- "RNA"
DEGS_CD8TCM <- FindMarkers(Baseline.CD8TCM.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8TCM <- DEGS_CD8TCM %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD8TCM)

write.csv(DEGS_CD8TCM, "CD8 TCM Differentially expressed Genes Baseline cohort vs control.CSV")   

#Volcano Plot CD8 TCM cells
plot1 <- EnhancedVolcano(DEGS_CD8TCM, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8TCM$Gene, FCcutoff=0.1)
  ggsave(plot1,filename="Volcano Plot CD8 TCM cells DEGs Baseline.png", height = 7.5, width = 7.5)

#CD8 Naive (177 cells)
Baseline.CD8naive.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("CD8 Naive"))
Idents(Baseline.CD8naive.no142) <- "Population"
DefaultAssay(Baseline.CD8naive.no142) <- "RNA"
DEGS_CD8naive<- FindMarkers(Baseline.CD8naive.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_CD8naive <- DEGS_CD8naive %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_CD8naive)

write.csv(DEGS_CD8naive, "CD8 Naive Differentially expressed Genes Baseline cohort vs control.CSV")   

#Volcano Plot CD8 naive cells
plot5 <- EnhancedVolcano(DEGS_CD8naive, x="avg_log2FC", y="p_val_adj", lab=DEGS_CD8naive$Gene, FCcutoff=0.1)
  ggsave(plot5,filename="Volcano Plot CD8 Naive cells DEGs Baseline.png", height = 7.5, width = 7.5)

#gdT cells
Baseline.gdT.no142 <- subset(x=Baseline.no142, subset=predicted.celltype.l2==c("gdT"))
Idents(Baseline.gdT.no142) <- "Population"
DefaultAssay(Baseline.gdT.no142) <- "RNA"
DEGS_gdT <- FindMarkers(Baseline.gdT.no142, ident.1="Cohort", ident.2="Control", vars.to.regress=c("Pool", "Sex"))
DEGS_gdT <- DEGS_gdT %>%
            tibble::rownames_to_column(var="Gene")
head(DEGS_gdT)

write.csv(DEGS_gdT, "gdT Differentially expressed Genes Baseline cohort vs control.CSV")   

#Volcano Plot gdT naive cells
plot4 <- EnhancedVolcano(DEGS_gdT, x="avg_log2FC", y="p_val_adj", lab=DEGS_gdT$Gene, FCcutoff=0.1)
  ggsave(plot4,filename="Volcano Plot gdT cells DEGs Baseline.png", height = 7.5, width = 7.5)


