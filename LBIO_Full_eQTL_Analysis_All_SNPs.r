#eQTL Analysis loop all SNPs and genes

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


#import pairs CSV
# Option A: base‐R
pairs <- read.csv("SNP_with_Gene.csv", stringsAsFactors = FALSE)

# Option B: data.table’s fread (if you prefer speed)
library(data.table)
pairs <- fread("SNP_with_Gene.csv", stringsAsFactors = FALSE)

meta.dt <- as.data.table(LBtox@meta.data)
Idents(LBtox) <- "MajoritySinglet_Individual_Assignment"

# 2) Define the function (do NOT change snp_col/gene_name here):
run_one_eQTL <- function(snp_col, gene_name) {
  # … your code for χ² + boxplots …
  #
  dt_unique <- unique(LBtox@meta.data[, c("Population", snp_col)])
  tbl2      <- table(dt_unique$Population, dt_unique[[snp_col]])
  chi2_res  <- chisq.test(tbl2)
  
  out_csv <- paste0("eQTL_results/chi2_", snp_col, "_", gene_name, ".csv")
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(
    data.frame(
      SNP     = snp_col,
      Gene    = gene_name,
      Chi2    = chi2_res$statistic,
      df      = chi2_res$parameter,
      p_value = chi2_res$p.value
    ),
    out_csv, row.names = FALSE
  )
  
  LBtox_baseline <- subset(LBtox, subset = Timepoint == "Baseline")
  meta_base      <- LBtox_baseline@meta.data
  
  genes      <- c("RPS26","RPS9","FOS","JUN","HLA-C","HLA-DRB5")
  expr_matrix <- LBtox_baseline@assays$SCT@counts[genes, ]
  for (g in genes) {
    meta_base[[g]] <- as.numeric(expr_matrix[g, ])
  }
  
  df <- meta_base %>%
    as.data.frame() %>%
    select(
      celltype = predicted.celltype.l2,
      genotype = all_of(snp_col),
      expr     = all_of(gene_name)
    ) %>%
    filter(!is.na(genotype) & !is.na(celltype))
  df$expr <- as.numeric(df$expr)
  
  alleles      <- sort(unique(df$genotype))
  palette_vals <- setNames(
    RColorBrewer::brewer.pal(n = length(alleles), name = "Set2")[seq_along(alleles)],
    alleles
  )
  
  celltypes <- unique(df$celltype)
  plot_list <- lapply(celltypes, function(ct) {
    sub <- df %>% filter(celltype == ct)
    if (nrow(sub) < 3) return(NULL)
    p <- ggplot(sub, aes(x = genotype, y = expr, color = genotype)) +
      geom_boxplot(width = 0.5, outlier.shape = NA, coef = 0) +
      geom_quasirandom(width = 0.1, alpha = 0.6) +
      scale_color_manual(values = palette_vals) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = paste(gene_name, "by", snp_col, "in", ct),
        x     = snp_col,
        y     = paste0(gene_name, " expression")
      )
    if (length(alleles) >= 2 && all(alleles %in% sub$genotype)) {
      maxy        <- max(sub$expr, na.rm = TRUE)
      comparisons <- combn(alleles, 2, simplify = FALSE)
      y_pos       <- seq(from = maxy * 1.05, by = maxy * 0.05, length.out = length(comparisons))
      p <- p + geom_signif(
        comparisons      = comparisons,
        map_signif_level = TRUE,
        y_position       = y_pos
      )
    }
    return(p)
  })
  
  pdf_dir <- file.path("eQTL_results", paste0("plots_", snp_col, "_", gene_name))
  dir.create(pdf_dir, showWarnings = FALSE)
  for (i in seq_along(plot_list)) {
    if (is.null(plot_list[[i]])) next
    ct    <- celltypes[i]
    fname <- file.path(
      pdf_dir,
      paste0(
        gene_name, "_", snp_col, "_",
        gsub("[[:space:][:punct:]]+", "_", ct),
        ".pdf"
      )
    )
    ggsave(
      filename = fname,
      plot     = plot_list[[i]],
      device   = "pdf",
      width    = 10,
      height   = 8,
      units    = "cm"
    )
  }
  return(NULL)
}

# --------------------------------------------------------
# 3) NOW you call the function with real SNP/gene names:
# --------------------------------------------------------

# Example 1: single, interactive call
run_one_eQTL("rs1131017", "RPS26")

# --------------------------------------------------------
# 4) If you have a 'pairs' data.frame or CSV, loop over it:
# --------------------------------------------------------

# Suppose pairs.csv has columns: SNP,Gene
pairs <- read.csv("SNP_with_Gene.csv", stringsAsFactors = FALSE)

for (i in seq_len(nrow(pairs))) {
  snp_i  <- pairs$SNP[i]
  gene_i <- pairs$Gene[i]
  message("Running eQTL for ", snp_i, " → ", gene_i, " …")
  run_one_eQTL(snp_i, gene_i)
}