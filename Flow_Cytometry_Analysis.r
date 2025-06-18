# -------------------------------------------
# MFI Calculations for c-FOS and c-JUN
# -------------------------------------------

# Load libraries
library(tidyverse)

# Set working directory
setwd("/Users/jenniferlim/Dropbox (Personal)/PhD/PBMC Study/IRAE Substudy/Validation Experiment Planning/Formal Cytek Experiment/Data")

# -------------------------------------------
# Define reusable plotting + test function
# -------------------------------------------
run_mfi_analysis <- function(data, markers, title_prefix = "Geometric Mean") {
  results <- list()

  for (marker in markers) {
    # Run t-test
    test_result <- t.test(reformulate("Toxicity", marker), data = data, var.equal = TRUE)
    print(paste0("T-test for ", marker))
    print(test_result)

    # Boxplot
    boxplot(
      data[[marker]] ~ data$Toxicity,
      main = paste(marker, title_prefix),
      xlab = "Toxicity",
      ylab = "Geometric Mean",
      col = c("lightblue", "lightcoral")
    )

    # Store result
    results[[marker]] <- test_result
  }

  return(results)
}

# -------------------------------------------
# Analysis: c-FOS dataset
# -------------------------------------------

cfos_data <- read.csv("CFOS_Geometric_Mean.csv")

# Preview
head(cfos_data)

# Define markers to analyze
cfos_markers <- c("CD4_TCM", "CD4_TEM", "CD8_TCM", "B_NAIVE", "MONOCYTES")

# Run analysis
cfos_results <- run_mfi_analysis(cfos_data, cfos_markers, "Geometric Mean (c-FOS)")

# -------------------------------------------
# Analysis: c-JUN dataset
# -------------------------------------------

cjun_data <- read.csv("CJUN_Geometric_Mean.csv")
cjun_markers <- c("CD4_TCM", "CD4_TEM", "CD8_TCM", "T_REG")

# Run analysis
cjun_results <- run_mfi_analysis(cjun_data, cjun_markers, "Geometric Mean (c-JUN)")

pdf("cFOS_MFI_Boxplots.pdf", width = 6, height = 5)
cfos_results <- run_mfi_analysis(cfos_data, cfos_markers, "Geometric Mean (c-FOS)")
dev.off()

pdf("cJUN_MFI_Boxplots.pdf", width = 6, height = 5)
cjun_results <- run_mfi_analysis(cjun_data, cjun_markers, "Geometric Mean (c-JUN)")
dev.off()
