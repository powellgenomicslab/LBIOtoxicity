# ---------------------------------------------
# HLA Analysis: LBIO and EnRICH Cohorts
# ---------------------------------------------

# Load libraries
if (!require("BIGDAWG")) install.packages("BIGDAWG")
library(BIGDAWG)
library(tidyverse)

# ----------------------------
# Set working directories
# ----------------------------
setwd("/Users/jenniferlim/Dropbox (Personal)/PhD/PBMC Study/IRAE Substudy/10x Data/HLA Analysis/")

# ----------------------------
# Load and prepare LBIO data
# ----------------------------
HLA <- read.csv("LBIO_HLA.csv")
HLA$Toxicity <- factor(HLA$Toxicity, levels = c(0, 1), labels = c("Control", "Case"))

# ----------------------------
# Run BIGDAWG: LBIO cohort
# ----------------------------
results_LBIO <- BIGDAWG(
  Data = HLA,
  Run.Tests = "H",
  Loci.Set = list(c("C", "DRB1")),
  Missing = 0,
  Output = FALSE,
  Return = TRUE
)

# ----------------------------
# Extract and save frequency table
# ----------------------------
freq_table <- results_LBIO$H$Set1$freq
colnames(freq_table) <- c("Haplotype", "Count_Control", "Count_Case")
write.csv(freq_table, "LBIO_HLA_Frequency_Table.csv", row.names = FALSE)

# ----------------------------
# Summarize results with ORs and CIs
# ----------------------------
or_list <- results_LBIO$H$Set1$OR
summary_df <- map_dfr(or_list, ~{
  data.frame(
    OR = as.numeric(.x["OR"]),
    CI.L = as.numeric(.x["CI.L"]),
    CI.U = as.numeric(.x["CI.U"]),
    stringsAsFactors = FALSE
  )
}, .id = "Haplotype")

summary_df <- bind_cols(freq_table, summary_df)
write.csv(summary_df, "LBIO_HLA_haplotype_summary.csv", row.names = FALSE)

# ----------------------------
# Fisher test helper function
# ----------------------------
run_fisher_analysis <- function(freq, N_case, N_ctrl, outfile) {
  total_case_haps <- N_case * 2
  total_ctrl_haps <- N_ctrl * 2

  freq$Raw_Case_Count <- round(freq$Count_Case * total_case_haps)
  freq$Raw_Ctrl_Count <- round(freq$Count_Control * total_ctrl_haps)

  freq$Fisher_p <- map_dbl(seq_len(nrow(freq)), function(i) {
    a <- freq$Raw_Case_Count[i]
    b <- total_case_haps - a
    c <- freq$Raw_Ctrl_Count[i]
    d <- total_ctrl_haps - c
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  })

  freq$Fisher_adj_p <- p.adjust(freq$Fisher_p, method = "BH")
  write.csv(freq, outfile, row.names = FALSE)
  return(freq)
}

# ----------------------------
# EnRICH HLA analysis
# ----------------------------
setwd("/Users/jenniferlim/Dropbox (Personal)/PhD/EnRICH SNP IRAE Study/EnRICH SNP Substudy/Data Folder/HLA analysis/")
Enrich_HLA <- read.csv("Enrich_HLA_Pivoted.csv")

# Define samples
N_case <- 32
N_ctrl <- 185 - N_case

# ----------------------------
# Run BIGDAWG on full HLA loci
# ----------------------------
results_enrich_full <- BIGDAWG(
  Data = Enrich_HLA,
  Run.Tests = "H",
  Loci.Set = list(c("A", "B", "C", "DQB1", "DRB1")),
  Missing = 0,
  Output = FALSE,
  Return = TRUE
)

freq_full <- results_enrich_full$H$Set1$freq
colnames(freq_full) <- c("Haplotype", "Count_Control", "Count_Case")
write.csv(freq_full, "Enrich_HLA_Frequency_Table.csv", row.names = FALSE)

# Fisher's test + adjustment
freq_with_p <- run_fisher_analysis(freq_full, N_case, N_ctrl, "Enrich_HLA_Frequency_Table_with_pvalues.csv")
sig_haplos <- filter(freq_with_p, Fisher_adj_p < 0.05)
print(sig_haplos)

# ----------------------------
# Repeat for HLA-C~DRB1
# ----------------------------
results_enrich_CDRB1 <- BIGDAWG(
  Data = Enrich_HLA,
  Run.Tests = "H",
  Loci.Set = list(c("C", "DRB1")),
  Missing = 0,
  Output = FALSE,
  Return = TRUE
)

freq_cdrb1 <- results_enrich_CDRB1$H$Set1$freq
colnames(freq_cdrb1) <- c("Haplotype", "Count_Control", "Count_Case")
write.csv(freq_cdrb1, "Enrich_HLA_C_DRB1_Frequency_Table.csv", row.names = FALSE)

freq_cdrb1_with_p <- run_fisher_analysis(freq_cdrb1, N_case, N_ctrl, "Enrich_HLA_C_DRB1_with_pvalues.csv")
sig_haplos_cdrb1 <- filter(freq_cdrb1_with_p, Fisher_adj_p < 0.05)
print(sig_haplos_cdrb1)

# ----------------------------
# Repeat for HLA-C only (single locus)
# ----------------------------
results_enrich_C <- BIGDAWG(
  Data = Enrich_HLA,
  Run.Tests = "L",
  Loci.Set = list(c("C")),
  Missing = 0,
  Output = FALSE,
  Return = TRUE
)

freq_c <- results_enrich_C$H$Set1$freq
colnames(freq_c) <- c("Haplotype", "Count_Control", "Count_Case")
write.csv(freq_c, "Enrich_HLA_C_Frequency_Table.csv", row.names = FALSE)

freq_c_with_p <- run_fisher_analysis(freq_c, N_case, N_ctrl, "Enrich_HLA_C_with_pvalues.csv")
sig_haplos_c <- filter(freq_c_with_p, Fisher_adj_p < 0.05)
print(sig_haplos_c)
