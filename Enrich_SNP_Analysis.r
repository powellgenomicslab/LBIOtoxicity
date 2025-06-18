# -----------------------------------------------
# EnRICH SNP and HLA Analysis
# -----------------------------------------------

# Load libraries
library(tidyverse)
library(ggplot2)
library(janitor)
library(patchwork)
library(ggpubr)
library(survival)
library(survminer)
library(poolr)
extrafont::loadfonts(device = "win")

# Set working directory
setwd("/Users/jenniferlim/Dropbox (Personal)/PhD/EnRICH SNP IRAE Study/EnRICH SNP Substudy/Data Folder")

# -----------------------------------------------
# Read in EnRICH Data
# -----------------------------------------------
enrich_data <- read_csv("Enrich_IO_Data_1K.csv", col_types = cols(
  .default = col_double(),
  Biobank_ID = col_character(),
  Treatment = col_character(),
  Treatment_Start_Date = col_date(format = "%d/%m/%Y"),
  G3_IRAE = col_logical(),
  IRAE_Type = col_character(),
  Date_IRAE_Onset = col_date(format = "%d/%m/%Y"),
  PD_Date = col_date(format = "%d/%m/%Y"),
  Death_Date = col_date(format = "%d/%m/%Y"),
  .delim = ","
))

# -----------------------------------------------
# Derive Responder Variable
# -----------------------------------------------
enrich_data <- enrich_data %>%
  mutate(Responder = if_else(CR | PR, "Responder", "Non-Responder"))

# -----------------------------------------------
# Chi-square: Responder vs G3 IRAE
# -----------------------------------------------
tbl_responder <- table(enrich_data$G3_IRAE, enrich_data$Responder)
chi_responder <- chisq.test(tbl_responder)
print(chi_responder)

# Stacked Bar Plot
ggplot(enrich_data, aes(x = Responder, fill = G3_IRAE)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("lightblue", "salmon"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "G3 IRAE by Responder Status", x = "Responder Status", y = "Proportion", fill = "G3 IRAE") +
  theme_minimal() +
  theme(legend.position = "top")

# -----------------------------------------------
# Chi-square: Toxicity vs SNPs (looped)
# -----------------------------------------------
snp_list <- c(
  "rs1131017", "rs6032664", "rs4899554", "rs11171739", "rs158366", "rs12047412",
  "rs34629529", "rs3759094", "rs2442752", "rs3869120", "rs2395471", "rs116812356",
  "rs3749946", "rs805262", "rs9268914", "rs2844463", "rs877636", "rs2249741"
)

snp_results <- map_df(snp_list, function(snp) {
  tbl <- table(enrich_data$G3_IRAE, enrich_data[[snp]])
  test <- suppressWarnings(chisq.test(tbl))
  tibble(
    SNP = snp,
    p_value = test$p.value,
    statistic = test$statistic,
    df = test$parameter
  )
})

print(snp_results %>% arrange(p_value))

# -----------------------------------------------
# HLA Logistic Regression Analysis
# -----------------------------------------------
hla_data <- read_csv("HLA_Short_Data.csv") %>%
  mutate(G3_IRAE = as.factor(G3_IRAE))

hla_model <- function(haplotype_var) {
  model <- glm(G3_IRAE ~ .data[[haplotype_var]], data = hla_data, family = binomial)
  summary_model <- summary(model)
  tibble(
    Haplotype = names(coef(model)),
    Odds_Ratio = exp(coef(model)),
    P_Value = coef(summary_model)[, "Pr(>|z|)"]
  ) %>%
    filter(P_Value < 0.05)
}

significant_HLAC <- hla_model("HLAC_Short_Haplotype")
significant_HLADRB1 <- hla_model("HLADRB1_Short_Haplotype")

print(significant_HLAC)
print(significant_HLADRB1)

# -----------------------------------------------
# Multiple Testing Correction
# -----------------------------------------------
p_val <- read_csv("Combined_Pvalue_LBIO_Enrich.csv")

p_val <- p_val %>%
  mutate(
    LBIO_P_Val_Bonferroni = p.adjust(LBIO_P_Val, method = "bonferroni"),
    ENRICH_P_Val_Bonferroni = p.adjust(ENRICH_P_Val, method = "bonferroni"),
    LBIO_P_Val_BH = p.adjust(LBIO_P_Val, method = "BH"),
    ENRICH_P_Val_BH = p.adjust(ENRICH_P_Val, method = "BH")
  )

print(p_val)
