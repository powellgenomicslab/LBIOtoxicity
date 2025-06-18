# -----------------------------------------------
# EnRICH Polygenic Risk Score Analysis
# -----------------------------------------------

# Load libraries
library(tidyverse)
library(patchwork)
library(ggrepel)
library(extrafont)
library(survival)
library(survminer)
library(pROC)
library(ggpubr)

# Load fonts for Windows
extrafont::loadfonts(device = "win")

# Set working directory
setwd("/Users/jenniferlim/Dropbox (Personal)/PhD/EnRICH SNP IRAE Study/EnRICH SNP Substudy/SNP Data")

# -----------------------------------------------
# Import data
# -----------------------------------------------
enrich_polygenic <- read.csv("Polygenic_Risk_Enrich.csv")

# -----------------------------------------------
# Calculate GRS (unweighted sum of 6 SNPs post-LD pruning)
# -----------------------------------------------
snp_list6 <- c("rs1131017", "rs705704", "rs158366", "rs3759094", "rs6032664", "rs11171739")
enrich_polygenic$GRS_reduced6 <- rowSums(enrich_polygenic[, snp_list6], na.rm = TRUE)

# -----------------------------------------------
# Calculate PRS (weighted sum using chi² values)
# -----------------------------------------------

# Define SNP column range
snp_cols <- names(enrich_polygenic)[3:ncol(enrich_polygenic)]

# Ensure length consistency
LBIO_chi2 <- c(4.49, 7.64, 1.552, 0.29, 0.46, 4.1, 2.04, 0.42, 1.91, 0.97, 0, 1.37, 0, 1.44, 0.8, 0.97, 0.92, 0.04)
stopifnot(length(snp_cols) == length(LBIO_chi2))

# Calculate full PRS
enrich_polygenic$PRS <- as.numeric(as.matrix(enrich_polygenic[, snp_cols]) %*% LBIO_chi2)

# LD-pruned SNPs
snps_to_drop <- c("rs1131017", "rs877636")
snp_cols_pruned <- setdiff(snp_cols, snps_to_drop)
chi2_pruned <- LBIO_chi2[!snp_cols %in% snps_to_drop]

# Calculate LD-pruned PRS
enrich_polygenic$PRS_pruned <- as.numeric(as.matrix(enrich_polygenic[, snp_cols_pruned]) %*% chi2_pruned)

# -----------------------------------------------
# Statistical Testing
# -----------------------------------------------

# Boxplot with Wilcoxon test
ggplot(enrich_polygenic, aes(x = as.factor(Toxicity), y = PRS_pruned, fill = as.factor(Toxicity))) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("blue", "red"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "LD-Pruned PRS by Toxicity Status", x = "Toxicity", y = "PRS") +
  theme_minimal() +
  theme(legend.position = "none")

# Wilcoxon test
print(wilcox.test(PRS_pruned ~ Toxicity, data = enrich_polygenic))

# -----------------------------------------------
# ROC Curve + Threshold (Youden index)
# -----------------------------------------------
roc_obj <- roc(enrich_polygenic$Toxicity, enrich_polygenic$PRS_pruned)
threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden") %>% as.numeric()
print(threshold)

# Assign risk group
enrich_polygenic$Risk_group <- ifelse(enrich_polygenic$PRS_pruned >= threshold, "High Risk", "Low Risk")

# Chi-squared test
print(chisq.test(table(enrich_polygenic$Risk_group, enrich_polygenic$Toxicity)))

# -----------------------------------------------
# Plot: ROC Curve
# -----------------------------------------------
plot(roc_obj, col = "darkblue", lwd = 2, main = "ROC Curve for PRS Prediction")
legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 2)), col = "darkblue", lwd = 2)

# -----------------------------------------------
# Visualizations
# -----------------------------------------------

# Boxplot by Risk Group
ggplot(enrich_polygenic, aes(x = Risk_group, y = PRS_pruned, fill = Risk_group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "PRS by Risk Group", x = "Risk Group", y = "PRS") +
  theme_minimal() +
  theme(legend.position = "none")

# Boxplot with jitter for toxicity overlay
ggplot(enrich_polygenic, aes(x = Risk_group, y = PRS_pruned, fill = Risk_group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(aes(color = as.factor(Toxicity)), width = 0.2, height = 0.1) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "PRS by Risk Group with Toxicity Overlay", x = "Risk Group", y = "PRS") +
  theme_minimal()

# Histogram of PRS by toxicity
ggplot(enrich_polygenic, aes(x = PRS_pruned, fill = as.factor(Toxicity))) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("blue", "red"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "Histogram of PRS Scores by Toxicity", x = "PRS Score", y = "Count") +
  theme_minimal()

# Boxplot by toxicity with jitter
ggplot(enrich_polygenic, aes(x = as.factor(Toxicity), y = PRS_pruned)) +
  geom_boxplot(aes(fill = as.factor(Toxicity)), alpha = 0.7) +
  geom_jitter(aes(color = as.factor(Toxicity)), width = 0.2, height = 0.1) +
  scale_fill_manual(values = c("blue", "red"), labels = c("No Toxicity", "Toxicity")) +
  scale_color_manual(values = c("blue", "red"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "PRS by Toxicity Status", x = "Toxicity", y = "PRS Score") +
  theme_minimal()
