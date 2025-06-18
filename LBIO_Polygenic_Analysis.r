# Polygenic Risk Score (PRS) Analysis - LBIO Pilot

# Load required libraries
library(ggplot2)
library(Seurat)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(magrittr)
library(stringr)
library(janitor)
library(survival)
library(survminer)
library(readr)
library(ggsci)
library(ggpubr)
library(stats)

# Set working directory
my_dir <- "/Users/jenniferlim/Dropbox (Personal)/PhD/PBMC Study/IRAE Substudy/10x Data/SNP analysis"
setwd(my_dir)

# Import data
polygenic_data <- read.csv("LBIO_Pilot_Polygenic v2.csv")

# Define chi-squared values and SNP columns
chi2_values <- c(4.49, 7.64, 3.13, 1.552, 0.29, 0.46, 4.1, 2.04, 0.42,
                 1.91, 0.97, 0, 1.37, 0, 1.44, 0.8, 0.97, 0.92, 0.04)
snp_cols <- names(polygenic_data)[3:ncol(polygenic_data)]
stopifnot(length(chi2_values) == length(snp_cols))

# Calculate PRS
polygenic_data$PRS <- as.numeric(as.matrix(polygenic_data[, snp_cols]) %*% chi2_values)

# PRS Histogram
ggplot(polygenic_data, aes(x = PRS)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of PRS", x = "PRS", y = "Count") +
  theme_minimal()

# PRS Boxplot by Toxicity
boxplot(PRS ~ Toxicity, data = polygenic_data,
        xlab = "Toxicity", ylab = "Polygenic Risk Score",
        main = "PRS by Toxicity Status", col = c("lightblue", "salmon"))

# Logistic regression
model <- glm(Toxicity ~ PRS, data = polygenic_data, family = binomial)
polygenic_data$predicted_prob <- predict(model, type = "response")

# PRS vs Toxicity with logistic curve
plot(polygenic_data$PRS, polygenic_data$Toxicity,
     xlab = "Polygenic Risk Score", ylab = "Toxicity (actual)",
     main = "PRS vs Toxicity", pch = 16, col = "gray")
curve(predict(model, data.frame(PRS = x), type = "response"),
      add = TRUE, col = "blue", lwd = 2)

# ROC Curve
roc_obj <- pROC::roc(polygenic_data$Toxicity, polygenic_data$predicted_prob)
plot(roc_obj, col = "darkblue", lwd = 2, main = "ROC Curve for PRS Model")
legend("bottomright", legend = paste("AUC =", round(pROC::auc(roc_obj), 2)), col = "darkblue", lwd = 2)

# LD Pruning
snps_to_drop2 <- c("rs1131017", "rs877636", "rs34629529", "rs3869120", "rs116812356")
snp_cols_pruned <- setdiff(snp_cols, snps_to_drop2)
chi2_pruned <- chi2_values[!snp_cols %in% snps_to_drop2]
polygenic_data$PRS_pruned <- as.numeric(as.matrix(polygenic_data[, snp_cols_pruned]) %*% chi2_pruned)

# Histogram of pruned PRS
ggplot(polygenic_data, aes(x = PRS_pruned)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Pruned PRS", x = "Pruned PRS", y = "Count") +
  theme_minimal()

# Boxplot of pruned PRS by Toxicity
boxplot(PRS_pruned ~ Toxicity, data = polygenic_data,
        xlab = "Toxicity", ylab = "LD-pruned PRS",
        main = "LD-Pruned PRS by Toxicity", col = c("lightblue", "salmon"))
stripchart(PRS_pruned ~ Toxicity, data = polygenic_data,
           method = "jitter", pch = 16, vertical = TRUE, col = "darkgray", add = TRUE)

# Wilcoxon test
wilcox.test(PRS_pruned ~ Toxicity, data = polygenic_data)

# Logistic regression on pruned PRS
model_pruned <- glm(Toxicity ~ PRS_pruned, data = polygenic_data, family = binomial)
polygenic_data$predicted_prob_pruned <- predict(model_pruned, type = "response")

# ROC Curve for pruned model
roc_obj_pruned <- pROC::roc(polygenic_data$Toxicity, polygenic_data$predicted_prob_pruned)
plot(roc_obj_pruned, col = "darkblue", lwd = 2, main = "ROC Curve for LD-Pruned PRS")
legend("bottomright", legend = paste("AUC =", round(pROC::auc(roc_obj_pruned), 2)), col = "darkblue", lwd = 2)

# Risk stratification
threshold <- as.numeric(pROC::coords(roc_obj_pruned, "best", ret = "threshold", best.method = "youden"))
polygenic_data$Risk_group <- ifelse(polygenic_data$PRS_pruned >= threshold, "High Risk", "Low Risk")

# Risk group vs toxicity
print(threshold)
table(polygenic_data$Risk_group, polygenic_data$Toxicity)
chisq.test(table(polygenic_data$Risk_group, polygenic_data$Toxicity))

# Boxplot of PRS by Risk Group
ggplot(polygenic_data, aes(x = Risk_group, y = PRS_pruned, fill = Risk_group)) +
  geom_boxplot() +
  labs(title = "LD-Pruned PRS by Risk Group", x = "Risk Group", y = "LD-Pruned PRS") +
  scale_fill_manual(values = c("Low Risk" = "blue", "High Risk" = "red")) +
  theme_minimal() +
  theme(legend.position = "none")

# Histogram by Toxicity
ggplot(polygenic_data, aes(x = PRS_pruned, fill = as.factor(Toxicity))) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("blue", "red"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "Histogram of PRS by Toxicity", x = "PRS Score", y = "Count") +
  theme_minimal()

# Boxplot of PRS with jitter by toxicity
ggplot(polygenic_data, aes(x = Risk_group, y = PRS_pruned, fill = Risk_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = as.factor(Toxicity)), width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("Low Risk" = "blue", "High Risk" = "red")) +
  scale_color_manual(values = c("0" = "gray30", "1" = "red"), labels = c("No Toxicity", "Toxicity")) +
  labs(title = "PRS by Risk Group and Toxicity", x = "Risk Group", y = "LD-Pruned PRS") +
  theme_minimal() +
  theme(legend.position = "right")