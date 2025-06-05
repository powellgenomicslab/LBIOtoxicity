#single cell proportion test

#install packages
library(Seurat) 
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(magrittr)
library(jsonlite)
library(data.table)
library(Azimuth)
library(sctransform)
library(scProportionTest)
library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(speckle)
devtools::install_github("rpolicastro/scProportionTest")

#set wd
#set working directory
my_dir <- "/directflow/SCCGGroupShare/projects/jenli3/LBTOX"
setwd(my_dir)

#Read RDS
LBtox <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/LBTox_Merged_Metadata_IncSex.RDS")
LBtox
Idents(LBtox) <- "Population"

colnames(LBtox[[]])
Idents(LBtox) <- "predicted.celltype.l2"
DefaultAssay(LBtox) <- "RNA"

#baseline no LB142
LBtox.no142 <- subset(x=LBtox, subset= MajoritySinglet_Individual_Assignment == c("LB142-C1"), invert=TRUE)
Baseline.no142 <- subset(x=LBtox.no142, subset=Timepoint==c("Baseline"))
Idents(Baseline.no142) <- "predicted.celltype.l2"
DefaultAssay(Baseline) <- "RNA"

#Summary of cell clusters and numbers - baseline
colnames(LBtox[[]])

#Metadata table
LBtox.metadatasummary <- Baseline.no142@meta.data %>% 
group_by(MajoritySinglet_Individual_Assignment, predicted.celltype.l2, Population, Timepoint) %>%
summarise(counts=n()) %>%
group_by(MajoritySinglet_Individual_Assignment, Population, Timepoint) %>%
mutate(Percentage=counts/sum(counts)*100)
head(LBtox.metadatasummary)

#prop test vignette
prop_test <- sc_utils(LBtox)
prop_test <- permutation_test(prop_test, cluster_identity="predicted.celltype.l2", sample_1="Cohort", sample_2="Control", sample_identity="Population")
plot1 <- permutation_plot(prop_test)
ggsave(plot1)


#Speckle Test
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(version = "3.17")

BiocManager::install("speckle")
library(speckle)

#From Walter

library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(patchwork)
​
​


​
# Create plot -------------------------------------------------------------
​
# Define a colour palette for each cohort
Population.palette <- c("Cohort" = "#F58800",
                    "Control" = "#2c7e8c")
​
# Plot
plot.list <- lapply(unique(LBtox.metadatasummary$predicted.celltype.l2), function(i){
    p1 <- LBtox.metadatasummary %>%
        filter(predicted.celltype.l2  ==  i) %>%
        ggplot(
            mapping = aes(
                x = Population,
                y = Percentage,
                color = as.factor(Population)
            ),
            fill = as.factor(Population)
        ) +
        geom_boxplot(
            aes(color = as.factor(Population)),
            width = 0.45,
            position = position_dodge(0.8),
            outlier.shape = NA,
            coef = 0
        ) +
        coord_cartesian(ylim=c(0,40)
        ) +
        geom_quasirandom(dodge.width = .8,
                         width = 0.04,
                         alpha = 0.75) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = Population.palette,
                            name = "") +
        theme(
              axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            ) +
        xlab(i) +
        theme(axis.title.x = element_text(size=7)) +
        ylab("Percentage of cells (%)") +
        ylim(c(min(LBtox.metadatasummary$Percentage),
          1.1*max(LBtox.metadatasummary$Percentage)))
        
        # Apply the Wilcoxon Rank Sum Test between tumour pairs
        sigFunc = function(x){
            if(x < 0.001){"***"} 
            else if(x < 0.01){"**"}
            else if(x < 0.05){"*"}
            else{NA}}
    
        p1 <- p1 + geom_signif(
            comparisons = list(c("Cohort", "Control")),
            test = "wilcox.test",
            map_signif_level  =  sigFunc,
            colour  =  "black",
            y_position = 30
        )
    
        # Remove y axis except for the first plot
    if (i != "B intermediate") {
        p1 <- p1 + theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank()
        )
    }
    
    return(p1)
})
​
ggsave(filename = "Cell Proportion Significance v6.pdf",
       plot = wrap_plots(plot.list, ncol = 4, guides="collect"),
       device = "pdf", width = 28, height = 45, units = "cm")


#calculate p value
result <- sapply(unique(LBtox.metadatasummary$predicted.celltype.l2), function(i){
    x <- filter(LBtox.metadatasummary, predicted.celltype.l2==i)
    group1 <- x$Percentage[x$Population=="Cohort"]
    group2 <- x$Percentage[x$Population=="Control"]
    result <- wilcox.test(group1, group2)
    return(result$p.value)
})


# Print the result
print(result)

#stacked barchart for ALL individuals and cell proportions

metadata.summary <- readRDS("/directflow/SCCGGroupShare/projects/jenli3/LBTOX/SeuratObjects/Baseline_Metadata.RDS") %>%
    group_by(Timepoint, predicted.celltype.l2) %>%
    summarise(counts = n()) %>%
    group_by(Timepoint) %>%
    mutate(percentage = counts/sum(counts)*100)
​
​
​
# Plot stacked barchart ---------------------------------------------------


p1 <- C2.metadatasummary %>%
    ggplot(aes(fill=predicted.celltype.l2, y=Percentage, x=MajoritySinglet_Individual_Assignment)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab("Individuals") +
    ylab("Percentage") +
    theme_bw() +
    ggtitle("Cell type percentages: Individuals") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_discrete("predicted.celltype.l2")
​
ggsave(filename = "Individual cell proportions bar chart C2.pdf",
       plot = p1,
       device = "pdf", width = 20, height = 20, units = "cm")


#Cycle 2 cohort vs control 
#Metadata table
C2.metadatasummary <- C2.no142@meta.data %>% 
group_by(MajoritySinglet_Individual_Assignment, predicted.celltype.l2, Population, Timepoint) %>%
summarise(counts=n()) %>%
group_by(MajoritySinglet_Individual_Assignment, Population, Timepoint) %>%
mutate(Percentage=counts/sum(counts)*100)
head(C2.metadatasummary)

# Define a colour palette for each cohort
Population.palette <- c("Cohort" = "#F58800",
                    "Control" = "#2c7e8c")
​
# Plot
C2.plot.list <- lapply(unique(C2.metadatasummary$predicted.celltype.l2), function(i){
    p2 <- C2.metadatasummary %>%
        filter(predicted.celltype.l2  ==  i) %>%
        ggplot(
            mapping = aes(
                x = Population,
                y = Percentage,
                color = as.factor(Population)
            ),
            fill = as.factor(Population)
        ) +
        geom_boxplot(
            aes(color = as.factor(Population)),
            width = 0.45,
            position = position_dodge(0.8),
            outlier.shape = NA,
            coef = 0
        ) +
        geom_quasirandom(dodge.width = .8,
                         width = 0.04,
                         alpha = 0.75) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        scale_colour_manual(values = Population.palette,
                            name = "") +
        theme(
              axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            ) +
        xlab(i) +
        theme(axis.title.x = element_text(size=7)) +
        ylab("Percentage of cells (%)") +
        ylim(c(min(C2.metadatasummary$Percentage),
          1.1*max(C2.metadatasummary$Percentage)))
        
        # Apply the Wilcoxon Rank Sum Test between tumour pairs
        sigFunc = function(x){
            if(x < 0.001){"***"} 
            else if(x < 0.01){"**"}
            else if(x < 0.05){"*"}
            else{NA}}
    
        p2 <- p2 + geom_signif(
            comparisons = list(c("Cohort", "Control")),
            test = "wilcox.test",
            map_signif_level  =  sigFunc,
            colour  =  "black",
            y_position = 30
        )
    
    
    return(p2)
})
​
ggsave(filename = "Cell Proportion Significance C2.pdf",
       plot = wrap_plots(C2.plot.list, ncol = 4, guides="collect"),
       device = "pdf", width = 28, height = 45, units = "cm")


#calculate p value
result <- sapply(unique(C2.metadatasummary$predicted.celltype.l2), function(i){
    x <- filter(C2.metadatasummary, predicted.celltype.l2==i)
    group1 <- x$Percentage[x$Population=="Cohort"]
    group2 <- x$Percentage[x$Population=="Control"]
    result <- wilcox.test(group1, group2)
    return(result$p.value)
})