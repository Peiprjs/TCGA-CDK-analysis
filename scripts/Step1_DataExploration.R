# Script:       Step1_DataExploration.R
# Description:  In this script, we will explore the differential gene 
#               expression dataset comparing lung cancer vs. healthy tissue
#               samples. The RNA-sequencing dataset was retrieved from 
#               TCGA (The Cancer Genome Atlas) and pre-processed in R. 
#               Differential gene expression analysis was performed with the 
#               DESeq2 R-package. 
# Version: 2.0
# Last updated: 2025-06-08
# Author: mkutmon

# #############################################
# R INSTRUCTIONS
# #############################################

# * Lines that start with a # are comments
# * You can run a code line by placing the cursor in the line and clicking 
#   CTRL/Command + Enter
# * In between you will see the ??? Questions which refers to a question 
#   in the question document provided on Canvas.


# #############################################
# R SETUP
# #############################################

# Here we install and load all required packages. 
# This might take a few minutes
if (!("BiocManager" %in% installed.packages())) { install.packages("BiocManager", update=FALSE) }
if (!("rstudioapi" %in% installed.packages())) { BiocManager::install("rstudioapi", update=FALSE) }
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db", update=FALSE) }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr", update=FALSE) }
if (!("EnhancedVolcano" %in% installed.packages())) { BiocManager::install("EnhancedVolcano", update=FALSE) }
if (!("readxl" %in% installed.packages())) { BiocManager::install("readxl", update=FALSE) }
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler", update=FALSE) }
if (!("enrichplot" %in% installed.packages())) { BiocManager::install("enrichplot", update=FALSE) }
if (!("Rgraphviz" %in% installed.packages())) { BiocManager::install("Rgraphviz", update=FALSE) }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3", update=FALSE) }
if (!("msigdbr" %in% installed.packages())) { BiocManager::install("msigdbr",update=FALSE) }
if (!("RColorBrewer" %in% installed.packages())) { BiocManager::install("RColorBrewer",update=FALSE) }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr",update=FALSE) }
if (!("rWikiPathways" %in% installed.packages())) { BiocManager::install("rWikiPathways",update=FALSE) }

library(rstudioapi)
library(org.Hs.eg.db)
library(dplyr)
library(EnhancedVolcano)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(Rgraphviz)
library(RCy3)
library(msigdbr)
library(RColorBrewer)
library(readr)
library(rWikiPathways)
library(stringr)

# We will set the working directory to the location where the current 
# script is located. This way, we can use relative file path locations. 
setwd(dirname(str_replace(rstudioapi::getActiveDocumentContext()$path, "Step1_DataExploration.R", "")
))

# We will create an output folder where all figures and files will be stored
out.folder <- "output/"
dir.create(out.folder)

# #############################################
# IMPORT DATA
# #############################################

# Next we will import the dataset
data <- read_excel("data/data-breast-cancer.xlsx")

# ??? Question 1 - answer in document


# #############################################
# DIFFERENTIALLY EXPRESSED GENES
# #############################################

# Let's see how many genes are differentially expressed in our dataset

log2fc.cutoff <- 1
pvalue.cutoff <- 0.05
degs <- data[abs(data$log2FC) > log2fc.cutoff & data$adj.P.Value < pvalue.cutoff,]

# let's write out the table with all differentially expressed genes
write.table(degs, file=paste0(out.folder,"degs.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

# ??? Question 2 - answer in document


# Based on the code in line 68, can you adapt the code to select only up- or
# down-regulated genes? 

#####################################################################
# CODING TASK - REPLACE ... WITH ADAPTED CODE BASED ON LINE 85
#####################################################################

log2fc.cutoff <- -1
genes.down <- data[data$log2FC < log2fc.cutoff & data$adj.P.Value < pvalue.cutoff,]
write.table(degs, file=paste0(out.folder,"degs_down.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

log2fc.cutoff <- 1
genes.up <- data[data$log2FC > log2fc.cutoff & data$adj.P.Value < pvalue.cutoff,]
write.table(degs, file=paste0(out.folder,"degs_up.tsv"), row.names = FALSE, sep="\t", quote = FALSE)


# ??? Question 3 - answer in document

# #############################################
# VOLCANO PLOT
# #############################################

# Let's create a Volcano plot to get a better understand of the intensity and
# direction of the changed genes

EnhancedVolcano(data, title = paste0("Breast cancer vs. Healthy (",nrow(degs), " DEGs)"), lab = data$GeneName, x = "log2FC", y = "adj.P.Value", pCutoff = pvalue.cutoff, FCcutoff = log2fc.cutoff, labSize = 3, xlim = c(-15,15), ylim=c(0,4))

# the code below saves the figure in a file in our output folder
filename <- paste0(out.folder,"volcano-plot.png")
png(filename , width = 2000, height = 1500, res = 150)
EnhancedVolcano(data, title = paste0("Breast cancer vs. Healthy (",nrow(degs), " DEGs)"), lab = data$GeneName, x = "log2FC", y = "adj.P.Value", pCutoff = pvalue.cutoff, FCcutoff = log2fc.cutoff, labSize = 3, xlim = c(-15,15), ylim=c(0,4))
dev.off()

# ??? Question 4 - answer in document

