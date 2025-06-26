# Script:       Step0_Preamble.R
# Description:  In this script, we import the packages required, and set up all needed environment variables.
# Version: 1.0
# Last updated: 2025-06-08
# Author: peiprJS
###############################################################################################################
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

###############################################################################################################

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
library(ggplot2)

###############################################################################################################

cytoscapePing()

if (!"name: WikiPathways, version: 3.3.10, status: Installed" %in% RCy3::getInstalledApps()) {RCy3::installApp("WikiPathways")}
if (!"name: stringApp, version: 2.2.0, status: Installed" %in% RCy3::getInstalledApps()) {RCy3::installApp("stringApp")}
if (!"name: clusterMaker2, version: 2.3.4, status: Installed" %in% RCy3::getInstalledApps()) {RCy3::installApp("clusterMaker2")}
if (!"name: CyTargetLinker, version: 4.1.0, status: Installed" %in% RCy3::getInstalledApps()) {RCy3::installApp("CyTargetLinker")}


###############################################################################################################

# /!\ Uncomment setwd line if you're running manually
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
out.folder <- "../output/"
dir.create(out.folder) 

###############################################################################################################

#setwd(dirname(str_replace(rstudioapi::getActiveDocumentContext()$path, "Step0_Preamble.R", "")))
data <- read_excel("../data/data-breast-cancer.xlsx")

###############################################################################################################

log2fc.cutoff <- 1
pvalue.cutoff <- 0.05

###############################################################################################################

pw.id <- "WP45" # Pathway of interest
interest_genes <- c('CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D')
interest_cluster <- ""
go.term <- "GO:0002768"

# Possible pathways of interest:
# WP179 - Cell Cycle
# WP5497 - Cyclin-dependent kinase 4/6 inhibitors in breast cancer
# WP5213 - Amino acid metabolism in triple-negative breast cancer cells
# WP4262 - Breast cancer pathway
# WP5211 - Glucose metabolism in triple-negative breast cancer cells
# WP1984 - Integrated breast cancer pathway
# WP5215 - Targeted agents in triple negative breast cancer
# WP5496 - Targeted therapy in breast cancer
# WP5353 - Macrophage-stimulating protein (MSP) signaling
