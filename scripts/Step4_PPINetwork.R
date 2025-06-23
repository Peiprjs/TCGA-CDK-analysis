# Script:       Step4_PPINetwork.R
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

# Make sure you ran the scripts for step 1-3 before this script

# ==================================================================
# PPI network creation with the stringApp for Cytoscape
# ==================================================================
# make sure Cytoscape is running
RCy3::cytoscapePing()

library(ggplot2)

# Check if WikiPathways app is installed
if(!"name: stringApp, version: 2.2.0, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("stringApp")
}

query <- format_csv(as.data.frame(degs$GeneName), col_names=F, escape = "double", eol =",")
commandsPOST(paste0('string protein query cutoff=0.7 newNetName="PPI network" query="',query,'" limit=0 species="Homo sapiens"'))

exportImage(paste0(out.folder,'PPI-network.png'), type='PNG', zoom=500) 

# ??? Question 10 - answer in document

# ==================================================================
# Analysis of hub nodes and high betweenness nodes
# ==================================================================

# Let's analyze the network and plot the degree distribution
RCy3::analyzeNetwork()
hist(RCy3::getTableColumns(columns = "Degree")$Degree, breaks=100, main = "PPI Degree distribution", xlab = "Degree")

filename <- paste0(out.folder,"PPI-degree-distribution.png")
png(filename , width = 500, height = 1200, res = 150)
dev.off()

# ??? Question 11 - answer in document

# Let's create a visualization based on the node degree (node color) and 
# betweenness (node size)

RCy3::createVisualStyle("centrality")
RCy3::setNodeLabelMapping("display name", style.name = "centrality")
control.points <- c (0, 0.3)
colors <-  c ('#FFFFFF', '#DD8855')
setNodeColorMapping("Degree", c(0,60), colors, style.name = "centrality", default.color = "#C0C0C0")
setNodeSizeMapping("BetweennessCentrality", table.column.values =  c (0, 0.5), sizes = c(50,100), mapping.type = "c", style.name = "centrality", default.size = 10)
RCy3::setVisualStyle("centrality")
toggleGraphicsDetails()

exportImage(paste0(out.folder,'PPI-centrality.png'), type='PNG', zoom=500) 

# ??? Question 12 - answer in document

# ==================================================================
# Visualization of data
# ==================================================================

RCy3::loadTableData(data=data, data.key.column = "GeneName", table = "node", table.key.column = "query term")
RCy3::createVisualStyle("log2FC vis")
RCy3::setNodeLabelMapping("display name", style.name = "log2FC vis")
control.points <- c (-3.0, 0.0, 3.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FC", control.points, colors, style.name = "log2FC vis", default.color = "#C0C0C0")
RCy3::setVisualStyle("log2FC vis")
RCy3::lockNodeDimensions("TRUE", "log2FC vis")

exportImage(paste0(out.folder,'PPI-with-data.png'), type='PNG', zoom=500) 

# ??? Question 13 - answer in document
