# Script:       Step3_PathwayVisualization.R
# Description:  In this script, we will explore the differential gene 
#               expression dataset comparing lung cancer vs. healthy tissue
#               samples. The RNA-sequencing dataset was retrieved from 
#               TCGA (The Cancer Genome Atlas) and pre-processed in R. 
#               Differential gene expression analysis was performed with the 
#               DESeq2 R-package. 
# Version: 3.0
# Last updated: 2025-06-24
# Author: mkutmon & peiprJS

# #############################################
# PATHWAY VISUALIZATION
# #############################################

# Open Pathway of interest - based on the res.wp.up.df and res.wp.down.df, you can
# select pathways of interest
# Find the associated pathway identifier
# https://www.wikipathways.org/browse/table.html
# Make sure you select the ID of the human pathway
# Example: Cell cycle pathway

#####################################################################
# CODING TASK - CHANGE PATHWAY ID TO PATHWAY OF INTEREST
#####################################################################
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',pw.id)) 

toggleGraphicsDetails()

# load the data into Cytoscape (columns get added at the bottom)
loadTableData(data, data.key.column = "GeneID", table.key.column = "Ensembl")

# visualize the log2FC as a node fill color gradient
RCy3::setNodeColorMapping(table.column = 'log2FC', mapping.type = 'c', table.column.values = c(-1,0,1), colors = paletteColorBrewerRdBu, default.color = '#FFFFFF', style.name = 'WikiPathways')

# Select significant genes and change border color
x <- RCy3::createColumnFilter('adj.P.Value', 'adj.P.Value', 0.05, "LESS_THAN")
RCy3::setNodeBorderColorBypass(x$nodes, new.colors = "#009900")
RCy3::setNodeBorderWidthBypass(x$nodes, new.sizes = 7)
RCy3::clearSelection()

exportImage(paste0(out.folder,'pathway-',pw.id,'.png'), type='PNG', zoom=500) 

# ??? Question 9 - answer in document

