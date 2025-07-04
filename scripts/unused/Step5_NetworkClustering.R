# Script:       Step5_NetworkClustering.R
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

# Make sure you ran the scripts for step 1-4 before this script

# ==================================================================
# PPI network creation with the stringApp for Cytoscape
# ==================================================================
# make sure Cytoscape is running
RCy3::cytoscapePing()

# Check if clusterMaker2 app is installed
if(!"name: clusterMaker2, version: 2.3.4, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("clusterMaker2")
}

# let's use community clustering to find clusters/modules in the network
commandsRun(paste0('cluster glay createGroups=TRUE'))

# ??? Question 14 - answer in document

# ==================================================================
# Visualize clusters in large network
# ==================================================================

cluster_data <- getTableColumns('node', columns=c('__glayCluster'))
unique_clusters <- sort(unique(cluster_data$`__glayCluster`))

cluster_counts <- table(cluster_data$`__glayCluster`)
large_clusters <- names(cluster_counts[cluster_counts > 15])
colors <- rainbow(length(large_clusters))  

RCy3::createVisualStyle("clustering")
RCy3::setNodeLabelMapping("display name", style.name = "centrality")
setNodeColorMapping(table.column = "__glayCluster", table.column.values = large_clusters, mapping.type = "d", colors = colors, style.name = "clustering", default.color = "#D9D9D9")
RCy3::setVisualStyle("clustering")
toggleGraphicsDetails()
exportImage(paste0(out.folder,'clustered-network.png'), type='PNG', zoom=0, width=4096) #.png; use zoom or width args to increase size/resolution

# ==================================================================
# Function of up-regulated cluster
# ==================================================================

RCy3::setVisualStyle("log2FC vis")

# In the node table in Cytoscape, scroll to the right to see the "__glayCluster" column
# Select a node in the up-regulated cluster and check which number it has
# Add the number below with the cluster number you are interested in

# let's first create a subnetwork with only the up-regulated cluster
####################################################
# CODING TASK - ADD CLUSTER NUMBER INSTEAD OF ...
####################################################
cluster <- 1339
nodes.cluster <- RCy3::createColumnFilter('__glayCluster', '__glayCluster', cluster, predicate = "IS")
RCy3::createSubnetwork(nodes = nodes.cluster$nodes, nodes.by.col = "shared name", subnetwork.name = paste0("PPI-cluster-", cluster))
exportImage(paste0(out.folder,'cluster-',cluster,'.png'), type='PNG', zoom=500) #.png; use zoom or width args to increase size/resolution

# Functional analysis
# You could then run an enrichment analysis as we did in step 2:
genes <- RCy3::getTableColumns(columns = c("query term"))$`query term`
res.go <- clusterProfiler::enrichGO(genes, OrgDb = "org.Hs.eg.db", 
                                       keyType="SYMBOL", universe=data$GeneName, ont = "BP", 
                                       pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                                       minGSSize = 20, maxGSSize = 400)
res.go.df <- as.data.frame(res.go)


# ??? Question 15 - answer in document


# ==================================================================
# Extend the cluster with known drug-target interactions from DrugBank
# ==================================================================

# Check if CyTargetLinker app is installed
if(!"name: CyTargetLinker, version: 4.1.0, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("CyTargetLinker")
}

unzip(system.file("extdata","drugbank-5.1.0.xgmml.zip", package="rWikiPathways"), exdir = getwd())
drugbank <- file.path(getwd(), "drugbank-5.1.0.xgmml")

commandsRun(paste0('cytargetlinker extend idAttribute="GeneID" linkSetFiles="', drugbank, '"') )
commandsRun('cytargetlinker applyLayout network="current"')

my.drugs <- selectNodes("drug", by.col = "CTL.Type", preserve = FALSE)$nodes #easy way to collect node SUIDs by column value
clearSelection()
setNodeColorBypass(my.drugs, "#DD99FF")
setNodeShapeBypass(my.drugs, "hexagon")

drug.labels <- getTableColumns(columns=c("SUID","CTL.label"))
drug.labels <- na.omit(drug.labels)
mapply(function(x,y) setNodeLabelBypass(x,y), drug.labels$SUID, drug.labels$CTL.label)

# Try different layouts (e.g. yFiles organic layout) if nodes are overlapping too much 
# Cytoscape > Layout menu!

exportImage(paste0(out.folder,'cluster-',cluster,'-with-drugs.png'), type='PNG', zoom=500) #.png; use zoom or width args to increase size/resolution

# ??? Question 16 - answer in document


# ==================================================================
# SAVING CYTOSCAPE SESSION
# ==================================================================

saveSession(paste0(out.folder,'lung-cancer-example.cys'))

