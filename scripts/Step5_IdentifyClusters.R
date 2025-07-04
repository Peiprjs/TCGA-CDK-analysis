# Script:       Step5_NetworkClustering.R
# Description:  In this script, we will explore the differential gene 
#               expression dataset comparing lung cancer vs. healthy tissue
#               samples. The RNA-sequencing dataset was retrieved from 
#               TCGA (The Cancer Genome Atlas) and pre-processed in R. 
#               Differential gene expression analysis was performed with the 
#               DESeq2 R-package. 
# Version: 3.0
# Last updated: 2025-06-26
# Author: mkutmon & peiprjs

# Make sure you ran the scripts for step 1-4 before this script

# ==================================================================
# PPI network creation with the stringApp for Cytoscape
# ==================================================================

# make sure Cytoscape is running
RCy3::cytoscapePing()

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
exportImage(paste0(out.folder,'clustered-network.svg'), type='SVG', zoom=500) #.png; use zoom or width args to increase size/resolution

# ==================================================================
  # Function of interest cluster
# ==================================================================

RCy3::setVisualStyle("log2FC vis")

genes.full <- RCy3::getTableColumns(columns = "GeneName,__glayCluster")
genes.interest <- genes.full %>% filter_at(vars(GeneName), any_vars(. %in% c(interest_genes)))    

genes.interest

cluster <- as.integer(interest_cluster)
nodes.cluster <- RCy3::createColumnFilter('__glayCluster', '__glayCluster', cluster, predicate = "IS")


# if (interest_cluster == "") {
#   nodes.cluster <- RCy3::createColumnFilter('__glayCluster', '__glayCluster', unique(genes.interest[,2]), predicate = "IS")
# } else {
#   cluster <- interest_cluster
#   nodes.cluster <- RCy3::createColumnFilter('__glayCluster', '__glayCluster', cluster, predicate = "IS")
# }
  
RCy3::createSubnetwork(nodes = nodes.cluster$nodes, nodes.by.col = "shared name", subnetwork.name = paste0("PPI-cluster-Interest"))
exportImage(paste0(out.folder,'cluster.svg'), type='SVG', zoom=500)

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

exportImage(paste0(out.folder,cluster,'cluster-interest-with-drugs.svg'), type='SVG', zoom=500) #.png; use zoom or width args to increase size/resolution

# ==================================================================
# SAVING CYTOSCAPE SESSION
# ==================================================================

#saveSession(paste0(out.folder,'breast-cancer-example.cys'))