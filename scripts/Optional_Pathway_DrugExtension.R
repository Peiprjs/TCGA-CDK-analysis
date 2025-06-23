# Make sure you run at least step 1 + step 2 before this
# usually you would also do the pathway visualization step 3
# before doing the drug extension

# #############################################
# PATHWAY VISUALIZATION
# #############################################

# Check if Cytoscape is running
cytoscapePing()

# Check if WikiPathways app is installed
if(!"name: WikiPathways, version: 3.3.10, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("WikiPathways")
}

# Open Pathway of interest - based on the res.wp.up.df and res.wp.down.df, you can
# select pathways of interest
# Find the associated pathway identifier
# https://www.wikipathways.org/browse/table.html
# Make sure you select the ID of the human pathway
# Example: Cell cycle pathway

# difference to step 3 --> pathway will be imported as a network
# so all layout information is removed
pw.id <- "WP5497"
RCy3::commandsRun(paste0('wikipathways import-as-network id=',pw.id)) 

toggleGraphicsDetails()

# load the data into Cytoscape (columns get added at the bottom)
loadTableData(data, data.key.column = "GeneID", table.key.column = "Ensembl")

# visualize the log2FC as a node fill color gradient
RCy3::setNodeColorMapping(table.column = 'log2FC', mapping.type = 'c', table.column.values = c(-1,0,1), colors = paletteColorBrewerRdBu, default.color = '#FFFFFF', style.name = 'WikiPathways-As-Network')

# Select significant genes and change border color
x <- RCy3::createColumnFilter('P.Value', 'P.Value', 0.05, "LESS_THAN")
RCy3::setNodeBorderColorBypass(x$nodes, new.colors = "#009900")
RCy3::setNodeBorderWidthBypass(x$nodes, new.sizes = 7)
RCy3::clearSelection()


# ==================================================================
# Extend the pathway with known drug-target interactions from DrugBank
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



