# Make sure you run at least step 1 + step 2 before this
# the input for this optional step is one of the GO processes
# that you find in res.go.up.df or res.go.down.df
# the idea is that you find all genes that are associated
# with the process, make a PPI network, visualize the data 
# and find drug targets

# so instead of finding the module from the large network,
# you define the module based on a process of interest

# #############################################
# GET GO TERM
# #############################################

# Check if Cytoscape is running
cytoscapePing()

# Check if WikiPathways app is installed
if(!"name: stringApp, version: 2.1.1, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("stringApp")
}

# example "DNA replication"
# look up GO ID in the enrichment result tables
go.term <- "GO:0006260"
genes <- res.go@geneSets[go.term]


query <- format_csv(as.data.frame(genes[1]), col_names=F, escape = "double", eol =",")
commandsRun(paste0('string protein query cutoff=0.9 newNetName="PPI network" query="',query,'" limit=0 species="Homo sapiens"'))

analyzeNetwork()
RCy3::loadTableData(data=data, data.key.column = "GeneName", table = "node", table.key.column = "query term")
RCy3::createVisualStyle("log2FC vis")
RCy3::setNodeLabelMapping("display name", style.name = "log2FC vis")
control.points <- c (-3.0, 0.0, 3.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FC", control.points, colors, style.name = "log2FC vis", default.color = "#C0C0C0")
RCy3::setVisualStyle("log2FC vis")
RCy3::lockNodeDimensions("TRUE", "log2FC vis")

# Select significant genes and change border color
x <- RCy3::createColumnFilter('P.Value', 'P.Value', 0.05, "LESS_THAN")
RCy3::setNodeBorderColorBypass(x$nodes, new.colors = "#009900")
RCy3::setNodeBorderWidthBypass(x$nodes, new.sizes = 7)
RCy3::clearSelection()

# ==================================================================
# Extend the process PPI with known drug-target interactions from DrugBank
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

# ==================================================================
# SAVING CYTOSCAPE SESSION
# ==================================================================

saveSession(paste0(out.folder,'lung-cancer-example.cys'))



## OPTIONAL - FILTER network to only include drugs + targets

selectNodes(my.drugs)
selectFirstNeighbors()
selected <- getSelectedNodes()

createSubnetwork("selected", "all")
