# Make sure you run at least step 1 + step 2 before this
# the input for this optional step is a list of pathway names
# that you find in res.wp.df table
# the idea is that you find all genes that are associated
# with the pathways, make one PPI network, visualize the data 
# and find drug targets

# so instead of finding the module from the large network,
# you define the module based on the pathways of interest

# #############################################
# GET GO TERM
# #############################################

library(tidyr)

# Check if Cytoscape is running
cytoscapePing()

# look up pathway names in the enrichment result tables - you can add more
# in the vector if you need to
# these pathways are from the lung cancer practical - replace with pathways
# from your breast cancer analysis!
pathways <- pathways.stepeight
geneSets <- res.wp@geneSets
genes <- geneSets[names(geneSets) %in% pathways]

# long format first: Gene - Pathway
gene_pathway_df <- do.call(rbind, lapply(names(genes), function(pw) {
  data.frame(Pathway = pw, Gene = genes[[pw]], stringsAsFactors = FALSE)
}))

pathway.genes <- gene_pathway_df %>%
  dplyr::mutate(Present = 1) %>%
  tidyr::pivot_wider(names_from = Pathway, values_from = Present, values_fill = 0)

# add NumPathways column: count how many pathways each gene appears in
pathway.genes <- pathway.genes %>%
  mutate(NumPathways = rowSums(across(all_of(pathways))))

# reorder columns: Gene, N_Pathways, then pathways
pathway.genes <- pathway.genes %>%
  select(Gene, NumPathways, all_of(pathways))

# view result
head(pathway.genes)

query <- format_csv(as.data.frame(pathway.genes$Gene), col_names=F, escape = "double", eol =",")
commandsPOST(paste0('string protein query cutoff=0.9 newNetName="PPI network" query="',query,'" limit=0 species="Homo sapiens"'))

analyzeNetwork()
RCy3::loadTableData(data=data, data.key.column = "GeneID", table = "node", table.key.column = "query term")
RCy3::loadTableData(data=pathway.genes, data.key.column = "Gene", table = "node", table.key.column = "query term")

# VISUALIZATION 1 
# fill color = log2FC gradient
# green border = significance

RCy3::createVisualStyle("log2FC vis")
RCy3::setNodeLabelMapping("display name", style.name = "log2FC vis")
control.points <- c (-3.0, 0.0, 3.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FC", control.points, colors, style.name = "log2FC vis", default.color = "#C0C0C0")
RCy3::setEdgeColorDefault("#CCCCCC", style.name = "log2FC vis")
RCy3::setVisualStyle("log2FC vis")
RCy3::lockNodeDimensions("TRUE", "log2FC vis")

# Select significant genes and change border color
x <- RCy3::createColumnFilter('adj.P.Value', 'adj.P.Value', 0.05, "LESS_THAN")
RCy3::setNodeBorderColorBypass(x$nodes, new.colors = "#009900")
RCy3::setNodeBorderWidthBypass(x$nodes, new.sizes = 7)
RCy3::clearSelection

exportImage(paste0(out.folder,'selected-pathways-PPI.png'), type='PNG', zoom=500) #.png; use zoom or width args to increase size/resolution


# ==================================================================
# Extend the process PPI with known drug-target interactions from DrugBank
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

exportImage(paste0(out.folder,'selected-pathways-PPI-with-drugs.png'), type='PNG', zoom=500) #.png; use zoom or width args to increase size/resolution

# ==================================================================
# SAVING CYTOSCAPE SESSION
# ==================================================================

saveSession(paste0(out.folder,'selected-pathways-PPI.cys'))

## OPTIONAL - FILTER network to only include drugs + targets

selectNodes(my.drugs)
selectFirstNeighbors()
selected <- getSelectedNodes()

createSubnetwork("selected", "all")
RCy3::layoutNetwork()
