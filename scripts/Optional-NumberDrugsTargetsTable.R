
# make sure you have the extended CTL network open in Cytoscape!

my.drugs <- selectNodes("drug", by.col = "CTL.Type", preserve = FALSE)$nodes #easy way to collect node SUIDs by column value
selectFirstNeighbors()
selected <- getSelectedNodes()

createSubnetwork("selected", "all")

selectEdges("pp", by.col="interaction")
invertEdgeSelection()
e <- getSelectedEdges()
selectAllNodes()
n <- getSelectedNodes()
createSubnetwork(n, e)

RCy3::layoutNetwork()

RCy3::analyzeNetwork(directed=TRUE)

res <- getTableColumns(table="node", columns = c("display name", "GeneID", "CTL.Type", "Indegree", "Outdegree"))
colnames(res) <- c("Name", "ID", "Type", "NumDrugs", "NumTargets")

write.table(res, file="drug-target-info.tsv", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
