library(dplyr)

interest_genes <- c('CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D')

cluster <- RCy3::getTableColumns(columns = "GeneName,__glayCluster")
cluster %>% filter_at(vars(GeneName), any_vars(. %in% c(interest_genes)))                                                                                                     

RCy3::createSubnetwork(nodes = nodes.cluster$nodes, nodes.by.col = "shared name", subnetwork.name = paste0("PPI-cluster-", cluster))

