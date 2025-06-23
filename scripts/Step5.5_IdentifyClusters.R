library(dplyr)

interest_genes <- c('CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D')

genes.full <- RCy3::getTableColumns(columns = "GeneName,__glayCluster")
genes.interest <- genes.full %>% filter_at(vars(GeneName), any_vars(. %in% c(interest_genes)))                                                                                                     

nodes.cluster <- RCy3::createColumnFilter('__glayCluster', '__glayCluster', unique(genes.interest[,2]), predicate = "IS")

RCy3::createSubnetwork(nodes = nodes.cluster$nodes, nodes.by.col = "shared name", subnetwork.name = paste0("PPI-cluster-Interest"))
exportImage(paste0(out.folder,'cluster-',cluster,'.png'), type='PNG', zoom=500)