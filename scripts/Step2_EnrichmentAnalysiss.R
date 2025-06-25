# Script:       Step2_EnrichmentAnalysis.R
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
# GENE ONTOLOGY ENRICHMENT ANALYSIS
# #############################################

# We will first explore what kind of biological processes are affected by 
# performing a Gene Ontology enrichment analysis. 

# We will start by looking at processes that are up-regulated
res.go.up <- clusterProfiler::enrichGO(genes.up$GeneID, OrgDb = "org.Hs.eg.db", 
                   keyType="ENSEMBL", universe=data$GeneID, ont = "BP", 
                   pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                   minGSSize = 20, maxGSSize = 400)
res.go.up.df <- as.data.frame(res.go.up)
write.table(res.go.up.df, file=paste0(out.folder,"go-up.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

res.go.down <- clusterProfiler::enrichGO(genes.down$GeneID, OrgDb = "org.Hs.eg.db", 
                                       keyType="ENSEMBL", universe=data$GeneID, ont = "BP", 
                                       pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                                       minGSSize = 20, maxGSSize = 400)
res.go.down.df <- as.data.frame(res.go.down)
write.table(res.go.down.df, file=paste0(out.folder,"go-down.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Both lists are very long and difficult to interpret. Let's see if the 
# treeplots (discussed in workshop 1) help with the interpretation:
# we first calculate the similarity between the processes
res.go.up.sim <- enrichplot::pairwise_termsim(res.go.up)
res.go.down.sim <- enrichplot::pairwise_termsim(res.go.down)

# then we visualize them in a treeplot (this will only look at the top 100 terms)
treeplot(res.go.up.sim, label_format = 0.5, showCategory = 100, cluster.params = list(n = 15, label_words_n = 0))
treeplot(res.go.down.sim, label_format = 0.5, showCategory = 100, cluster.params = list(n = 15, label_words_n = 0))

# we will also save files that has a nicely readable figure
filename <- paste0(out.folder,"GO_Upregulated_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(treeplot(res.go.up.sim, label_format = 0.5, showCategory = 100, cluster.params = list(n = 15, label_words_n = 0)))
dev.off()

filename <- paste0(out.folder,"GO_Downregulated_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(treeplot(res.go.down.sim, label_format = 0.5, showCategory = 100, cluster.params = list(n = 15, label_words_n = 0)))
dev.off()

# #############################################
# PATHWAY ENRICHMENT ANALYSIS
# #############################################

# We will now explore what kind of pathways are affected by 
# performing a pathway enrichment analysis. 

# Let's retrieve information about the human pathways in WikiPathways
genesets.wp <- msigdbr(species = "Homo sapiens", subcollection = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, ensembl_gene)

# Let's first look for altered pathways (up- and down-regulated genes together)
# Most often when looking at pathways, we look for up- and down-regulated genes together
res.wp <- clusterProfiler::enricher(degs$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 400)
res.wp.df <- as.data.frame(res.wp)
write.table(res.wp.df, file=paste0(out.folder,"wp-degs.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# For exploration, we can also check if there are pathways that are enriched for up- or down-regulated genes separately
res.wp.up <- clusterProfiler::enricher(genes.up$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 400)
res.wp.up.df <- as.data.frame(res.wp.up)
write.table(res.wp.up.df, file=paste0(out.folder,"wp-up.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

res.wp.down <- clusterProfiler::enricher(genes.down$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 400)
res.wp.down.df <- as.data.frame(res.wp.down)
write.table(res.wp.down.df, file=paste0(out.folder,"wp-down.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Let's see if the treeplots (discussed in workshop 1) help with the interpretation:
# we first calculate the similarity between the processes
res.wp.sim <- enrichplot::pairwise_termsim(res.wp)
res.wp.up.sim <- enrichplot::pairwise_termsim(res.wp.up)
res.wp.down.sim <- enrichplot::pairwise_termsim(res.wp.down)

# then we visualize them in a treeplot
treeplot(res.wp.sim, label_format = 0.5, showCategory = 100, cluster.params = list(label_words_n = 0))
treeplot(res.wp.up.sim, label_format = 0.5, showCategory = 100, cluster.params = list(label_words_n = 0))
treeplot(res.wp.down.sim, label_format = 0.5, showCategory = 100, cluster.params = list(label_words_n = 0))

# we will also save files that has a nicely readable figure

filename <- paste0(out.folder,"WP_DEGs_Treeplot.png")
png(filename , width = 3000, height = 2000, res = 150)
plot(treeplot(res.wp.sim, label_format = 0.5, showCategory = 100, cluster.params = list(label_words_n = 0)))
dev.off()

filename <- paste0(out.folder,"WP_Upregulated_Treeplot.png")
png(filename , width = 3000, height = 2000, res = 150)
plot(treeplot(res.wp.up.sim, label_format = 0.5, showCategory = 100, cluster.params = list(label_words_n = 0)))
dev.off()

filename <- paste0(out.folder,"WP_Downregulated_Treeplot.png")
png(filename , width = 3000, height = 2000, res = 150)
plot(treeplot(res.wp.down.sim, label_format = 0.5, showCategory = 100, cluster.params = list(label_words_n = 0)))
dev.off()



