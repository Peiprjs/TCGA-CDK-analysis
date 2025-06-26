# Script:       Step1_DataExploration.R
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
# DIFFERENTIALLY EXPRESSED GENES
# #############################################

# Let's see how many genes are differentially expressed in our dataset
degs <- data[abs(data$log2FC) > log2fc.cutoff & data$adj.P.Value < pvalue.cutoff,]

# let's write out the table with all differentially expressed genes
write.table(degs, file=paste0(out.folder,"degs.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

log2fc.cutoff <- -1
genes.down <- data[data$log2FC < log2fc.cutoff & data$adj.P.Value < pvalue.cutoff,]
write.table(degs, file=paste0(out.folder,"degs_down.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

log2fc.cutoff <- 1
genes.up <- data[data$log2FC > log2fc.cutoff & data$adj.P.Value < pvalue.cutoff,]
write.table(degs, file=paste0(out.folder,"degs_up.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

# #############################################
# VOLCANO PLOT
# #############################################

# Let's create a Volcano plot to get a better understand of the intensity and
# direction of the changed genes

EnhancedVolcano(data, title = paste0("Breast cancer vs. Healthy (",nrow(degs), " DEGs)"), lab = data$GeneName, x = "log2FC", y = "adj.P.Value", pCutoff = pvalue.cutoff, FCcutoff = log2fc.cutoff, labSize = 3, xlim = c(-15,15), ylim=c(0,4.5))

# the code below saves the figure in a file in our output folder
filename <- paste0(out.folder,"volcano-plot.png")
png(filename , width = 2000, height = 1500, res = 150)
EnhancedVolcano(data, title = paste0("Breast cancer vs. Healthy (",nrow(degs), " DEGs)"), lab = data$GeneName, x = "log2FC", y = "adj.P.Value", pCutoff = pvalue.cutoff, FCcutoff = log2fc.cutoff, labSize = 3, xlim = c(-15,15), ylim=c(0,4.5))
dev.off()
