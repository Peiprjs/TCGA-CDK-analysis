if (!("rmarkdown" %in% installed.packages())) { install.packages("rmarkdown", update=FALSE) }
if (!("knitr" %in% installed.packages())) { install.packages("knit", update=FALSE) }

setwd("~/gits/BBS2061-breast-cancer/")
rmarkdown::render("scripts/Step0_Preamble.R")
rmarkdown::render("scripts/Step1_DataExploration.R")
rmarkdown::render("scripts/Step2_EnrichmentAnalysiss.R")
rmarkdown::render("scripts/Step3_PathwayVisualization.R")
rmarkdown::render("scripts/Step4_PPINetwork.R")
interest_cluster <- "3297" #3297 and 3283
rmarkdown::render("scripts/Step5_IdentifyClusters.R")
interest_cluster <- "3283" #3297 and 3283
rmarkdown::render("scripts/Step5_IdentifyClusters.R")
# rmarkdown::render("scripts/Step6-NumberDrugsTargetsTable.R")

rmarkdown::render("scripts/Step0_Preamble.R")
rmarkdown::render("scripts/Step1_DataExploration.R")
rmarkdown::render("scripts/Step2_EnrichmentAnalysiss.R")
rmarkdown::render("scripts/Step6_GO-Process_PPI_DrugExtension.R")
rmarkdown::render("scripts/Step7_Pathways2PPI_DrugExtension.R")

system("mv scripts/*.html ../BBS2061-breast-cancer-pages")