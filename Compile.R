if (!("rmarkdown" %in% installed.packages())) { install.packages("rmarkdown", update=FALSE) }
if (!("knitr" %in% installed.packages())) { install.packages("knit", update=FALSE) }

rmarkdown::render("scripts/Step0_Preamble.R")
rmarkdown::render("scripts/Step1_DataExploration.R")
rmarkdown::render("scripts/Step2_EnrichmentAnalysiss.R")
rmarkdown::render("scripts/Step3_PathwayVisualization.R")
rmarkdown::render("scripts/Step4_PPINetwork.R")
rmarkdown::render("scripts/Step5_IdentifyClusters.R")
rmarkdown::render("scripts/Step6-NumberDrugsTargetsTable.R")
rmarkdown::render("scripts/Step7_GO-Process_PPI_DrugExtension.R")
rmarkdown::render("scripts/Step8_Pathways2PPI_DrugExtension.R")

system("mv scripts/*.html ../BBS2061-breast-cancer-pages")
