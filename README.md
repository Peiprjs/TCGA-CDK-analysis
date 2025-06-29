# BBS2061 - 2A - Topic A (Network Biology) 

[![DOI](https://zenodo.org/badge/1006748923.svg)](https://doi.org/10.5281/zenodo.15767005)\

For the course BBS2061, we were instructed to modify/create a code and run a set of computational analyses, in order to answer a set research question. 
The project involved groups of 3-4 people, which were given a certain topic regarding the course information/content. 
This respective project involves 4 individuals; Roca Cugat M.(*i6351071*) ; Sapsai I.(*i6355626*) ; Rijk IR.(*i6342444*) ; Mirensky Roffe DS.(*i6323580*).
The focus/topic of the group is the Topic A (Network Biology)

## Instructions and data. 

### Research Question
> How is the G1 to S Cell Cycle Altered in Stage IA Breast Cancer Cells, Specifically Cyclins and its Inhibitors Involved

The publicly available RNA sequencing data from the TCGA breast cancer project was used for this analysis [1]. Twenty patients were selected based on age, sex, ethnicity, and tumor stage. All patients were female, Caucasian, between 60 and 75 years old, and diagnosed with stage IA breast cancer. Information on hormonal receptor expression was not available for these patients. Gene expression data for both the primary tumor and matched normal tissue were extracted. The RNA-sequencing count data was pre-processed to remove lowly expressed genes. A paired differential gene expression analysis between tumor and matched normal tissue was performed using the DESeq2 package in R [2].

After careful analysis of the data, Wikipathways was utilized to differentiate and to find a specific pathway that related highly to the provided data. 
The chosen pathways (provided/explained in the written report) focused on the cell cyle, CDK and cyclin inhibitors which mainly focused on modulating CDK4/6 responses. 
The cell cycle pathways were chosen due to their important role in the tumor growth/proliferation, which peaked the interest of the respective group memnbers. 
Previously aquired knowledge led to the foundation of a possible hypotehsis, which can be found in the line 19. 

### Hypothesis
> The Cylcin-Dependent Kinase 4/6 Inhibtors will be Downregulated. 

You can find the execution of all the code steps in [https://peiprjs.github.io/BBS2061-breast-cancer/](https://peiprjs.github.io/BBS2061-breast-cancer/)

## References
[1] Cancer Genome Atlas Network (2012) “Comprehensive molecular portraits of human breast tumours.” Nature, 490(7418):61-70. doi: 10.1038/nature11412
[2] Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140

# Steps

[Step 0](Step0_Preamble.html)

[Step 1](Step1_DataExploration.html)

[Step 2](Step2_EnrichmentAnalysiss.html)

[Step 3](Step3_PathwayVisualization.html)

[Step 4](Step4_PPINetwork.html)

[Step 5](Step5_IdentifyClusters.html)

[Step 6](Step6-NumberDrugsTargetsTable.html)

[Step 7](Step7_GO-Process_PPI_DrugExtension.html)

[Step 8](Step8_Pathways2PPI_DrugExtension.html)

[Tools used](Toolsused.html)
