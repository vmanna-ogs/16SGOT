# 16SGOT
Repo for the article:Annual recurrence of prokaryotic climax communities in shallow waters of the North Mediterranean
https://doi.org/10.1111/1462-2920.16595

Scripts are numbered sequentially:
1-data-import-and-cleaning performs the upload of environmental and sequence data (as phyloseq object), located in the input_data folder.
2-multivariate-analyses gathers the scripts used for PCoA and NMDS analysis and the definition of the recurrent clusters.
3-indicator-species-analysis performs the indicspecies analysis based on the recurrent clusters identified by multivariate analyses. 
4-cluster-nested-barplot is the script needed to generate the supplementary figures depicting nested barplot for the relative abundance of climactic taxa.
5-seasonality collects the scripts needed to perform a seasonality check on individual ASVs as well as their response to environmental variables.

The folder outpt_files collects intermediate files generated and needed to perform different analyses if you don't want to launch all the scripts.
