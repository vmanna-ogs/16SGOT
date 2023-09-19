---
title: "Companion scripts for 16Smetabarcoding in northeast MedSea paper Celussi, Manna et al. - 3"
author: "vmanna@OGS"
date: "2023-09-18"
output:
  html_document:
    keep_md: yes
---
## Indicator species analysis
### Packages needed

```r
#Data handling
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.0.10     ✔ readr     2.1.4 
## ✔ forcats   1.0.0      ✔ stringr   1.5.0 
## ✔ ggplot2   3.4.1      ✔ tibble    3.1.8 
## ✔ lubridate 1.9.0      ✔ tidyr     1.3.0 
## ✔ purrr     1.0.1      
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the ]8;;http://conflicted.r-lib.org/conflicted package]8;; to force all conflicts to become errors
```

```r
#Phyloseq object handling
library(phyloseq)
library(microbiome)
```

```
## 
## microbiome R package (microbiome.github.com)
##     
## 
## 
##  Copyright (C) 2011-2022 Leo Lahti, 
##     Sudarshan Shetty et al. <microbiome.github.io>
## 
## 
## Attaching package: 'microbiome'
## 
## The following object is masked from 'package:ggplot2':
## 
##     alpha
## 
## The following object is masked from 'package:base':
## 
##     transform
```

```r
library(microViz)
```

```
## Warning: package 'microViz' was built under R version 4.2.3
```

```
## microViz version 0.10.10 - Copyright (C) 2023 David Barnett
## ! Website: https://david-barnett.github.io/microViz
## ✔ Useful?  For citation details, run: `citation("microViz")`
## ✖ Silence? `suppressPackageStartupMessages(library(microViz))`
```

```r
#Statistics
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.6-4
## 
## Attaching package: 'vegan'
## 
## The following object is masked from 'package:microbiome':
## 
##     diversity
```

```r
library(indicspecies)
```

```
## Warning: package 'indicspecies' was built under R version 4.2.3
```

```
## 
## Attaching package: 'indicspecies'
## 
## The following object is masked from 'package:microbiome':
## 
##     coverage
```

```r
#Plotting
library(Cairo)
```

```
## Warning: package 'Cairo' was built under R version 4.2.3
```

```r
library(viridis)
```

```
## Loading required package: viridisLite
```

```r
library(ggpubr)
```

```
## 
## Attaching package: 'ggpubr'
## 
## The following object is masked from 'package:microViz':
## 
##     stat_chull
```
Here the goal is to identify ASVs associated with each group of samples. We will do this by checking each cluster against all the other samples.

```r
#import data
taxa_greater1<-readRDS("C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/taxa_greater1.rds")
#Replace all the NAs (i.e., samples that have not been assigned to any of the three clusters) with "all"
sample_data(taxa_greater1)$cluster<-ifelse(is.na(sample_data(taxa_greater1)$cluster),
                                           "all",
                                           as.character(sample_data(taxa_greater1)$cluster))
#create new factors to check each level (i.e., apr, hot, and dark) against all the other samples

sample_data(taxa_greater1)$cl_apr<-ifelse(sample_data(taxa_greater1)$cluster=="apr",
                                          "apr",
                                          "all")
sample_data(taxa_greater1)$cl_hot<-ifelse(sample_data(taxa_greater1)$cluster=="hot",
                                          "hot",
                                          "all")
sample_data(taxa_greater1)$cl_dark<-ifelse(sample_data(taxa_greater1)$cluster=="dark",
                                           "dark",
                                           "all")

#prepare for indicator species analysis
#community matrix
com_mat_rel1_vs<-as.matrix(otu_table(taxa_greater1))
#vectors of partitioning
#Spring cluster - apr
clustvec_apr<-sample_data(taxa_greater1)[,"cl_apr"]
names_apr<-rownames(clustvec_apr)
clustvec_apr<-factor(clustvec_apr$cl_apr)
names(clustvec_apr)<-names_apr

#Summer cluster - hot
clustvec_hot<-sample_data(taxa_greater1)[,"cl_hot"]
names_hot<-rownames(clustvec_hot)
clustvec_hot<-factor(clustvec_hot$cl_hot)
names(clustvec_hot)<-names_hot

#Winter cluster - dark
clustvec_dark<-sample_data(taxa_greater1)[,"cl_dark"]
names_dark<-rownames(clustvec_dark)
clustvec_dark<-factor(clustvec_dark$cl_dark)
names(clustvec_dark)<-names_dark
```


```r
indics_aprvsall<-multipatt(com_mat_rel1_vs,clustvec_apr,func = "IndVal.g",
                           duleg=TRUE,control = how(nperm=9999))
indics_hotvsall<-multipatt(com_mat_rel1_vs,clustvec_hot,func = "IndVal.g",
                           duleg=TRUE,control = how(nperm=9999))
indics_darkvsall<-multipatt(com_mat_rel1_vs,clustvec_dark,func = "IndVal.g",
                            duleg=TRUE,control = how(nperm=9999))
```


```r
#Spring cluster vs all samples
indsp_net_apr <- c()
for (i in unique(clustvec_apr)) {
  i_stats <- indics_aprvsall$sign[, c(which(names(indics_aprvsall$sign) == paste0("s.", i)),
                                      (length(unique(clustvec_apr))+1):ncol(indics_aprvsall$sign))]
  i_stats$specificity <- sapply(1:nrow(i_stats), function(k) indics_aprvsall$A[k, i_stats[k,2]]) # available only with IndVal.g
  i_stats$fidelity <- sapply(1:nrow(i_stats), function(k) indics_aprvsall$B[k, i_stats[k,2]]) # available only with IndVal.g
  i_stats_clean <- i_stats[intersect(which(i_stats[,1] == 1), which(i_stats$p.value < 0.05)), ] # it would creates separated islands in the network
  i_stats_clean2 <- i_stats_clean[which(i_stats_clean$p.value < 0.05), ]
  
  if (is.null(unlist(i_stats_clean2$specificity))) { # r.g
    i_df <- data.frame(Cluster = rep(i, nrow(i_stats_clean2)),
                       Traits = row.names(i_stats_clean2),
                       IndSP = i_stats_clean2$stat,
                       pval = i_stats_clean2$p.value)
  } else { # IndVal.g
    i_df <- data.frame(Cluster = rep(i, nrow(i_stats_clean2)),
                       ASVs = row.names(i_stats_clean2),
                       IndSP = i_stats_clean2$stat,
                       Specificity = i_stats_clean2$specificity,
                       Fidelity = i_stats_clean2$fidelity,
                       pval = i_stats_clean2$p.value)
  }
  
  indsp_net_apr <- rbind(indsp_net_apr, i_df)
}

#Summer cluster vs all samples
indsp_net_hot <- c()
for (i in unique(clustvec_hot)) {
  i_stats <- indics_hotvsall$sign[, c(which(names(indics_hotvsall$sign) == paste0("s.", i)),
                                      (length(unique(clustvec_hot))+1):ncol(indics_hotvsall$sign))]
  i_stats$specificity <- sapply(1:nrow(i_stats), function(k) indics_hotvsall$A[k, i_stats[k,2]]) # available only with IndVal.g
  i_stats$fidelity <- sapply(1:nrow(i_stats), function(k) indics_hotvsall$B[k, i_stats[k,2]]) # available only with IndVal.g
  i_stats_clean <- i_stats[intersect(which(i_stats[,1] == 1), which(i_stats$p.value < 0.05)), ] # it would creates separated islands in the network
  i_stats_clean2 <- i_stats_clean[which(i_stats_clean$p.value < 0.05), ]
  
  if (is.null(unlist(i_stats_clean2$specificity))) { # r.g
    i_df <- data.frame(Cluster = rep(i, nrow(i_stats_clean2)),
                       Traits = row.names(i_stats_clean2),
                       IndSP = i_stats_clean2$stat,
                       pval = i_stats_clean2$p.value)
  } else { # IndVal.g
    i_df <- data.frame(Cluster = rep(i, nrow(i_stats_clean2)),
                       ASVs = row.names(i_stats_clean2),
                       IndSP = i_stats_clean2$stat,
                       Specificity = i_stats_clean2$specificity,
                       Fidelity = i_stats_clean2$fidelity,
                       pval = i_stats_clean2$p.value)
  }
  
  indsp_net_hot <- rbind(indsp_net_hot, i_df)
}

#Winter cluster vs all samples
indsp_net_dark <- c()
for (i in unique(clustvec_dark)) {
  i_stats <- indics_darkvsall$sign[, c(which(names(indics_darkvsall$sign) == paste0("s.", i)),
                                       (length(unique(clustvec_dark))+1):ncol(indics_darkvsall$sign))]
  i_stats$specificity <- sapply(1:nrow(i_stats), function(k) indics_darkvsall$A[k, i_stats[k,2]]) # available only with IndVal.g
  i_stats$fidelity <- sapply(1:nrow(i_stats), function(k) indics_darkvsall$B[k, i_stats[k,2]]) # available only with IndVal.g
  i_stats_clean <- i_stats[intersect(which(i_stats[,1] == 1), which(i_stats$p.value < 0.05)), ] # it would creates separated islands in the network
  i_stats_clean2 <- i_stats_clean[which(i_stats_clean$p.value < 0.05), ]
  
  if (is.null(unlist(i_stats_clean2$specificity))) { # r.g
    i_df <- data.frame(Cluster = rep(i, nrow(i_stats_clean2)),
                       Traits = row.names(i_stats_clean2),
                       IndSP = i_stats_clean2$stat,
                       pval = i_stats_clean2$p.value)
  } else { # IndVal.g
    i_df <- data.frame(Cluster = rep(i, nrow(i_stats_clean2)),
                       ASVs = row.names(i_stats_clean2),
                       IndSP = i_stats_clean2$stat,
                       Specificity = i_stats_clean2$specificity,
                       Fidelity = i_stats_clean2$fidelity,
                       pval = i_stats_clean2$p.value)
  }
  
  indsp_net_dark <- rbind(indsp_net_dark, i_df)
}

#Bind all together in one dataframe
indsp_res_clusters<-rbind(indsp_net_apr[which(indsp_net_apr$Cluster=="apr"),],
                          indsp_net_hot[which(indsp_net_hot$Cluster=="hot"),],
                          indsp_net_dark[which(indsp_net_dark$Cluster=="dark"),])

#Filter ASVs based on fidelity and specificity >0.8
indsp_res_clusters_filt<-indsp_res_clusters%>%filter(Fidelity>=0.8, Specificity>=0.8)

#Save this data frame
write.csv(indsp_res_clusters_filt,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/indsp_res_clust_filt.csv")
```


```r
saveRDS(taxa_greater1,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/taxa_greater1.rds")
rm(com_mat_rel1_vs)
rm(clustvec_apr)
rm(clustvec_dark)
rm(clustvec_hot)
rm(i)
rm(i_df)
rm(i_stats)
rm(i_stats_clean)
rm(i_stats_clean2)
rm(indics_aprvsall)
rm(indics_darkvsall)
rm(indics_hotvsall)
rm(indsp_res_clusters)
```
