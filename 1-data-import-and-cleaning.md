---
title: "Companion scripts for 16Smetabarcoding in northeast MedSea paper Celussi, Manna et al. - 1"
author: "vmanna@OGS"
date: "2023-09-19"
output:
  html_document:
    keep_md: yes
---

## Data import and cleaning

### Loading packages

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
```

### Import files

```r
Sharemed<-readRDS("C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/input_data/Sharemed.rds")
Sharemed_envdata_20032023<-read.csv("C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/input_data/sharemed_envdata_20032023.csv")
metadata<-Sharemed_envdata_20032023
```

### Replace metadata inside phyloseq object

```r
#Assign sample names as rownames
rownames(metadata)<-metadata$sample_name
#Define some factors useful for later analyses/visualization
metadata$area<-factor(metadata$area,levels = c("C1","BF"),ordered = TRUE)
metadata$depth<-factor(metadata$depth,levels = c("surface","bottom"),ordered = TRUE)
metadata$mm_text<-factor(metadata$mm_text,levels = month.abb,ordered = TRUE)
metadata$mm_num<-factor(metadata$mm_num,levels = as.character(seq(1:12)),ordered = TRUE)
metadata$mm_yy<-factor(metadata$mm_yy,levels = c("10_2018","11_2018","12_2018","01_2019","02_2019","03_2019",
                                                 "04_2019","05_2019","06_2019","07_2019","08_2019","09_2019",
                                                 "10_2019","11_2019","12_2019","01_2020","02_2020","03_2020",
                                                 "04_2020","05_2020","06_2020","07_2020","08_2020","09_2020",
                                                 "10_2020","11_2020","12_2020","01_2021","02_2021","03_2021",
                                                 "04_2021","05_2021","06_2021","07_2021","08_2021","09_2021",
                                                 "10_2021","11_2021","12_2021"),ordered = TRUE)
metadata$year<-factor(metadata$year,levels = c("2018","2019","2020","2021"),ordered = TRUE)
#Put back metadata in the pseq object
sample_data(Sharemed)<-metadata[,-c(1)]
```

### Remove samples other than surface and bottom

```r
#Define a pattern to identify samples to remove. Including just "5m" since the pattern would be detected also in strings containing "15m"
pattern<-c("5m","10m")
samples_to_remove<-map(pattern,str_subset,string=sample_names(Sharemed))
samples_to_remove<-unlist(samples_to_remove)
#Subset the phyloseq object to remove the unwanted samples
Sharemed #186 samples
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 13790 taxa and 186 samples ]
## sample_data() Sample Data:       [ 186 samples by 42 sample variables ]
## tax_table()   Taxonomy Table:    [ 13790 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 13790 tips and 13751 internal nodes ]
## refseq()      DNAStringSet:      [ 13790 reference sequences ]
```

```r
#Remove the samples
sharemed_noint<-subset_samples(Sharemed, !sample_names(Sharemed) %in% samples_to_remove)
sharemed_noint #150 samples
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 13790 taxa and 150 samples ]
## sample_data() Sample Data:       [ 150 samples by 42 sample variables ]
## tax_table()   Taxonomy Table:    [ 13790 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 13790 tips and 13751 internal nodes ]
## refseq()      DNAStringSet:      [ 13790 reference sequences ]
```
### Revise taxonomy table

```r
#If a cell contains NA, this loop replace the NA with the prefix "unassigned_" followed by the previously known taxonomic level so that we have all formatted in the right way and avoid problems in merging by taxonomic levels

tax_tab_torevise<-tax_table(sharemed_noint)

for (col in 2:ncol(tax_tab_torevise)) {
  for (row in 1:nrow(tax_tab_torevise)) {
    if (is.na(tax_tab_torevise[row,col])) {
      if (!grepl("feature_id", tax_tab_torevise[row,col-1]) & !grepl("unassigned", tax_tab_torevise[row,col-1])) {
        tax_tab_torevise[row,col] <- paste0("unassigned_", tax_tab_torevise[row,col-1])
      } else {
        tax_tab_torevise[row,col] <- tax_tab_torevise[row,col-1]
      }
    }
  }
}
#Replace tax table in the pseq object
tax_table(sharemed_noint)<-tax_tab_torevise
#by removing some samples, there are probably a number of taxa which sum up to zero because they were only present in those samples removed. Check for this:
any(taxa_sums(sharemed_noint)==0)
```

```
## [1] TRUE
```

```r
#ok, this is TRUE. Remove those taxa
sharemed_noint<-prune_taxa(taxa_sums(sharemed_noint)>0,sharemed_noint) 
#before 13790, after 13755
#Save the resulting phyloseq object
saveRDS(sharemed_noint,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/sharemed_noint.rds")
```
### Cleanup 

```r
rm(pattern)
rm(col)
rm(row)
rm(samples_to_remove)
rm(tax_tab_torevise)
rm(Sharemed)
rm(Sharemed_envdata_20032023)
```
