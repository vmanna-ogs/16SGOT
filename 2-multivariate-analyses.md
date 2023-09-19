---
title: "Companion scripts for 16Smetabarcoding in northeast MedSea paper Celussi, Manna et al. - 2"
author: "vmanna@OGS"
date: "2023-09-19"
output:
  html_document:
    keep_md: yes
---

## Ordination analyses


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
library(Hmisc)
```

```
## 
## Attaching package: 'Hmisc'
## 
## The following objects are masked from 'package:dplyr':
## 
##     src, summarize
## 
## The following objects are masked from 'package:base':
## 
##     format.pval, units
```

```r
library(corrplot)
```

```
## corrplot 0.92 loaded
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
library(ggrepel)
```

```
## Warning: package 'ggrepel' was built under R version 4.2.3
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

```r
library(gridExtra)
```

```
## 
## Attaching package: 'gridExtra'
## 
## The following object is masked from 'package:dplyr':
## 
##     combine
```

### Further data cleaning

```r
#import pseq object of interest
sharemed_noint<-readRDS("C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/sharemed_noint.rds")
#There are a bunch of samples wich have NAs in the metadata. Create a new phyloseq object without those samples to proceed with the RDA
sharemed_noint_NA<-ps_drop_incomplete(sharemed_noint, 
                                      vars = c("temp","sal","doc","don",
                                               "tdn","ammon","nitrites",
                                               "nitrates","phosphate","chla"), 
                                      verbose = FALSE)
#save this as well
saveRDS(sharemed_noint_NA,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/sharemed_noint_NA.rds")
```

### RDA


```r
#Create a df storing env variables used in multivariate analisys
colnames(sample_data(sharemed_noint_NA))
```

```
##  [1] "area"        "date"        "julian_date" "depth"       "mm_text"    
##  [6] "mm_num"      "mm_yy"       "month_year"  "year"        "day_length" 
## [11] "prochl"      "prok_tot"    "hcp"         "syn"         "hb"         
## [16] "lna"         "hna"         "v1"          "v2"          "v3"         
## [21] "v4"          "hnf"         "pnf"         "peuk1"       "peuk2"      
## [26] "peuk3"       "crypto"      "temp"        "sal"         "density"    
## [31] "chla"        "poc"         "doc"         "ammon"       "nitrites"   
## [36] "nitrates"    "phosphate"   "don"         "tdn"         "dop"        
## [41] "tdp"         "tpn"
```

```r
env_data<-sample_data(sharemed_noint_NA)[,c(1,2,4,5,6,10,28,29,31,33,34,35,36,37,38,39)]
colnames(env_data)
```

```
##  [1] "area"       "date"       "depth"      "mm_text"    "mm_num"    
##  [6] "day_length" "temp"       "sal"        "chla"       "doc"       
## [11] "ammon"      "nitrites"   "nitrates"   "phosphate"  "don"       
## [16] "tdn"
```

```r
#check for collinearity
cormat<-cor(env_data[,6:16])
corrplot((cormat),"circle")
```

![](2-multivariate-analyses_files/figure-html/Collinearity checking-1.png)<!-- -->

```r
#Temperature vs. DOC and TDN vs DON are highly correlated.
#For the analysis we will keep only temperature and DON
#remove DOC and TDN
colnames(env_data)
```

```
##  [1] "area"       "date"       "depth"      "mm_text"    "mm_num"    
##  [6] "day_length" "temp"       "sal"        "chla"       "doc"       
## [11] "ammon"      "nitrites"   "nitrates"   "phosphate"  "don"       
## [16] "tdn"
```

```r
env_data<-env_data[,-c(10,16)]
#Scale environmental variables
env_data_scaled<-decostand(env_data[,6:14],method = "standardize")
#Identify scaled variable appending the suffix "scaled" in their colnames
newnames<-paste(colnames(env_data_scaled),"scaled",sep="_")
colnames(env_data_scaled)<-newnames
#check if rownames are in the same order prior binding
all.equal(rownames(env_data_scaled),rownames(sample_data(sharemed_noint_NA)))
```

```
## [1] TRUE
```

```r
ncol(sample_data(sharemed_noint_NA))
```

```
## [1] 42
```

```r
sample_data(sharemed_noint_NA)[,43:51]<-env_data_scaled
#cleanup
rm(cormat)
rm(newnames)
rm(env_data)
```
### Community data transformation
We will proceed with multivariate analyses using a subset of "abundant" taxa, meaning those with a relative abundance greater than 1% in at least one sample

```r
#Transform in relative abundance (0-1 range)
sharemed_noint_rel_NA<-microbiome::transform(sharemed_noint_NA,"compositional")
saveRDS(sharemed_noint_rel_NA,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/sharemed_noint_rel_NA.rds")
#Filtering out rare taxa. Keeping only those >1% in at least one sample
taxa_greater1_NA<-filter_taxa(sharemed_noint_rel_NA,function(x) max(x)>0.01,TRUE)
taxa_greater1_NA #260 taxa against 13755
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 260 taxa and 142 samples ]
## sample_data() Sample Data:       [ 142 samples by 51 sample variables ]
## tax_table()   Taxonomy Table:    [ 260 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 260 tips and 256 internal nodes ]
## refseq()      DNAStringSet:      [ 260 reference sequences ]
```

```r
saveRDS(taxa_greater1_NA,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/taxa_greater1_NA.rds")
#Give a look at how much of the community (in proportion) are we retaining with this filter
mean(rowSums(taxa_greater1_NA@otu_table@.Data))#78.41%
```

```
## [1] 0.7841881
```

```r
max(rowSums(taxa_greater1_NA@otu_table@.Data))#94.84%
```

```
## [1] 0.9468708
```

```r
min(rowSums(taxa_greater1_NA@otu_table@.Data))#61.17%
```

```
## [1] 0.6117122
```

```r
#on average, 12% of the community is composed of ASVs that have a relative abundance lower than 1% in at least one sample
```
### RDA on "abundant" subset

```r
#CAP ordination (pseq-based RDA) with 1% subset
cap_ord_scaled_1<-ordinate(physeq = taxa_greater1_NA,  
                           method = "CAP", 
                           distance = "bray", 
                           formula = ~ day_length_scaled+
                                       temp_scaled+
                                       sal_scaled+
                                       chla_scaled+
                                       ammon_scaled+
                                       nitrites_scaled+
                                       nitrates_scaled+
                                       phosphate_scaled+
                                       don_scaled)
#Testing the model
model_rda<-anova.cca(cap_ord_scaled_1,step=1000)
model_rda
```

```
## Permutation test for capscale under reduced model
## Permutation: free
## Number of permutations: 999
## 
## Model: capscale(formula = OTU ~ day_length_scaled + temp_scaled + sal_scaled + chla_scaled + ammon_scaled + nitrites_scaled + nitrates_scaled + phosphate_scaled + don_scaled, data = data, distance = distance)
##           Df SumOfSqs      F Pr(>F)    
## Model      9   18.702 14.965  0.001 ***
## Residual 132   18.329                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#Plot RDA trough phyloseq
cap_plot_scaled_rel_1<-
  plot_ordination(physeq = taxa_greater1_NA, 
                  ordination = cap_ord_scaled_1, 
                  color = "mm_text", 
                  shape = "depth",
                  axes = c(1,2)) + 
  geom_point(aes(colour = mm_text), size = 5) + 
  coord_fixed()+
  scale_color_discrete()
#Add environmental vectors as arrows
arrowmat_scaled_rel_1<- vegan::scores(cap_ord_scaled_1, display = "bp")
labels_vectors_rda<-c("Day Length","Temperature","Salinity",
                       "Chl a","Ammonia","Nitrites",
                       "Nitrates","Phosphate","DON")
#Create a df storing aestethic information for arrows
arrowdf_scaled_rel_1<- data.frame(labels_scaled_rel_1=labels_vectors_rda,
                                  arrowmat_scaled_rel_1)
#Define the arrow aesthetic mapping 
arrow_map_scaled_rel_1<-aes(xend = CAP1, 
                            yend = CAP2, 
                            x = 0, y = 0, 
                            shape = NULL, 
                            color = NULL, 
                            label = labels_scaled_rel_1)
label_map_scaled_rel_1<- aes(x = 1.3 * CAP1, 
                             y = 1.3 * CAP2, 
                             shape = NULL, 
                             color = NULL,
                             label = labels_scaled_rel_1)
arrowhead_scaled_rel_1= arrow(length = unit(0.02, "npc")) 
#Final plot
RDA_plot<-cap_plot_scaled_rel_1 + 
  scale_color_viridis_d(direction = -1)+
  scale_shape_discrete(breaks=c("surface","bottom"),
                       labels=c("Surface","Bottom"))+
  geom_segment(mapping = arrow_map_scaled_rel_1,
               size = 1.5, 
               data = arrowdf_scaled_rel_1,
               color = "black", 
               arrow = arrowhead_scaled_rel_1)+ 
  geom_text(mapping = label_map_scaled_rel_1,
            size = 5,  
            data = arrowdf_scaled_rel_1, 
            show.legend = F) + 
  theme(panel.background=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.5,
                                    color = "black",
                                    fill=NA),
        axis.text = element_text(size=12,
                                 colour = "black"),
        axis.ticks = element_line(linewidth = 0.8,
                                  color = "black"),
        axis.title = element_text(size=14, face="plain",
                                  colour = "black"),
        axis.title.x = element_text(vjust=-1.5),
        axis.title.y = element_text(vjust=+1.5),
        legend.text = element_text(size=12,
                                   colour = "black"),
        legend.title = element_text(size=12,
                                    colour = "black"),
        legend.key = element_blank(),
        legend.box = "vertical",
        legend.position = "bottom")+
  labs(color="Month",
       shape="Depth")+
  coord_fixed(ratio=1)
```

```
## Scale for colour is already present.
## Adding another scale for colour, which will replace the existing scale.
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
```

```
## Warning in geom_segment(mapping = arrow_map_scaled_rel_1, size = 1.5, data =
## arrowdf_scaled_rel_1, : Ignoring unknown aesthetics: label
```

```
## Coordinate system already present. Adding new coordinate system, which will
## replace the existing one.
```

```r
RDA_plot
```

![](2-multivariate-analyses_files/figure-html/RDA output -1.png)<!-- -->


### PCoA
Since for this we do not need environmental data, we will use all the samples. Naturally, there will be the subset at 1% relative abundance in at least one sample.

```r
#Transform in relative abundance (0-1 range)
sharemed_noint_rel<-microbiome::transform(sharemed_noint,"compositional")
#save this object for later
saveRDS(sharemed_noint_rel,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/sharemed_noint_rel.rds")
##Filtering out rare taxa
#1% in at least one sample
taxa_greater1<-filter_taxa(sharemed_noint_rel,function(x) max(x)>0.01,TRUE)
#BC dissimilarities
BC_rel1<-phyloseq::distance(taxa_greater1, method = "bray")

#Compute PCoA
pcoa_rel1<-ordinate(taxa_greater1,
                    method = "PCoA",
                    distance = BC_rel1)
```

```r
plot_ordination(taxa_greater1,pcoa_rel1,
                type="samples",
                shape="depth",
                color="mm_text")+
  scale_shape_discrete(breaks=c("surface","bottom"),
                       labels=c("Surface","Bottom"))+
  geom_point(size=5)+
  theme(panel.background=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.5,
                                    color = "black",
                                    fill=NA),
        axis.text = element_text(size=12,
                                 colour = "black"),
        axis.ticks = element_line(linewidth = 0.8,
                                  color = "black"),
        axis.title = element_text(size=14, face="plain",
                                  colour = "black"),
        axis.title.x = element_text(vjust=-1.5),
        axis.title.y = element_text(vjust=+1.5),
        legend.text = element_text(size=12,
                                   colour = "black"),
        legend.title = element_text(size=12,
                                    colour = "black"),
        legend.key = element_blank(),
        legend.box = "vertical",
        legend.position = "bottom")+
  labs(color="Month",
       shape="Depth")+
  coord_fixed(ratio = 1)
```

![](2-multivariate-analyses_files/figure-html/PCoA output-1.png)<!-- -->
On this plot we will identify the most distant points to define the three clusters of climax samples.
This was done outside R. Clusters were defined starting from the most extreme points (i.e., triangle vertices) and then taking the neighbouring samples falling into the 10th percentile of the Bray-Curtis dissimilarity indices.
The result is:

```r
rel1_apr1<-c("C1_Apr19_bott","C1_Apr19_surf","C1_Apr20_bott","C1_Apr20_surf","C1_Mar19_surf",
             "C1_Apr21_surf","C1_May20_bott","C1_May21_bott","BF_May19_surf","BF_Apr20_bott",
             "BF_Apr19_surf","BF_Apr19_bott","BF_May20_bott","BF_Apr21_surf","BF_Apr21_bott",
             "C1_Apr21_bott")
rel1_hot1<-c("C1_Jul19_surf","C1_Ago20_surf", "C1_Jul21_surf","C1_Sep20_surf","C1_Aug21_bott",
             "C1_Aug21_surf","C1_Sep21_bott","C1_Sep21_surf","BF_Oct18_surf","BF_Jul19_surf",
             "BF_Aug20_surf","BF_Aug19_surf","BF_Oct19_surf","BF_Jul20_surf","BF_Aug21_surf",
             "BF_Sep21_surf")
rel1_dark1<-c("C1_Dec18_bott","C1_Dec18_surf","C1_Dec19_bott","C1_Dec20_bott","C1_Dec20_surf",
              "C1_Jan21_bott","C1_Jan21_surf","C1_Dec21_bott","C1_Dec21_surf","BF_Dec21_bott",
              "BF_Dec20_surf","BF_Dec19_bott","BF_Dec21_surf","BF_Dec18_bott","BF_Dec20_bott")
#Create a factor in the pseq object metadata storing the information about clusters
sample_data(taxa_greater1)$cluster<-ifelse(sample_names(taxa_greater1)%in%rel1_apr1,"apr",
                                           ifelse(sample_names(taxa_greater1)%in%rel1_hot1,"hot",
                                                  ifelse(sample_names(taxa_greater1)%in%rel1_dark1,"dark",NA)))
#Transform as factor
sample_data(taxa_greater1)$cluster<-factor(sample_data(taxa_greater1)$cluster,
                                           levels = c("apr","hot","dark"),
                                           ordered = TRUE)

#Colors for cluster factor
grouping_cols<-c("hot"="#E74C3C",
                 "dark"="#3498DB",
                 "apr"="#1ABC9C")
```

```r
#Plot PCoA highlighting climax communities clusters
PCoA_plot<-plot_ordination(taxa_greater1,
                pcoa_rel1,
                type="samples",
                color="cluster")+
  geom_point(size=6,aes(shape=depth))+
  scale_color_manual(values=grouping_cols,
                     breaks=c("apr","hot","dark"),
                     labels=c("Spring",
                              "High temperature",
                              "Low light"))+
  scale_shape_discrete(breaks=c("surface","bottom"),
                       labels=c("Surface","Bottom"))+
  geom_text_repel(label =taxa_greater1@sam_data$mm_num, 
                  size = 4, 
                  nudge_y =-0.025, 
                  color = "black",
                  max.overlaps = 20)+
  labs(color="Cluster",shape="Depth")+
  theme(panel.background=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.5,
                                    color = "black",
                                    fill=NA),
        axis.text = element_text(size=12,
                                 colour = "black"),
        axis.ticks = element_line(linewidth = 0.8,
                                  color = "black"),
        axis.title = element_text(size=14, face="plain",
                                  colour = "black"),
        axis.title.x = element_text(vjust=-1.5),
        axis.title.y = element_text(vjust=+1.5),
        legend.text = element_text(size=12,
                                   colour = "black"),
        legend.title = element_text(size=12,
                                    colour = "black"),
        legend.key = element_blank(),
        legend.box = "vertical",
        legend.position = "bottom")+
  coord_fixed(ratio=1)
PCoA_plot
```

![](2-multivariate-analyses_files/figure-html/PCoA output clusters -1.png)<!-- -->


```r
#Join RDA and PCoA on the same plot
RDA_PCoA_composite<-ggarrange(RDA_plot,PCoA_plot,ncol =2,common.legend = FALSE,legend = "bottom",align = "hv")
RDA_PCoA_composite
```

![](2-multivariate-analyses_files/figure-html/RDA+PCoA plotting-1.png)<!-- -->

### Do some cleaning 

```r
saveRDS(taxa_greater1,"C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/taxa_greater1.rds")
rm(labels_vectors_rda)
rm(BC_rel1)
rm(taxa_greater1_NA)
rm(sharemed_noint_NA)
rm(sharemed_noint_rel_NA)
rm(RDA_plot)
rm(PCoA_plot)
rm(pcoa_rel1)
rm(RDA_PCoA_composite)
rm(cap_ord_scaled_1)
rm(cap_plot_scaled_rel_1)
rm(model_rda)
rm(label_map_scaled_rel_1)
rm(arrow_map_scaled_rel_1)
rm(arrowmat_scaled_rel_1)
rm(arrowdf_scaled_rel_1)
rm(arrowhead_scaled_rel_1)
rm(env_data_scaled)
```
