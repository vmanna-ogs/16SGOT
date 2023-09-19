---
title: "Companion scripts for 16Smetabarcoding in northeast MedSea paper Celussi, Manna et al. - 4"
author: "vmanna@OGS"
date: "2023-09-19"
output:
  html_document:
    keep_md: yes
---

### Relative abundance barplots
This script generates relative abundance barplots for each of the climactic taxon identified through indicator species analysis (Companion script 3) ad indicator of one of the three climax communities


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
library(readr)

#Phyloseq object handling
library(phyloseq)

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
library(ggnested)
```



```r
#Load pseq object
taxa_greater1<-readRDS("C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/taxa_greater1.rds")
indsp_res_clusters_filt<-read_csv("C:/Users/manna/OneDrive - OGS/SHAREMED_JOINT/Sharemed_code/output_data/indsp_res_clust_filt.csv")
```

```
## New names:
## Rows: 70 Columns: 7
## ── Column specification
## ──────────────────────────────────────────────────────── Delimiter: "," chr
## (2): Cluster, ASVs dbl (5): ...1, IndSP, Specificity, Fidelity, pval
## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
## Specify the column types or set `show_col_types = FALSE` to quiet this message.
## • `` -> `...1`
```

```r
#melt the pseq object with relative abundance of taxa greater than 1%
rel1_melt<-psmelt(taxa_greater1)
rel1_melt$Abundance<-rel1_melt$Abundance*100

#Subset by cluster keeping variable of interest

dark_df<-rel1_melt%>%
  filter(OTU%in%indsp_res_clusters_filt[which(indsp_res_clusters_filt$Cluster=="dark"),"ASVs"]$ASVs,
         cl_dark=="dark")%>%
  arrange(desc(Abundance),
          by_group=TRUE)%>%
  mutate(Genus = factor(Genus, 
                        levels = unique(Genus),
                        ordered = TRUE))

apr_df<-rel1_melt%>%
  filter(OTU%in%indsp_res_clusters_filt[which(indsp_res_clusters_filt$Cluster=="apr"),"ASVs"]$ASVs,
         cl_apr=="apr")%>%
  arrange(desc(Abundance),
          by_group=TRUE)%>%
  mutate(Genus = factor(Genus, 
                        levels = unique(Genus),
                        ordered = TRUE))

hot_df<-rel1_melt%>%
  filter(OTU%in%indsp_res_clusters_filt[which(indsp_res_clusters_filt$Cluster=="hot"),"ASVs"]$ASVs,
         cl_hot=="hot")%>%
  arrange(desc(Abundance),
          by_group=TRUE)%>%
  mutate(Genus = factor(Genus, 
                        levels = unique(Genus),
                        ordered = TRUE))




#bind together
clust_df<-rbind(dark_df,apr_df,hot_df)

#Arrange for descending abundance and transform Genus as ordered factor to keep the 
#ordering effective in plots
clust_df<-clust_df%>%
  arrange(desc(Abundance))%>%
  mutate(Genus = factor(Genus, levels = unique(Genus),ordered = TRUE))
```

The next step is to build the nested palette. We will use ggnested which will create a palette using on tints/shades of a base color, assigned to a given Genus, to color the ASVs within that Genus. This can be used for every nested factor as needed.


```r
#Take out the first most abundant 20 Genera. This is possible because we ordered this df according to descending abundance in the previous chunk
top20_genera<-levels(clust_df$Genus)[1:20]

#Group the other Genera as "Others" while keeping the others as they are
clust_df$Genus_top20<-ifelse(clust_df$Genus%in%top20_genera,
                             as.character(clust_df$Genus),"Other")

#Since we created a new facotr, we should make sure to refactor this to keep the order of appearance in plots
clust_df<-clust_df%>%
  mutate(Genus_top20 = factor(Genus_top20, 
                              levels = unique(Genus_top20),
                              ordered = TRUE))
#Final check
levels(clust_df$Genus_top20)
```

```
##  [1] "Ascidiaceihabitans"                 "Candidatus Nitrosopumilus"         
##  [3] "Synechococcus CC9902"               "unassigned_Marine Group II"        
##  [5] "NS4 marine group"                   "NS5 marine group"                  
##  [7] "Cohaesibacter"                      "DS001"                             
##  [9] "unassigned_Balneolaceae"            "SUP05 cluster"                     
## [11] "unassigned_AEGEAN-169 marine group" "Fluviicola"                        
## [13] "Cyanobium PCC-6307"                 "unassigned_NS11-12 marine group"   
## [15] "HIMB11"                             "OM60(NOR5) clade"                  
## [17] "Croceicoccus"                       "Vibrio"                            
## [19] "unassigned_NS9 marine group"        "unassigned_Clade III"              
## [21] "Other"
```

```r
#Create a factor storing the AVs belonging to the top20 genera. This will be the subgroups for shading
ASV_top20<-unique(clust_df%>%select(OTU,Genus_top20)%>%filter(!Genus_top20=="Other"))
clust_df$ASV_top20<-ifelse(clust_df$OTU%in%ASV_top20$OTU,
                          clust_df$OTU,"Other")
```

```r
#Manually create a plette with 20 distinct base colors for the 20 genera.
#ASVs within this genera will be assigned shades of the same color
pal_top_20_TT<-c("#eccc68", #Ascidiaceihabitans
                 "#B53471", #Nitrosopumilus
                 "#ff6348", #Synechococcus #ff6348
                 "#2ed573", #MGII
                 "#833471", #NS4
                 "#9980FA", #NS5
                 "#009432", #Stappiaceae
                 "#ff4757", #Coaesibacter
                 "#12CBC4", #DS001
                 "#1289A7", #Balneolaceae
                 "#0652DD", #SUP05
                 "#1B1464", #Aegean169
                 "#FDA7DF", #Fluviicola
                 "#D980FA", #Cyanobium
                 "#3742fa", #NS11 NS12
                 "#5758BB", #HIMB11
                 "#ffa502", #OM60
                 "#B53471", #Croceicoccus
                 "#ff7f50", #Vibrio
                 "#6F1E51", #NS9
                 "#a4b0be") #Other

#Name the palette with the Genus top 20 factor
#First rename genera with refined taxonomy to avoid all those "unassigned_x" etc.
levels(clust_df$Genus_top20)<-c("Ascidiaceihabitans", "Ca. Nitrosopumilus", "Synechococcus",
                                "Marine Group II", "NS4","NS5","Stappiaceae","Cohaesibacter",
                                "DS001","Balneolaceae","SUP05","Aegean-169","Fluviicola",
                                "Cyanobium","NS11-12","HIMB11","OM60","Croceicoccus","Vibrio",
                                "NS9","Other")
names(pal_top_20_TT)<-levels(clust_df$Genus_top20)
#Labels for facets
depth_label<-c("Surface","Bottom")
names(depth_label)<-c("surface","bottom")

HT_label<-"High temperature"
names(HT_label)<-"hot"

LL_label<-"Low light"
names(LL_label)<-"dark"

SP_label<-"Spring"
names(SP_label)<-"apr"

ALL_label<-c("High temperature","Low light","Spring")
names(ALL_label)<-c("hot","dark","apr")
```

The plots will be generated separately for each cluster plus the legend. We will pu together all the info later either through R or vector handling softwares


```r
#High temperature cluster
hot_plot_nested<-
  clust_df%>%filter(cluster=="hot")%>%
  ggnested(aes_string(main_group="Genus_top20",
                                sub_group="ASV_top20",
                                x="month_year",
                                y="Abundance"),
                     main_palette=pal_top_20_TT,
                     gradient_type = "tints")+
  geom_bar(stat="identity")+
 #ylab("Relative Abundance (%)\n")+
  ylim(0,60)+
  facet_grid(cols=vars(area,depth),
             rows = vars(cluster),
             scales = "free_x",
             space = "free",
             labeller = labeller(depth=depth_label,
                                 cluster=HT_label))+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
        axis.text = element_text(size=11),
        axis.title = element_text(size=15),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size=11),
        strip.text.y = element_text(size=13,face="bold"),
        panel.background = element_rect(fill=NA),
        axis.line = element_line(linewidth = 1.5,lineend = "round"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2,"mm"),
        axis.title.x = element_blank())
```

```
## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
## ℹ Please use tidy evaluation ideoms with `aes()`
```

```r
hot_plot_nested
```

![](4-cluster-nested-barplot_files/figure-html/Plotting/1-1.png)<!-- -->

```r
dark_plot_nested<-
  clust_df%>%filter(cluster=="dark")%>%
  ggnested::ggnested(aes_string(main_group="Genus_top20",
                                sub_group="ASV_top20",
                                x="month_year",
                                y="Abundance"),
                     main_palette=pal_top_20_TT,
                     gradient_type = "tints")+
  ylim(0,60)+
  xlab("Month_Year")+
  geom_bar(stat="identity")+
  #ylab("Relative Abundance (%)\n")+
  facet_grid(cols=vars(area,depth),
             rows=vars(cluster),
             scales = "free_x",
             space = "free",
             labeller = labeller(depth=depth_label,
                                 cluster=LL_label))+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
        axis.text = element_text(size=11),
        axis.title = element_text(size=15),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size=11),
        strip.text.y = element_text(size=13,face="bold"),
        panel.background = element_rect(fill=NA),
        axis.line = element_line(linewidth = 1.5,lineend = "round"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2,"mm"),
        axis.title.x = element_blank())
dark_plot_nested
```

![](4-cluster-nested-barplot_files/figure-html/Plotting/2-1.png)<!-- -->

```r
#Spring cluster
spring_plot_nested<-
  clust_df%>%filter(cluster=="apr")%>%
  ggnested::ggnested(aes_string(main_group="Genus_top20",
                                sub_group="ASV_top20",
                                x="month_year",
                                y="Abundance"),
                     main_palette=pal_top_20_TT,
                     gradient_type = "tints")+
  geom_bar(stat="identity")+
  #ylab("Relative Abundance (%)\n")+
  ylim(0,60)+
  facet_grid(cols=vars(area,depth),
             rows=vars(cluster),
             scales = "free_x",
             space = "free",
             labeller = labeller(depth=depth_label,
                                 cluster=SP_label))+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
        axis.text = element_text(size=11),
        axis.title = element_text(size=15),
        strip.text.x = element_text(size=11),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size=13,face="bold"),
        panel.background = element_rect(fill=NA),
        axis.line = element_line(linewidth = 1.5,lineend = "round"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2,"mm"),
        axis.title.x = element_blank())
spring_plot_nested
```

![](4-cluster-nested-barplot_files/figure-html/Plotting/3-1.png)<!-- -->

```r
legend_all<-
  get_legend(clust_df%>%
                         ggnested::ggnested(aes_string(main_group="Genus_top20",
                                                       sub_group="ASV_top20",
                                                       x="month_year",
                                                       y="Abundance"),
                                            main_palette=pal_top_20_TT,
                                            gradient_type = "tints",
                                            legend_title = "")+
                         geom_bar(stat="identity")+
                         ylab("Relative Abundance (%)\n")+
                         facet_grid(cols=vars(area,depth),rows=vars(cluster),scales = "free_x",space = "free")+
                         theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
                               panel.spacing.x=unit(-0.25,"lines"),
                               panel.border = element_rect(fill=NA,linewidth = 0.5),
                               strip.background = element_blank(),
                               strip.text = element_blank(),
                               legend.position = "bottom",
                               legend.key.size = unit(6,"mm"),
                               legend.key = element_blank(),
                               legend.background = element_blank(),
                               legend.box.background = element_blank(),
                               legend.direction = "horizontal"))
plot(legend_all)
```

![](4-cluster-nested-barplot_files/figure-html/Plotting/4 Legend-1.png)<!-- -->

```r
cluster_barplot<-ggarrange(hot_plot_nested,
          dark_plot_nested,
          spring_plot_nested,
          legend_all,
          nrow = 4,
          ncol = 1,
          align = "hv")
```

```
## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
## Placing graphs unaligned.
```

```
## Warning: Graphs cannot be horizontally aligned unless the axis parameter is
## set. Placing graphs unaligned.
```

```r
cluster_barplot<-annotate_figure(cluster_barplot,
                                 left = text_grob("Relative Abundance (%)", 
                                                  color = "black", 
                                                  rot = 90,
                                                  face = "bold",
                                                  size = 15,
                                                  hjust = -0.3))

plot(cluster_barplot)
```

![](4-cluster-nested-barplot_files/figure-html/Arrange plots together-1.png)<!-- -->




