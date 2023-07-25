library(phyloseq)
library(microbiome)
library(microViz)
library(tidyverse)
library(propr)
library(lomb)
library(lubridate)
library(fdrtool)

#lomb-scargle periodogram
sharemed_noint
sharemed_transf<-data_transformation(sharemed_noint)
sharemed_clr<-sharemed_transf$clr


#melt the clr-transformed dataset and transform date as decimal
tsdf <- sharemed_clr %>% 
  psmelt() %>% 
  mutate(decimaldat = decimal_date(date)) 

#Arrange by date and subset, creating a list. Each ASV is an element of the list;
#inside there are abundance data for each sampling time
df.ts <- tsdf %>% 
  arrange(date) %>% 
  select(date, decimaldat, Abundance, OTU) %>% 
  split(.$OTU, drop = T)

#Run the lomb-scargle periodogram
lomb_res<-df.ts %>% 
  map(~lsp( x =.x$Abundance,
            times = .x$decimaldat,
            type = 'period',
            normalize="press",
            plot = F))

#put the results in a table, subsetting for only those ASVs with PN>10 and q<0.05, calculated through FDR correction
lomb_tab<-tibble( asv = names(lomb_res),
                  pval = map_dbl(lomb_res, ~.x$p.value),
                  peak = map_dbl(lomb_res, ~.x$peak),
                  interval  = map(lomb_res, ~.x$peak.at),
                  int.min = map_dbl(interval, ~.[[2]]),
                  int.max = map_dbl(interval, ~.[[1]])) %>% 
  mutate( qval = fdrtool::fdrtool(pval, statistic = 'pvalue')$qval) %>% 
  filter(qval <= 0.01, peak >= 10, int.max <= 2)

results_lomb <- lomb_res[lomb_tab %>% pull(asv)]


#now melt the relative abundance data
sharemed_rel_melt<-psmelt(sharemed_noint_rel)
colnames(sharemed_rel_melt)
#create a factor storing the information about the periodicity of the ASVs
sharemed_rel_melt$periodic<-ifelse(sharemed_rel_melt$OTU%in%lomb_tab$asv,
                                   "periodic",
                                   "non_periodic")
#The same thing for the clr dataframe
sharemed_clr_melt<-psmelt(sharemed_clr)
sharemed_clr_melt$periodic<-ifelse(sharemed_clr_melt$OTU%in%lomb_tab$asv,
                                   "periodic",
                                   "non_periodic")
#Create a vector in both RA and CLR DFs storing the season information. Here we use astronomical seasons
sharemed_clr_melt$season<-ifelse(sharemed_clr_melt$mm_num%in%c("3","4","5"),
                                 "Spring",
                                 ifelse(sharemed_clr_melt$mm_num%in%c("6","7","8"),
                                        "Summer",
                                        ifelse(sharemed_clr_melt$mm_num%in%c("9","10","11"),
                                               "Autumn",
                                               ifelse(sharemed_clr_melt$mm_num%in%c("12","1","2"),
                                                      "Winter",NA))))

sharemed_rel_melt$season<-ifelse(sharemed_rel_melt$mm_num%in%c("3","4","5"),
                                 "Spring",
                                 ifelse(sharemed_rel_melt$mm_num%in%c("6","7","8"),
                                        "Summer",
                                        ifelse(sharemed_rel_melt$mm_num%in%c("9","10","11"),
                                               "Autumn",
                                               ifelse(sharemed_rel_melt$mm_num%in%c("12","1","2"),
                                                      "Winter",NA))))

#Put relative abundance in 0-100% range
sharemed_rel_melt$Abundance<-sharemed_rel_melt$Abundance*100

#Periodicity statistics

length(unique(lomb_tab$asv)) #1872 seasonal ASVs out of 13755
min(lomb_tab$peak)#11.06
max(lomb_tab$peak)#62.14
#density plot
lomb_tab%>%
  ggplot()+
  geom_histogram(aes(x=peak))

periodicity_stats<-sharemed_rel_melt%>%
  filter(Abundance>0)%>%
  group_by(Sample,mm_text,area,depth,periodic)%>%
  summarise(num_seasonal_asvs=n_distinct(OTU),relab=sum(Abundance))%>%
  mutate(perc_periodic_asvs=num_seasonal_asvs/sum(num_seasonal_asvs)*100)

write.csv(periodicity_stats,"periodicity_stats.csv")

#Display seasonal % of periodic ASVs in surface and bottom samples
sharemed_rel_melt%>%
  filter(Abundance>0)%>%
  group_by(Sample,mm_text,area,depth,periodic)%>%
  summarise(count=n_distinct(OTU),relab=sum(Abundance))%>%
  mutate(perc_periodic_asvs=count/sum(count)*100)%>%
  group_by(mm_text,depth,periodic)%>%
  summarise(abundance=mean(relab),
            sd=sd(relab),
            perc_periodic_asvs=mean(perc_periodic_asvs))%>%
  #filter(periodic=="periodic")%>%
  ggplot(aes(x=mm_text,y=abundance))+
  geom_bar(stat = "identity",position = "stack",aes(fill=periodic))+
  #geom_point(aes(col=periodic))+
  facet_wrap(~depth,scales = "free_x")+
  #geom_smooth(method = "gam",
  # formula = y ~ s(x, k =12,bs="cp"),
  # aes(group=periodic,color=periodic))+
  ggtitle("Seasonal % of periodic ASVs")+
  ylab("Relative proportion (%)")+
  xlab("Month")+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5))

###Highlight seasonal patterns with heatmaps
library(pheatmap)
#outsourcing the sample order in Excel
sample_order <- read_csv("sample_order.csv")

#reordering the pseq object based on custom sample ordering
sharemed_clr<-ps_reorder(sharemed_clr,sample_order$Sample)
#Creating a factor in the pseq object to store the season variable
sample_data(sharemed_clr)$season<-ifelse(sample_data(sharemed_clr)$mm_num%in%c("3","4","5"),
                                         "Spring",
                                         ifelse(sample_data(sharemed_clr)$mm_num%in%c("6","7","8"),
                                                "Summer",
                                                ifelse(sample_data(sharemed_clr)$mm_num%in%c("9","10","11"),
                                                       "Autumn",
                                                       ifelse(sample_data(sharemed_clr)$mm_num%in%c("12","1","2"),
                                                              "Winter",NA))))

#Since the number of seasonal ASVs is very high, we did a subset of the strongly recurrent ones (i.e., PN>40)
strong_recurrent<-lomb_tab%>%filter(peak>=40)%>%select(asv)
#Subset a CLR-transformed pseq object with only the strongly recurrent ASVs
strong_seas_pseq<-prune_taxa(strong_recurrent$asv,sharemed_clr)
#create a matrix
matrix_strong_seas<-t(otu_table(strong_seas_pseq))

#Annotations for the heatmap

#Seasons
season_order<-c("Autumn","Winter","Spring","Summer")
season_anno<-sample_data(sharemed_clr) %>% data.frame() %>% 
  select(season)  %>% 
  rownames_to_column() %>% 
  mutate(season = factor(season, levels = season_order,ordered = TRUE)) %>% 
  column_to_rownames( )
season_pal<-c("#f44336","#2196f3","#8bc34a","#ffeb3b")
names(season_pal)<-season_order

#Taxonomy. We only display the Genera appearing more than 15 times. The others will be regarded as "Others".
tax_anno<- as(tax_table(strong_seas_pseq), 'matrix') %>% as_tibble(rownames = 'asv')  %>% 
  mutate(family = fct_lump(Family, n = 10)) %>% 
  mutate(genus = fct_lump(Genus, n = 15)) %>% 
  select(asv,genus)%>%
  column_to_rownames(var = 'asv') 

#Create a palette for those Genera and naming with the 15 entries plus Other
genus_pal<-c("#FBE183FF", "#F4C40FFF",
                        "#FE9B00FF", "#D8443CFF",
                        "#9B3441FF", "#DE597CFF", "#E87B89FF",
                        "#E6A2A6FF", "#AA7AA1FF", "#9F5691FF",
                        "#633372FF", "#1F6E9CFF", "#2B9B81FF",
                        "#92C051FF", "#2196f3", "#607d8b","grey60")
                        names(genus_pal)<-unique(tax_anno$genus)
                        
                        #Put together the annotation colors
                        anno_colors<-list(season=season_pal,
                                          genus=genus_pal)
                        
                        
                        pheatmap(matrix_strong_seas,
                                 cluster_cols = F,
                                 cluster_rows = T,
                                 annotation_col = season_anno,
                                 annotation_row = tax_anno,
                                 annotation_colors = anno_colors,
                                 gaps_col=c(39,78,114),
                                 cutree_rows = 3)

##Rhytmicity in clusters

#data frame storing ASVs belonging to clusters
indsp_res_clusters_filt
#DF with RA of ASVs >1%
rel1_melt

#Check how many of the ASVs significantly defining one of the three clusters are seasonal
length(indsp_res_clusters_filt$ASVs%in%lomb_tab$asv)
#All of them! 71

#heatmap of the abundance of the 71 within cluster ASV
seas_clust_ps<-prune_taxa(indsp_res_clusters_filt$ASVs,sharemed_clr)

seas_clust_ps<-ps_reorder(seas_clust_ps,sample_order$Sample)
seas_clust_mat<-t(otu_table(seas_clust_ps))

#Taxa annotations
#Keep only those taxa which are discussed in the text
taxa_to_keep_clust<-c("Synechococcus CC9902","Cyanobium PCC-63072",
                      "NS4 marine group","NS5 marine group","HIMB11",
                      "Candidatus Nitrosopumilus","SUP05 cluster",
                      "unassigned_Marine Group II","Fluviicola",
                      "Ascidiaceihabitans","Vibrio")

tax_anno_cluster<- as(tax_table(seas_clust_ps), 'matrix') %>% as_tibble(rownames = 'asv')  %>% 
  mutate(family = fct_lump(Family, n = 10)) %>% 
  mutate(genus = fct_other(Genus, taxa_to_keep_clust )) %>% 
  select(asv,genus)  %>% 
  column_to_rownames(var = 'asv') 

tax_anno_cluster$ASVs<-rownames(tax_anno_cluster)
tax_anno_cluster<-merge(tax_anno_cluster,indsp_res_clusters_filt[,1:2],by="ASVs")
rownames(tax_anno_cluster)<-tax_anno_cluster$ASVs
tax_anno_cluster<-tax_anno_cluster%>%select(!ASVs)
unique(tax_anno_cluster$genus)

#Rename clusters 
tax_anno_cluster$Cluster<-ifelse(tax_anno_cluster$Cluster=="dark","Low light",
                                 ifelse(tax_anno_cluster$Cluster=="apr","Spring",
                                        ifelse(tax_anno_cluster$Cluster=="hot","High temperature",NA)))

#Define a palette
genus_clust_pal<-c("#607d8b"  ,"#FBE183FF",
                              "#FE9B00FF","#D8443CFF",
                              "#9B3441FF","#1F6E9CFF",
                              "#E6A2A6FF","#AA7AA1FF",
                              "#633372FF","#2B9B81FF",
                              "#2196f3")
                              #Rename the palette
names(genus_clust_pal)<-unique(tax_anno_cluster$genus)

#The same for clusters annotation
cluster_pal<-c("#546e7a","#8bc34a","#f44336")
names(cluster_pal)<-c("Low light","Spring","High temperature")
#Put together
anno_cols_clust<-list(season=season_pal,
                      genus=genus_clust_pal,
                      Cluster=cluster_pal)

pheatmap(seas_clust_mat,
         cluster_cols = F,
         cluster_rows = T,
         annotation_col = season_anno,
         annotation_row = tax_anno_cluster,
         annotation_colors = anno_cols_clust,
         gaps_col=c(39,78,114),
         cutree_rows = 4)