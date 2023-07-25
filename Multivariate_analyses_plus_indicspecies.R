#Packages loading list

library(phyloseq)
library(microbiome)
library(microViz)
library(fantaxtic)
library(ggrepel)
library(ggpubr)
library(forcats)
library(indicspecies)
library(stringr)
library(tidyverse)
library(viridis)
library(viridisLite)
library(vegan)
library(Hmisc)
library(corrplot)
library(pheatmap)
library(svglite)
library(pals)

#Inventory
sharemed_noint    #original pseq object bearing only surface and bottom samples
metadata          #metadata extracted from the pseq object sharemed_noint
sharemed_noint_NA #samples with NA in temp, sal, chla, doc and nuts have been removed
sharemed_noint_rel_NA #the same as above but in relative proportions (0-1)


#loading files
Sharemed<-
  readRDS("G:/Shared drives/Sharemed_Joint-ITA-SLO/Sharemed_Joint-ITA-SLO/RstudioProject_Sharemed/Sharemed.rds")

summarize_phyloseq(Sharemed)
rowSums(Sharemed@otu_table@.Data)
metadata_sharemed<-sample_data(Sharemed)
View(metadata_sharemed)

####Remove intermediate samples (5m,10m,15m)####
#Including just "5m" since the pattern would be detected also in strings containing "15m"
pattern<-c("5m","10m")
#apply the str_subset functions to all the elements in the sample names vector
samples_to_remove<-map(pattern,str_subset,string=metadata_sharemed$MySampleName)
#the result of this function is a list. Unlist it
str(samples_to_remove)
samples_to_remove<-unlist(samples_to_remove)
###Subset the phyloseq object to remove intermediate samples###
#This line will return a pseq object containing the samples to remove
#Just a double check that everything is ok. 36 elements in the vector "sample_to_remove" and
#36 samples in the resulting phyloseq object
subset_samples(Sharemed, MySampleName %in% samples_to_remove)
#We can proceed to subset
sharemed_noint<-subset_samples(Sharemed, !(MySampleName %in% samples_to_remove))
#Store the new pseq object 
saveRDS(sharemed_noint,"sharemed_no_intermediate.rds")
#Extract again the metadata from the pseq object
metadata_sharemed<-sample_data(sharemed_noint)

####Replace sample data with updated data####
Sharemed_envdata_20032023 <- read_csv("G:/Shared drives/Sharemed_Joint-ITA-SLO/Sharemed_Joint-ITA-SLO/RstudioProject_Sharemed/Vincenzo/New folder/sharemed_clean/data/sharemed_envdata_20032023.csv")
#remove intermediate samples
pattern<-c("5m","10m")
samples_to_remove<-map(pattern,str_subset,string=Sharemed_envdata_20032023$sample_name)
samples_to_remove<-unlist(samples_to_remove)
metadata<-Sharemed_envdata_20032023[!(Sharemed_envdata_20032023$sample_name%in%samples_to_remove),]
#let's do some cleaning
rm(pattern)
rm(samples_to_remove)
rm(Sharemed)
rm(Sharemed_envdata_20032023)
#define some factors
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
#Check if everithyng is okay
View(metadata)
#put it in the phyloseq object
metadata<-as.data.frame(metadata)
rownames(metadata)<-metadata$sample_name
sample_data(sharemed_noint)<-metadata
View(metadata)
#cleanup
rm(metadata_sharemed)

####Revise taxonomy table####
#If a cell contains NA, this loop replace the NA with the prefix "unassigned_" followed by the previously known taxonomic level
#That is why we bothered to replace all the strange text with NA, so we have all formatted in the right way 
#and avoid problems in merging by taxonomic levels

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

View(tax_tab_torevise)
str(tax_tab_torevise)

#Replace taxonomy table
tax_table(sharemed_noint)<-tax_tab_torevise
#by removing some samples, there are probably a number of taxa which sum up to zero because they were only present in those samples removed.
#check this
any(taxa_sums(sharemed_noint)==0)
#ok, this is TRUE. Remove those taxa
sharemed_noint<-prune_taxa(taxa_sums(sharemed_noint)>0,sharemed_noint) #before 13790, after 13742


#There are a bunch of samples wich have NAs in the metadata.
#Create a new phyloseq object without those samples to proceed with the RDA
sharemed_noint_NA<-ps_drop_incomplete(sharemed_noint, 
                                      vars = c("temp","sal","doc","don",
                                               "tdn","ammon","nitrites",
                                               "nitrates","phosphate","chla"), 
                                      verbose = FALSE)

####Ordinations####
###Collinearity checking, scaling and renaming variables, creating a new Pseq object with scaled metadata for ordination

#Create a df storing env variables used in multivariate analisys
env_data<-sample_data(sharemed_noint_NA)[,c(1,2,5,6,8,10,4,11,29,30,32,34:40)]
#check for collinearity
cormat<-cor(env_data[,8:18])
corrplot((cormat),"circle")
#Temperature vs. DOC and TDN vs DON are highly correlated.
#For the analysis we will keep only temperature and DON
#remove DOC and TDN
colnames(env_data)
env_data<-env_data[,-c(12,18)]
#Scale environmental variables
env_data_scaled<-decostand(env_data[,8:16],method = "standardize")
#identify scaled variable appending the suffix "scaled" in their colnames
newnames<-colnames(env_data_scaled)
newnames<-paste(newnames,"scaled",sep="_")
newnames
colnames(env_data_scaled)<-newnames
#bind these with the original metadata
#check if rownames are in the same order prior binding
all.equal(rownames(env_data_scaled),rownames(sample_data(sharemed_noint_NA)))
sample_data(sharemed_noint_NA)[,44:59]<-env_data_scaled
##Cleanup
rm(cormat)

sample_data(sharemed_noint_NA)

#CAP ordination (RDA) with scaled environmental variables
cap_ord_scaled<-ordinate(physeq = sharemed_noint_NA,  
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
#Test the model
cap_ord_scaled$
  rda(cap_ord_scaled$Ybar)
model_rda<-anova.cca(cap_ord_scaled,step=1000,by="margin")
100-sum(model_rda$`F`[1:9])
RsquareAdj(cap_ord_scaled)
#           Df SumOfSqs      F Pr(>F)    
#Model      9   17.489 11.148  0.001 ***
#Residual 132   23.008   
#Model is significant, keep it

###Subset Pseq object to keep only taxa with a relative abundance >1% in at least one sample
#Relative proportions
#Transform in relative abundance (0-1 range)
sharemed_noint_rel_NA<-microbiome::transform(sharemed_noint_NA,"compositional")
##Filtering out rare taxa
#1% in at least one sample
taxa_greater1_NA<-filter_taxa(sharemed_noint_rel_NA,function(x) max(x)>0.01,TRUE)
taxa_greater1_NA
#26 taxa against 13790
#Give a look at how much of the community (in proportion) are we retaining with this filter
mean(rowSums(taxa_greater1_NA@otu_table@.Data))#78.41%
max(rowSums(taxa_greater1_NA@otu_table@.Data))#94.84%
min(rowSums(taxa_greater1_NA@otu_table@.Data))#61.17%
#on average, 12% of the community is composed of ASVs that have a relative abundance lower than 1%


###CAP ordination with Pseq subset t 1%
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
#Test the model

bc_rda<-distance(taxa_greater1_NA,"bray")
str(env_data_scaled)

summary(cap_ord_scaled_1)

rda(bc_rda~ day_length_scaled+
      temp_scaled+
      sal_scaled+
      chla_scaled+
      ammon_scaled+
      nitrites_scaled+
      nitrates_scaled+
      phosphate_scaled+
      don_scaled,env_data_scaled)

model_rda<-anova.cca(cap_ord_scaled_1,by="margin",step=1000)
model_rda
sum(model_rda$SumOfSqs)

#           Df SumOfSqs      F Pr(>F)    
#Model      9   18.702 14.965  0.001 ***
#Residual 132   18.329   
#Model is significant, keep it

#Plot
cap_plot_scaled_rel_1<-plot_ordination(physeq = taxa_greater1_NA, 
                                       ordination = cap_ord_scaled_1, 
                                       color = "mm_text", 
                                       shape = "depth",
                                       axes = c(1,2)) + 
  geom_point(aes(colour = mm_text), size = 5) + 
  geom_point(size = 3)+
  coord_fixed()+
  scale_color_discrete()
cap_plot_scaled_rel_1 + ggtitle("CAP_Plot")

# Now add the environmental variables as arrows 
arrowmat_scaled_rel_1<- vegan::scores(cap_ord_scaled_1, display = "bp")
# Add labels, make a data.frame 
arrowdf_scaled_rel_1<- data.frame(labels_scaled_rel_1=c("Day Length","Temperature","Salinity",
                                                        "Chl a","Ammonia","Nitrites",
                                                        "Nitrates","Phosphate","DON"), 
                                  arrowmat_scaled_rel_1)
# Define the arrow aesthetic mapping 
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

##now plot the arrow
cap_plot_scaled_rel_1<-cap_plot_scaled_rel_1 + 
  scale_color_viridis_d(direction = -1)+
  geom_segment(mapping = arrow_map_scaled_rel_1,
               size = .7, 
               data = arrowdf_scaled_rel_1,
               color = "black", 
               arrow = arrowhead_scaled_rel_1)+ 
  geom_text(mapping = label_map_scaled_rel_1,
            size = 5,  
            data = arrowdf_scaled_rel_1, 
            show.legend = F) + 
  ggtitle("dbRDA/bray/relab/1%/env_scaled")+
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
        legend.position = "bottom",
  )+
  labs(color="Month",
       shape="Depth")+
  coord_fixed()
###RDA final plot###  
cap_plot_scaled_rel_1
ggsave(file="RDA.svg",plot=cap_plot_scaled_rel_1)

####PCOA###
#use all the samples, environmental variables are not needed so we can use also those samples wich have NAs
#Transform in relative abundance (0-1 range)
sharemed_noint_rel<-microbiome::transform(sharemed_noint,"compositional")
##Filtering out rare taxa
#1% in at least one sample
taxa_greater1<-filter_taxa(sharemed_noint_rel,function(x) max(x)>0.01,TRUE)
#Stay at >1% abundance treshold
#BC dissimilarities
BC_rel1<-  phyloseq::distance(taxa_greater1, method = "bray")
#PCoA
pcoa_rel1<-ordinate(taxa_greater1,
                    method = "PCoA",
                    distance = BC_rel1)
plot_ordination(taxa_greater1,pcoa_rel1,
                type="samples",
                shape="depth",
                color="mm_text")+
  geom_point(size=5)+
  ggtitle("PCoA/relab/>1%")+
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
  )+
  labs(color="Month",
       shape="Depth")+
  coord_fixed()
#On this plot we identify the most distant points to define the three clusters
##indicator species analysis
#1-factors are defined based on bray curtis dissmilarity percentiles, starting from the most
#extreme points identifying the "triangle" on the PCoA plot. 

#define clusters
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

#create a factor in the pseq metadata storing the classification of samples in the three clusters
sample_data(taxa_greater1)$cluster<-ifelse(sample_data(taxa_greater1)$sample_name%in%rel1_apr1,
                                           "apr",
                                           NA)
sample_data(taxa_greater1)$cluster<-ifelse(sample_data(taxa_greater1)$sample_name%in%rel1_hot1,
                                           "hot",
                                           sample_data(taxa_greater1)$cluster)
sample_data(taxa_greater1)$cluster<-ifelse(sample_data(taxa_greater1)$sample_name%in%rel1_dark1,
                                           "dark",
                                           sample_data(taxa_greater1)$cluster)
#transform as factor
sample_data(taxa_greater1)$cluster<-factor(sample_data(taxa_greater1)$cluster)
#colors for factors
grouping_cols<-c("hot"="#E74C3C",
                 "dark"="#3498DB",
                 "apr"="#1ABC9C")

names(grouping_cols)
plot_ordination(taxa_greater1,
                pcoa_rel1,
                type="samples",
                color="cluster")+
  geom_point(size=6,aes(shape=depth))+
  scale_color_manual(values=grouping_cols,labels=c("Spring",
                                                   "Low Light",
                                                   "High Temperature"))+
  scale_shape_discrete(labels=c("Surface","Bottom"))+
  #ggtitle("PCoA/relab/>1%")+
  geom_text_repel(label =taxa_greater1@sam_data$mm_num, 
                  size = 4, 
                  nudge_y =-0.025, 
                  color = "black")+
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
        legend.position = "bottom")+
  coord_fixed()
####Indicator species analysis####
#HEre the goal is to identify ASVs associated with each group of samples
#We will do this by checking each cluster against all the other samples
#Replace all the NAs (i.e., samples that have not been assigned to any of the three clusters) with "all"
sample_data(taxa_greater1)$cluster<-ifelse(is.na(sample_data(taxa_greater1)$cluster),
                                           "all",
                                           as.character(sample_data(taxa_greater1)$cluster))
#Visual check that everything is in order
View(sample_data(taxa_greater1))

#create three new clusters to check each level (i.e., apr, hot, and dark) against all the other samples
sample_data(taxa_greater1)$cl_apr<-ifelse(sample_data(taxa_greater1)$cluster=="apr",
                                          "apr",
                                          "all")
sample_data(taxa_greater1)$cl_hot<-ifelse(sample_data(taxa_greater1)$cluster=="hot",
                                          "hot",
                                          "all")
sample_data(taxa_greater1)$cl_dark<-ifelse(sample_data(taxa_greater1)$cluster=="dark",
                                           "dark",
                                           "all")
#Visual check that everything is in order
View(sample_data(taxa_greater1))

#prepare for indicator species analysis
#community matrix
com_mat_rel1_vs<-as.matrix(otu_table(taxa_greater1))
#vector of partitioning
clustvec_apr<-sample_data(taxa_greater1)[,"cl_apr"]
names_apr<-rownames(clustvec_apr)
clustvec_apr<-factor(clustvec_apr$cl_apr)
names(clustvec_apr)<-names_apr

clustvec_hot<-sample_data(taxa_greater1)[,"cl_hot"]
names_hot<-rownames(clustvec_hot)
clustvec_hot<-factor(clustvec_hot$cl_hot)
names(clustvec_hot)<-names_hot

clustvec_dark<-sample_data(taxa_greater1)[,"cl_dark"]
names_dark<-rownames(clustvec_dark)
clustvec_dark<-factor(clustvec_dark$cl_dark)
names(clustvec_dark)<-names_dark

indics_aprvsall<-multipatt(com_mat_rel1_vs,clustvec_apr,func = "IndVal.g",
                           duleg=TRUE,control = how(nperm=9999))
indics_hotvsall<-multipatt(com_mat_rel1_vs,clustvec_hot,func = "IndVal.g",
                           duleg=TRUE,control = how(nperm=9999))
indics_darkvsall<-multipatt(com_mat_rel1_vs,clustvec_dark,func = "IndVal.g",
                            duleg=TRUE,control = how(nperm=9999))
#generate data frames storing the result of the analysis
#April cluster vs all samples
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

#Hot cluster vs all samples
rm(indsp_net_hot)
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

#dark cluster vs all samples
rm(indsp_net_dark)
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


indsp_res_clusters<-rbind(indsp_net_apr[which(indsp_net_apr$Cluster=="apr"),],
                          indsp_net_hot[which(indsp_net_hot$Cluster=="hot"),],
                          indsp_net_dark[which(indsp_net_dark$Cluster=="dark"),])

indsp_res_clusters_filt<-indsp_res_clusters%>%filter(Fidelity>=0.8, Specificity>=0.8)