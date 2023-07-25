####Relative Abundance plots for each cluster with top 20 ASVs####
#Subset by cluster keeping variable of interest
dark_df<-rel1_melt%>%
  filter(OTU%in%
           indsp_res_clusters_filt[which(indsp_res_clusters_filt$Cluster=="dark"),"ASVs"],
         cl_dark=="dark")%>%
  arrange(desc(Abundance),
          by_group=TRUE)%>%
  mutate(Genus = factor(Genus, 
                        levels = unique(Genus),
                        ordered = TRUE))

apr_df<-rel1_melt%>%
  filter(OTU%in%
           indsp_res_clusters_filt[which(indsp_res_clusters_filt$Cluster=="apr"),"ASVs"],
         cl_apr=="apr")%>%
  arrange(desc(Abundance),
          by_group=TRUE)%>%
  mutate(Genus = factor(Genus, 
                        levels = unique(Genus),
                        ordered = TRUE))

hot_df<-rel1_melt%>%
  filter(OTU%in%
           indsp_res_clusters_filt[which(indsp_res_clusters_filt$Cluster=="hot"),"ASVs"],
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

#Try what comes out with nested palettes
ggnested::ggnested(clust_df,
                   aes_string(main_group="Genus",
                              sub_group="OTU",
                              x="mm_yy",
                              y="Abundance"))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance (%)\n")+
  facet_wrap(cluster~area~depth,scales="free_x",nrow = 3)+
  theme(axis.text.x = element_text(angle = 90,vjust=0.3),
        panel.margin.x=unit(-0.25,"lines"),
        panel.border = element_rect(fill=NA,linewidth = 0.5),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(2,"mm"))


#Create the palette for the top 20 more abundant Genera
View(clust_df)
#Take out the first 20 Genera
top20_genera<-levels(clust_df$Genus)[1:20]

#Group the other Genera as "Others" while keeping the others as they are
clust_df$Genus_top20<-ifelse(clust_df$Genus%in%top20_genera,
                             as.character(clust_df$Genus),"Other")
#Be sure to refactor this to keep the order of appearance in plots
clust_df<-clust_df%>%
  mutate(Genus_top20 = factor(Genus_top20, 
                              levels = unique(Genus_top20),
                              ordered = TRUE))
#Check it
levels(clust_df$Genus_top20)

#Create a factor storing the AVs belonging to the top20 genera. This will be the subgroups for shading
ASV_top20<-unique(clust_df%>%select(OTU,Genus_top20)%>%filter(!Genus_top20=="Other"))
clust_df$ASVtop20<-ifelse(clust_df$OTU%in%ASV_top20$OTU,
                          clust_df$OTU,"Other")

#Manually create a plette with 20 distinct base colors for the 20 genera.
#ASVs within this genera will be assigned shades of the same color
pal_top_20<-c("#eccc68","#ffa502","#7bed9f","#ff7f50","#2ed573","#9980FA",
                       "#009432","#ff4757","#12CBC4","#1289A7","#0652DD","#1B1464",
                       "#FDA7DF","#D980FA","#3742fa","#5758BB","#ff6348","#B53471",
                       "#833471","#6F1E51","#a4b0be")
                       #Name the palette with the Genus top 20 factor
names(pal_top_20)<-levels(clust_df$Genus_top20)


#Labels for facets
depth_label<-c("Surface","Bottom")
names(depth_label)<-c("surface","bottom")

HT_label<-"High temperature"
names(HT_label)<-"hot"

LL_label<-"Low light"
names(LL_label)<-"dark"

SP_label<-"Spring"
names(SP_label)<-"apr"

#High temperature cluster
hot_plot_nested<-
  clust_df%>%filter(cluster=="hot")%>%
  ggnested::ggnested(aes_string(main_group="Genus_top20",
                                sub_group="ASVtop20",
                                x="mm_yy",
                                y="Abundance"),main_palette=pal_top_20)+
  geom_bar(stat="identity")+
  ylab("Relative Abundance (%)\n")+
  facet_grid(cols=vars(area,depth),
             rows = vars(cluster),
             scales = "free_x",
             space = "free",
             labeller = labeller(depth=depth_label,
                                 cluster=HT_label))+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
        axis.text = element_text(size=13),
        axis.title = element_text(size=15),
        strip.text.x = element_text(size=13),
        strip.text.y = element_text(size=15,face="bold"),
        panel.background = element_rect(fill=NA),
        axis.line = element_line(linewidth = 1.5,lineend = "round"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2,"mm"),
        axis.title.x = element_blank())
#Low light cluster
dark_plot_nested<-clust_df%>%filter(cluster=="dark")%>%
  ggnested::ggnested(aes_string(main_group="Genus_top20",
                                sub_group="ASVtop20",
                                x="mm_yy",
                                y="Abundance"),main_palette=pal_top_20)+
  ylim(0,60)+
  xlab("Month_Year")+
  geom_bar(stat="identity")+
  ylab("Relative Abundance (%)\n")+
  facet_grid(cols=vars(area,depth),
             rows=vars(cluster),
             scales = "free_x",
             space = "free",
             labeller = labeller(depth=depth_label,
                                 cluster=LL_label))+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
        axis.text = element_text(size=13),
        axis.title = element_text(size=15),
        strip.text.x = element_text(size=13),
        strip.text.y = element_text(size=15,face="bold"),
        panel.background = element_rect(fill=NA),
        axis.line = element_line(linewidth = 1.5,lineend = "round"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2,"mm"),
        axis.title.x = element_blank())

#Spring cluster
spring_plot_nested<-clust_df%>%filter(cluster=="apr")%>%
  ggnested::ggnested(aes_string(main_group="Genus_top20",
                                sub_group="ASVtop20",
                                x="mm_yy",
                                y="Abundance"),main_palette=pal_top_20)+
  geom_bar(stat="identity")+
  ylab("Relative Abundance (%)\n")+
  facet_grid(cols=vars(area,depth),
             rows=vars(cluster),
             scales = "free_x",
             space = "free",
             labeller = labeller(depth=depth_label,
                                 cluster=SP_label))+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
        axis.text = element_text(size=13),
        axis.title = element_text(size=15),
        strip.text.x = element_text(size=13),
        strip.text.y = element_text(size=15,face="bold"),
        panel.background = element_rect(fill=NA),
        axis.line = element_line(linewidth = 1.5,lineend = "round"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2,"mm"),
        axis.title.x = element_blank())

#Extract the legend
legend_all<-get_legend(clust_df%>%
                         ggnested::ggnested(aes_string(main_group="Genus_top20",
                                                       sub_group="ASVtop20",
                                                       x="Sample",
                                                       y="Abundance"),main_palette=pal_top_20)+
                         geom_bar(stat="identity")+
                         ylab("Relative Abundance (%)\n")+
                         facet_grid(cols=vars(area,depth),rows=vars(cluster),scales = "free_x",space = "free")+
                         theme(axis.text.x = element_text(angle = 45,vjust=1,hjust =1),
                               panel.margin.x=unit(-0.25,"lines"),
                               panel.border = element_rect(fill=NA,linewidth = 0.5),
                               strip.background = element_blank(),
                               strip.text = element_blank(),
                               legend.position = "bottom",
                               legend.key.size = unit(2,"mm")))


#Arrange plots together
barplot_all<-ggarrange(hot_plot_nested,dark_plot_nested,spring_plot_nested,
                       common.legend = F,nrow = 3,align = "hv")
barplot_all


plot(legend_all)