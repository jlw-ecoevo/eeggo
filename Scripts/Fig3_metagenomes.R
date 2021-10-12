# JLW - 2021
# Plot BIOGEOTRACES metagenome predictions from gRodon v2.0.0 

# Load Packages + Helper Functions ---------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(maps)
library(ggmap)
library(data.table)
library(randomForest)
library(missForest)
library(glmnet)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,
                                    by=.EACHI,allow.cartesian=TRUE)),
                       stringsAsFactors=F))
}

# Load BIOGEOTRACES growth predictions -----------------------------------------

setwd("~/eeggo/Data")
load("CodonStatistics_BIOGEOTRACES.RData")
bio_df <- bio_df %>% subset(nHE>=10)

# add BIOGEOTRACES metadata ----------------------------------------------------

setwd("~/eeggo/Data")
met_df <- read.csv("GEOTRACES_metadata.csv")
met_df$ffn <- paste0(met_df$NCBI.SRA.Accession.for.assembled.metagenome.contigs,
                     "_",
                     met_df$Sample.name,
                     ".ffn")
write.table(met_df %>% subset(select=c(ffn,NCBI.SRA.Accession)),
            col.names = F,
            row.names = F,
            quote = F,
            file = "srz_to_srr.tsv",
            sep="\t")

met_df$Lat <- met_df$Latitude.and.Longitude %>%
  gsub(pattern = " .*",replace= "") %>%
  as.numeric()
met_df$Lat[grepl("S",met_df$Latitude.and.Longitude)] <- 
  -1*met_df$Lat[grepl("S",met_df$Latitude.and.Longitude)]
met_df$Lon <- met_df$Latitude.and.Longitude %>%
  gsub(pattern = ".*N ",replace= "") %>%
  gsub(pattern = ".*S ",replace= "") %>%
  gsub(pattern = " .*",replace= "") %>%
  as.numeric()
met_df$Lon[grepl("W",met_df$Latitude.and.Longitude)] <- 
  -1*met_df$Lon[grepl("W",met_df$Latitude.and.Longitude)]
world_map <- map_data("world")

bio_df$NCBI.SRA.Accession.for.assembled.metagenome.contigs <- 
  bio_df$Assembly %>%
  gsub(pattern="metaeuk",replace="") %>%
  gsub(pattern="_.*",replace="")
plot_df <- merge.easy(bio_df,
                      met_df,
                      key="NCBI.SRA.Accession.for.assembled.metagenome.contigs")


# plot_df <- plot_df %>% subset(Cruise.series=="GEOTRACES") 
plot_df_srf <- plot_df %>% 
  subset(Depth..m.<100) %>% 
  subset(Cruise.series=="GEOTRACES")


# Plot BIOGEOTRACES  Growth ----------------------------------------------------

ppcr <- ggplot(plot_df_srf) +
  geom_polygon(data=world_map,
               aes(x=long,
                   y = lat,
                   group = group),
               fill="black",
               color="black") +
  geom_point(data=plot_df_srf,
             aes(x=Lon,
                 y=Lat,
                 size=log(2)/d.madin,
                 color=Cruise.ID)) +
  theme_void() +
  scale_color_manual(values=c("#332288",
                              "#88CCEE",
                              "#44AA99",
                              "#117733",
                              "#999933",
                              "#DDCC77",
                              "#CC6677",
                              "#882255")) +
  labs(color="Cruise", size="Max. Growth Rate") +
  theme(legend.position = c(0.1,0.66),
        legend.box.background = element_rect(fill = "white",color="white"),
        plot.title = element_text(vjust=1,hjust=0.05)) +
  guides(fill = guide_legend(order = 2),
         col = guide_legend(order = 1,override.aes = list(size=4))) +
  ggtitle("Prokaryotes")

plot_df_srf_pp <- plot_df_srf %>% 
  group_by(Latitude.and.Longitude) %>%
  filter(d.madin==min(d.madin))
plot_df_srf_pp$rate.madin <- log(2)/plot_df_srf_pp$d.madin
pp <- ggplot(plot_df_srf_pp) +
  geom_polygon(data=world_map,
               aes(x=long,
                   y = lat,
                   group = group),
               fill="black",
               color="black") +
  geom_point(data=plot_df_srf_pp,
             aes(x=Lon,
                 y=Lat,
                 size=rate.madin,
                 color=rate.madin) )+
  theme_void() +
  scale_size_continuous(limits=c(0.05,0.5), breaks=c(0.1,0.2,0.3,0.4)) +
  scale_color_gradient(low="#41ab5d",
                       high="#00441b",
                       guide = "legend",
                       limits=c(0.05,0.5), 
                       breaks=c(0.1,0.2,0.3,0.4)) +
  labs(color="Max. Growth Rate", size="Max. Growth Rate") +
  theme(legend.position = c(0.1,0.66),
        legend.box.background = element_rect(fill = "white",color="white"),
        plot.title = element_text(vjust=1,hjust=0.05)) +
  ggtitle("Prokaryotes")

pecr <- ggplot(plot_df_srf) +
  geom_polygon(data=world_map,
               aes(x=long,
                   y = lat,
                   group = group),
               fill="black",
               color="black") +
  geom_point(data=plot_df_srf,
             aes(x=Lon,y=Lat,
                 size=log(2)/d.eukaryotes,
                 color=Cruise.ID)) +
  theme_void() +
  scale_color_manual(values=c("#332288",
                              "#88CCEE",
                              "#44AA99",
                              "#117733",
                              "#999933",
                              "#DDCC77",
                              "#CC6677",
                              "#882255")) +
  scale_size(breaks=c(0.025,0.075,0.125)) +
  labs(color="Cruise", size="Max. Growth Rate") +
  theme(legend.position = c(0.1,0.67),
        legend.box.background = element_rect(fill = "white",color="white"),
        plot.title = element_text(vjust=1,hjust=0.05)) +
  guides(fill = guide_legend(order = 2),
         col = guide_legend(order = 1,override.aes = list(size=4))) +
  ggtitle("Eukaryotes")

plot_df_srf_pe <- plot_df_srf %>% 
  group_by(Latitude.and.Longitude) %>%
  filter(d.eukaryotes==min(d.eukaryotes))
plot_df_srf_pe$rate.eukaryotes <- log(2)/plot_df_srf_pe$d.eukaryotes
pe <- ggplot(plot_df_srf_pe) +
  geom_polygon(data=world_map,
               aes(x=long,
                   y = lat,
                   group = group),
               fill="black",
               color="black") +
  geom_point(data=plot_df_srf_pe,
             aes(x=Lon,y=Lat,
                 size=rate.eukaryotes,
                 color=rate.eukaryotes)) +
  theme_void() +
  scale_size_continuous(limits=c(0.01,0.2), breaks=c(0.025,0.05,0.075,0.1,0.125,0.15)) +
  scale_color_gradient(low="#41ab5d",
                       high="#00441b",
                       guide = "legend",
                       limits=c(0.01,0.2), 
                       breaks=c(0.025,0.05,0.075,0.1,0.125,0.15)) +
  labs(color="Max. Growth Rate", size="Max. Growth Rate") +
  theme(legend.position = c(0.1,0.66),
        legend.box.background = element_rect(fill = "white",color="white"),
        plot.title = element_text(vjust=1,hjust=0.05)) +
  ggtitle("Eukaryotes")


plot_df$depth <- factor(plot_df$Depth..m.<=100,
                        levels=c("TRUE","FALSE"),
                        labels=c("Depth < 100m","Depth > 100m"))
pd1 <- ggplot(plot_df %>% subset(Depth..m.<=100),
             aes(x=d.madin,y=d.eukaryotes,fill=Depth..m.)) + 
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3,pch=21) + 
  theme_pubclean() +
  theme(legend.position = c(0.1,0.8)) +
  labs(fill="Depth") +
  xlab("Min. Doubling Time Prokaryotes") +
  ylab("Min. Doubling Time Eukaryotes") + 
  scale_fill_gradient(low="white",high="black",limits=c(0,500)) +
  ggtitle("Depth < 100 meters")
pd2 <- ggplot(plot_df %>% subset(Depth..m.>100 & Depth..m. <500),
              aes(x=d.madin,y=d.eukaryotes,fill=Depth..m.)) + 
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3,pch=21) + 
  theme_pubclean() +
  theme(legend.position = "none") +
  labs(fill="Depth") +
  xlab("Min. Doubling Time Prokaryotes") +
  ylab("Min. Doubling Time Eukaryotes") + 
  scale_fill_gradient(low="white",high="black",limits=c(0,500)) +
  ggtitle("Depth > 100 meters")
pd <- ggarrange(pd1,pd2,nrow=2,labels=c("(a)","(b)"))

setwd("~/eeggo/Figures")
pdf("BIOGEOTRACES_panels.pdf",width=16,height=12)
ggarrange(ggarrange(pd1,pd2,nrow=2,labels=c("(a)","(b)")),
          ggarrange(pp,pe,nrow=2,labels=c("(c)","(d)")),
          ncol=2,
          widths=c(2,4))
dev.off()

setwd("~/eeggo/Figures")
pdf("BIOGEOTRACES_cruises.pdf",width=10,height=12)
ggarrange(ppcr,pecr,nrow=2,labels=c("(a)","(b)"))
dev.off()

summary(lm(d.madin~d.eukaryotes*Depth..m.,data=plot_df))
cor.test(plot_df[plot_df$Depth..m.<=100,]$d.madin,
         plot_df[plot_df$Depth..m.<=100,]$d.eukaryotes)
cor.test(plot_df[plot_df$Depth..m.<=100,]$d.madin,
         plot_df[plot_df$Depth..m.<=100,]$d.eukaryotes)$p.value
cor.test(plot_df[plot_df$Depth..m.>100,]$d.madin,
         plot_df[plot_df$Depth..m.>100,]$d.eukaryotes)

# Plot CUB ---------------------------------------------------------------------


pd1cub <- ggplot(plot_df %>% subset(Depth..m.<=100),
              aes(x=CUBHE.madin,y=CUBHE.eukaryotes,fill=Depth..m.)) + 
  # scale_x_log10() +
  # scale_y_log10() +
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3,pch=21) + 
  theme_pubclean() +
  theme(legend.position = c(0.8,0.8)) +
  # theme(legend.position = "right") +
  labs(fill="Depth") +
  xlab("CUB Ribosomal Proteins, Prokaryotes") +
  ylab("CUB Ribosomal Proteins, Eukaryotes") + 
  scale_fill_gradient(low="white",high="black",limits=c(0,500)) +
  ggtitle("Depth < 100 meters")
pd2cub <- ggplot(plot_df %>% subset(Depth..m.>100 & Depth..m. <500),
              aes(x=CUBHE.madin,y=CUBHE.eukaryotes,fill=Depth..m.)) + 
  # scale_x_log10() +
  # scale_y_log10() +
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3,pch=21) + 
  theme_pubclean() +
  theme(legend.position = "none") +
  # theme(legend.position = "right") +
  labs(fill="Depth") +
  xlab("CUB Ribosomal Proteins, Prokaryotes") +
  ylab("CUB Ribosomal Proteins, Eukaryotes") + 
  scale_fill_gradient(low="white",high="black",limits=c(0,500)) +
  ggtitle("Depth > 100 meters")
pdcub <- ggarrange(pd1cub,pd2cub,ncol=2,labels=c("(a)","(b)"))
pdcub


setwd("~/eeggo/Figures")
pdf("BIOGEOTRACES_euk_vs_prok_CUB.pdf",width=10,height=6)
pdcub
dev.off()


cor.test(plot_df[plot_df$Depth..m.<=100,]$CUBHE.madin,
         plot_df[plot_df$Depth..m.<=100,]$CUBHE.eukaryotes)
cor.test(plot_df[plot_df$Depth..m.>100,]$CUBHE.madin,
         plot_df[plot_df$Depth..m.>100,]$CUBHE.eukaryotes)


# Supplemental Depth Figure ----------------------------------------------------

plot_df$Collection.Date <- plot_df$Collection.Date %>% 
  gsub(pattern="T.*",replace="")
pdf <- plot_df %>% group_by(depth,Collection.Date) %>% 
  summarise(d.madin=mean(d.madin),d.eukaryotes=mean(d.eukaryotes)) %>%
  as.data.frame(stringsAsFacotrs = F)
pdfw <- reshape(pdf, direction = "wide", 
                timevar = "depth", 
                idvar = "Collection.Date")


p1dep <- ggplot(pdfw, aes(x=`d.madin.Depth < 100m`,
                          y=`d.madin.Depth > 100m`)) + 
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3) + 
  theme_pubclean() +
  xlab("Min. Doubling Time Prokaryotes, Depth < 100m") +
  ylab("Min. Doubling Time Prokaryotes, Depth > 100m") +
  xlim(0,15) +
  ylim(0,15)
p2dep <- ggplot(pdfw, aes(x=`d.madin.Depth < 100m`,
                          y=`d.eukaryotes.Depth > 100m`)) + 
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3) + 
  theme_pubclean() +
  xlab("Min. Doubling Time Prokaryotes, Depth < 100m") +
  ylab("Min. Doubling Time Eukaryotes, Depth > 100m")  +
  xlim(0,15) +
  ylim(0,55)
p3dep <- ggplot(pdfw, aes(x=`d.eukaryotes.Depth < 100m`,
                          y=`d.madin.Depth > 100m`)) + 
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3) + 
  theme_pubclean() +
  xlab("Min. Doubling Time Eukaryotes, Depth < 100m") +
  ylab("Min. Doubling Time Prokaryotes, Depth > 100m")  +
  xlim(10,40) +
  ylim(0,15)
p4dep <- ggplot(pdfw, aes(x=`d.eukaryotes.Depth < 100m`,
                          y=`d.eukaryotes.Depth > 100m`)) + 
  geom_smooth(color="black",fill="black",method="lm") +
  geom_point(size=3) + 
  theme_pubclean() +
  xlab("Min. Doubling Time Eukaryotes, Depth < 100m") +
  ylab("Min. Doubling Time Eukaryotes, Depth > 100m")  +
  xlim(10,40) +
  ylim(0,55)

setwd("~/eeggo/Figures")
pdf("BIOGEOTRACES_depthcomparison.pdf",width=12,height=12)
ggarrange(p1dep,
          p3dep,
          p2dep,
          p4dep,
          ncol=2,
          nrow=2,
          labels=c("(a)","(b)","(c)","(d)"))
dev.off()


cor.test(pdfw$`d.madin.Depth < 100m`,pdfw$`d.madin.Depth > 100m`)
cor.test(pdfw$`d.madin.Depth < 100m`,pdfw$`d.eukaryotes.Depth > 100m`)
cor.test(pdfw$`d.eukaryotes.Depth < 100m`,pdfw$`d.madin.Depth > 100m`)
cor.test(pdfw$`d.eukaryotes.Depth < 100m`,pdfw$`d.eukaryotes.Depth > 100m`)
