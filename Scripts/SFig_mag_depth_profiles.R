
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(ggridges)
library(tidyverse)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}


setwd("~/eeggo/Data/")
x <- read.csv("TOPAZ_MAG_abundance.csv")
alexander_met <- read.csv("TableS09_TOPAZ_HetScore.csv")
names(alexander_met)[1] <- "new_mag_name"
profiles <- merge.easy(x,alexander_met,key="new_mag_name")

depths <- read.csv("mag_sample_depths.csv")
sample_to_ers <- read.csv("sample_to_ers.csv")
names(err_to_ers)[1] <- "Station"
err_to_ers <- read.csv("err_to_ers.csv")
names(err_to_ers) <- c("Station","ERS")
depths <- merge.easy(depths,sample_to_ers,key="Sample")
depths <- merge.easy(depths,err_to_ers,key="ERS") %>%
  subset(!is.na(Depth)) %>%
  subset(!is.na(Station))

mag_profiles <- merge.easy(profiles,depths,key="Station")


profile_trophy <- mag_profiles %>% 
  group_by(Depth,PredictedTrophicMode,Station) %>%
  summarise(TPM=sum(TPM))

setwd("~/eeggo/Figures/")
pdf("MAG_depths.pdf",width=5,height=6)
ggplot(profile_trophy,aes(y=TPM,x=Depth,fill=PredictedTrophicMode)) + 
  geom_point(size=3,alpha=0.5,pch=21) +
  geom_smooth(color="black") +
  scale_y_log10() +
  theme_pubclean() + 
  geom_vline(xintercept=100,lty=2) +
  xlab("Depth (Meters)") + 
  labs(color="Predicted Trophic Mode") + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values=c("#4B0092","#1AFF1A"))
dev.off()



