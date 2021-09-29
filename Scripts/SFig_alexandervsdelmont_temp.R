library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(data.table)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,
                                    by=.EACHI,allow.cartesian=TRUE)),
                       stringsAsFactors=F))
}

setwd("~/eeggo/Data/")
cl <- read.delim("TableS12_Clustered99ANI.tsv",sep="\t",header=T)

cl_df <- data.frame(cluster = character(),
                    mag = character(),
                    stringsAsFactors = F)
for(i in 1:nrow(cl)){
  cl_df <- rbind(cl_df,
                 data.frame(cluster = cl$CLUSTER[i],
                            mag = strsplit(cl$MAGs[i],split=",") %>%
                              unlist() %>% 
                              trimws(),
                            stringsAsFactors = F))
}
cl_df$Source <- "Alexander"
cl_df$Source[grep("TARA",cl_df$mag)] <- "Delmont"


delmont <- read.delim("delmont_temp.tbl", header =F,sep = " ")
alexander <- read.delim("alexander_temp.tbl", header =F,sep = " ")
temps <- rbind(delmont,alexander)
names(temps) <- c("mag","OGT")

cl_ogt <- merge.easy(cl_df,temps,key="mag")

cl_sum <- cl_ogt %>% group_by(cluster,Source) %>% summarise(OGT=mean(OGT,na.rm=T))

cl_wide <- data.frame(cluster = character(),
                    OGT.alexander = character(),
                    OGT.delmont = character(),
                    stringsAsFactors = F)
cl_ids <- unique(cl_sum$cluster)
for(i in 1:length(cl_ids)){
  cl_wide <- rbind(cl_wide,
                   data.frame(cluster = cl_ids[i],
                              OGT.alexander = mean(cl_sum$OGT[cl_sum$Source=="Alexander" & cl_sum$cluster==cl_ids[i]]),
                              OGT.delmont = mean(cl_sum$OGT[cl_sum$Source=="Delmont" & cl_sum$cluster==cl_ids[i]]),
                              stringsAsFactors = F))
}

setwd("~/eeggo/Figures/")
pdf("TempPredCompare.pdf",width=7,height=7)
ggplot(cl_wide,aes(x=OGT.delmont,y=OGT.alexander)) + 
  geom_point() + 
  theme_pubclean() + 
  geom_abline(slope=1,intercept=0,lty=2) +
  xlim(-1,28) +
  ylim(-1,28) +
  xlab("Mean OGT of Cluster from ML Predictions on Delmont et al. MAGs")+
  ylab("Mean OGT of Cluster from 99% Predictions on Alexander et al. MAGs")
dev.off()

cor.test(cl_wide$OGT.alexander,cl_wide$OGT.delmont,na.rm=T)
