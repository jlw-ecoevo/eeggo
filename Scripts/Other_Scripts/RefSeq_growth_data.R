# JLW - 2021
# Cleanup and put together RefSeq growth data

# Load Packages
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(MASS)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}


setwd("~/eeggo/Data")
load("CodonStatistics_protozoa.RData")
assembly_df$Accession <- assembly_df$Assembly %>% 
  gsub(pattern="[.].*",replace="")
assembly_df <- assembly_df %>% subset(nHE>10)

x <- read.csv("euk_growth_ogt_Weissman_RefSeq.csv",stringsAsFactors = F) 
growth_rates <- data.frame(Accession=character(),
                Species.Name=character(),
                Doubling.Time=numeric(),
                OGT=numeric())
for(i in 1:nrow(x)){
  growth_rates <- rbind(growth_rates,
             data.frame(Accession=x[i,"Accession"] %>% 
                          strsplit(split=";") %>% 
                          unlist() %>% 
                          gsub(pattern="[.].*",replace=""),
                        Species.Name=x[i,"Species.Name"],
                        Doubling.Time=x[i,"Doubling.Time"],
                        OGT=x[i,"OGT"]))
}
growth_df <- merge.easy(assembly_df,growth_rates,key="Accession") %>%
  group_by(Species.Name) %>%
  summarise(Doubling.Time=mean(Doubling.Time),
            OGT=mean(OGT),
            CUBHE=mean(CUBHE))

setwd("~/eeggo/Data")
save(growth_df,file="growth_RefSeq.RData")

protozoa_growth_df <- merge.easy(assembly_df,growth_rates,key="Accession")
write.table(protozoa_growth_df%>%subset(select=c(Assembly,OGT)),
            file="assemblies_temps_protozoa.txt",
            sep="\t",
            row.names=F,
            col.names=F,
            quote=F)

