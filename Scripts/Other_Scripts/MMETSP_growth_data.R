# JLW - 2021
# Cleanup and put together MMETSP growth data

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

# After matching up to MMETSP by hand, load data from Thomas et al and Rose and Caron papers
setwd("~/eeggo/Data")
x1 <- read.csv("euk_growth_ogt_Thomas_MMETSP.csv",stringsAsFactors = F) %>%
  subset(select=c("MMETSP.Strain","Growth.rate","OGT"))
x1$Source <- "Thomas et al."
x2 <- read.csv("euk_growth_ogt_Rose_MMETSP.csv",stringsAsFactors = F)   %>%
  subset(select=c("MMETSP.Strain","Growth.rate","OGT"))
x2$Source <- "Rose and Caron"
y <- rbind(x1,x2) %>% subset(!is.na(MMETSP.Strain))
y$Doubling.Time <- 24*log(2)/y$Growth.rate
y$Growth.rate <- NULL

# Load values matched by literature search
x3 <- read.csv("euk_growth_ogt_Weissman_MMETSP.csv",stringsAsFactors = F)   #%>%
names(x3) <- c("MMETSP.Strain","Doubling.Time","OGT","Source")
x3 <- x3[,c(1,3,4,2)]
growth_rates <- rbind(y,x3)

# Extract species names from strain names
growth_rates$Species.Name <- gsub(";.*","",growth_rates$MMETSP.Strain)
growth_rates$Species.Name[!grepl("sp[.]",growth_rates$Species.Name)] <-
  growth_rates$Species.Name[!grepl("sp[.]",growth_rates$Species.Name)] %>%
  strsplit(.,split=" ") %>%
  lapply(.,"[",1:2) %>%
  lapply(.,paste,collapse=" ") %>%
  unlist() %>%
  gsub(pattern=",",replace="")

# Deduplicate entries (one per species)
growth_rates_dedup <- growth_rates %>%
  unique() %>%
  group_by(Species.Name) %>%
  slice(which.min(Doubling.Time))
growth_rates_dedup$MMETSP.Strain <- NULL

# Check Normality
qqnorm(growth_rates_dedup$Doubling.Time %>% log10())
qqline(growth_rates_dedup$Doubling.Time %>% log10())

# Test for outliers
test_outliers <- rosnerTest(growth_rates_dedup$Doubling.Time %>% log10())
outlier_ind <- test_outliers$all.stats$Obs.Num[test_outliers$all.stats$Outlier==T]
growth_rates_dedup_nooutliers <- growth_rates_dedup[-outlier_ind,]

# Save to CSV
write.csv(growth_rates_dedup,"MMETSP_all_growth_ogt_clean.csv")
write.csv(growth_rates_dedup_nooutliers,"MMETSP_all_growth_ogt_clean_nooutliers.csv")

# add codon usage
load("CodonStatistics_MMETSP.RData")
mmetsp_df$Accession <- gsub(".cds.filtered","",mmetsp_df$Assembly)

spp_id <- read.csv("mmetsp_spp.csv")
spp_id$Species.Name <- spp_id$Species
spp_id$Species.Name[!grepl("sp[.]",spp_id$Species.Name)] <-
  spp_id$Species.Name[!grepl("sp[.]",spp_id$Species.Name)] %>%
  strsplit(.,split=" ") %>%
  lapply(.,"[",1:2) %>%
  lapply(.,paste,collapse=" ") %>%
  unlist() %>%
  gsub(pattern=",",replace="")
spp_id <- spp_id %>% subset(Accession %in% mmetsp_df$Accession)
exclude <- readLines("eukprey_list.txt")

mmetsp_growth <- merge.easy(mmetsp_df,spp_id,key="Accession")

mmetsp_growth_df <- merge.easy(mmetsp_growth,growth_rates_dedup_nooutliers,key="Species.Name") %>%
  subset(!is.na(OGT)) %>%
  subset(nHE>10) %>%
  subset(!Species.Name %in% exclude)

#save data
setwd("~/eeggo/Data")
growth_df_mmetsp <- mmetsp_growth_df %>% subset(select=c("Species.Name","Doubling.Time","OGT","CUBHE"))
save(growth_df_mmetsp,file="growth_MMETSP.RData")

write.table(mmetsp_growth_df%>%subset(select=c(Assembly,OGT)),
            file="assemblies_temps_mmetsp.txt",
            sep="\t",
            row.names=F,
            col.names=F,
            quote=F)

