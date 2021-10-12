# JLW - 2021
# Plot BIOGEOTRACES metagenome predictions from gRodon v2.0.0 vs coverage of
# eukaryotic contigs

# Load Packages + Helper Functions ---------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,
                                    df2,
                                    all.x=TRUE,
                                    by=.EACHI,
                                    allow.cartesian=TRUE)),
                       stringsAsFactors=F))
}

# Load BIOGEOTRACES growth predictions -----------------------------------------

setwd("~/eeggo/Data")
load("CodonStatistics_BIOGEOTRACES.RData")
bio_df <- bio_df %>% subset(nHE>=10)
bio_df$ID <- bio_df$Assembly %>%
  gsub(pattern="metaeuk",replace="") %>%
  gsub(pattern="_.*",replace="")

# add BIOGEOTRACES coverage information ----------------------------------------

euk_cov <- read.delim("euk_cov.txt",head=F) %>% 
  mutate(V1=gsub("_.*","",V1))
names(euk_cov) <- c("ID","EukReads")

total_cov <- read.delim("biogeotraces_clean_fwd_read_counts.txt",head=F) %>%
  mutate(TotalReads=2*V2)
names(total_cov)[1:2] <- c("IDsrr","TotalFwd")

name_convert <- read.csv("GEOTRACES_metadata.csv") %>% 
  subset(select=c(NCBI.SRA.Accession,
                  NCBI.SRA.Accession.for.assembled.metagenome.contigs))
names(name_convert)[1:2] <- c("IDsrr","ID")

total_cov <- merge.easy(total_cov,name_convert,key="IDsrr")

coverage <- merge.easy(euk_cov,total_cov,key="ID")
coverage$EukRel <- coverage$EukReads/coverage$TotalReads


# Merge and plot ---------------------------------------------------------------

bio_cov <- merge.easy(bio_df,coverage,key="ID")

setwd("~/eeggo/Figures")
pdf("BIOGEOTRACES_coverage.pdf",width=7,height=5)
ggplot(bio_cov,aes(x=EukRel,y=d.eukaryotes)) + 
  scale_x_log10() +
  geom_point(size=3) +
  geom_smooth(color="black",fill="black",method="lm") +
  theme_pubclean() +
  xlab("Relative Abundance of Eukaryotes") +
  ylab("Min. Doubling Time Eukaryotes")
dev.off()

cor.test(log(bio_cov$d.eukaryotes),log(bio_cov$EukRel))
