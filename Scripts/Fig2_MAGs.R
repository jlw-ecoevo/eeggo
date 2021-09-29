# JLW - 2021
# Plot eukaryotic MAG predictions from gRodon v2.0.0 

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(ggridges)
library(tidyverse)

# Helper functions -------------------------------------------------------------

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

# Load Predictions -------------------------------------------------------------


setwd("~/eeggo/Data")
load("CodonStatistics_euk_eukmode.RData")
mmetsp_df <- assembly_df %>% subset(nHE>10)

load("CodonStatistics_alexander.RData")
alexander_df <- assembly_df %>% subset(nHE>10)
alexander_tax <- read.csv("TableS02_EukaryoticMAG.csv")
names(alexander_tax)[1] <- "ID"
alexander_met <- read.csv("TableS09_TOPAZ_HetScore.csv")
names(alexander_met)[1] <- "ID"
alexander_df$ID <- gsub(".cds","",alexander_df$Assembly)
alexander_df <- merge.easy(alexander_df,alexander_tax,key="ID") %>% 
  subset(groups!="Metazoa")
alexander_df <- merge.easy(alexander_df,alexander_met,key="ID")

load("CodonStatistics_delmont.RData")
delmont_df <- assembly_df %>% subset(nHE>10)
delmont_tax <- read.csv("Delmont_taxonomy.csv")
delmont_tax_eukulele <- read.table("delmont_groups.tbl")
names(delmont_tax_eukulele) <- c("ID","groups","percent")
names(delmont_tax)[1] <- "ID"
delmont_df$ID <- gsub(".cds","",delmont_df$Assembly)
delmont_df <- merge.easy(delmont_df,delmont_tax,key="ID") %>% 
  merge.easy(.,delmont_tax_eukulele,key="ID") %>%
  subset(Best_taxonomy_KINGDON!="Animalia")

tax_a <- alexander_df %>% subset(select=c(ID,d,groups))
tax_d <- delmont_df %>% subset(select=c(ID,d,groups))
tax_all <- rbind(tax_a,tax_d)

load("growth_MMETSP.RData")
names(growth_df_mmetsp)[2] <- "d"


# MAG vs Isolate Distribution --------------------------------------------------


pcomp1<- ggplot(NULL,aes(x=d)) + 
  geom_density(data=mmetsp_df,aes(fill=paste0("MMETSP (n=",nrow(mmetsp_df),")")),
               alpha=0.75,color="black") + 
  geom_density(data=tax_all,aes(fill=paste0("MAGs (n=",nrow(tax_all),")")),
               alpha=0.75,color="black") + 
  scale_x_log10(limits=c(.1,100)) +
  theme_pubclean() +
  scale_fill_manual(values=c("black","lightgray")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.key.size = unit(1, "lines")) + 
  geom_vline(xintercept = 10, lty = 2, color = "red") + 
  xlab("Predicted Minimal Doubling Time (Hours)") 
pcomp1

x <- data.frame(Fast=rep(c(F,T),3),
                Data=c("MMETSP","MMETSP",
                       "MMETSP Training Data","MMETSP Training Data",
                       "MAGs","MAGs"),
                Value=c(table(mmetsp_df$d<10)/nrow(mmetsp_df),
                        table(growth_df_mmetsp$d<10)/nrow(growth_df_mmetsp),
                        table(factor(tax_all$d<10,levels=c(FALSE,TRUE)))/nrow(tax_all)))
pcomp2 <- ggplot(x%>%subset(Fast==T),aes(x=Data,y=Value,fill=Data)) + 
  geom_bar(position="stack", stat="identity",alpha=0.75, color = "black") +
  theme_pubclean() + xlab("") + 
  ylab(expression("Proportion Very Fast Growers (<10hr)")) +
  # theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position = "none") +
  # scale_fill_brewer(palette = "Dark2",direction = -1)+
  theme(legend.position = "none") +
  # scale_fill_manual(values = brewer.pal(4,"Dark2")) + 
  scale_fill_manual(values=c("black","lightgray","white")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1,size=6),
        axis.title.y = element_text(size=9))
pcomp2

# MAG vs Isolate t-test
t.test(log10(tax_all$d),log10(mmetsp_df$d))


# Taxonomic Breakdown ----------------------------------------------------------


keep_grp <- names(table(tax_all$groups))[table(tax_all$groups)>2]
tax <- tax_all %>% subset(groups %in% keep_grp) 
y <- tax %>% group_by(groups) %>% summarize(x=n())
tax2 <- merge.easy(tax,y,key="groups")
tax2$label <- paste0(tax2$groups," (n=",tax2$x,")")

x <- tax2 %>% group_by(label) %>% summarize(x=median(d))
levels_vec <- x$label[order(x$x)]
tax2 <- tax2 %>%
  mutate(label = fct_relevel(label, levels =levels_vec))


ptax_all <- ggplot(tax2,aes(x=d, y=label, fill = factor(stat(quantile)))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    scale=0.9,
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.5
  ) +
  scale_fill_manual(values=c("white","lightgray","darkgray","black")) +
  labs(fill="Quartile") + 
  theme_pubclean() + 
  theme(legend.position="none") + 
  xlab("Minimum Doubling Time (Hours)") + 
  ylab("")
ptax_all


# Functional Breakdown ---------------------------------------------------------

adf <- alexander_df %>% subset(groups %in% keep_grp) 
keep_grp <- 
  names(table(alexander_df$PredictedTrophicMode))[table(alexander_df$PredictedTrophicMode)>2]
adf <- adf %>% subset(PredictedTrophicMode %in% keep_grp) 
y <- adf  %>% group_by(PredictedTrophicMode) %>% summarize(x=n())
adf3 <- merge.easy(adf,y,key="PredictedTrophicMode")
adf3$tlabel <- paste0(adf3$PredictedTrophicMode," (n=",adf3$x,")")

x <- adf3 %>% group_by(tlabel) %>% summarize(x=median(d))
levels_vec <- x$tlabel[order(x$x)]
adf4 <- adf3 %>%
  mutate(tlabel = fct_relevel(tlabel, levels =levels_vec))

ptroph <- ggplot(adf4,aes(x=d, y=tlabel, fill = factor(stat(quantile)))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    scale=0.9,
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.5
  ) +
  scale_fill_manual(values=c("white","lightgray","darkgray","black")) +
  labs(fill="Quartile") + 
  theme_pubclean() + 
  theme(legend.position="right") + 
  xlab("Minimum Doubling Time (Hours)") + 
  ylab("") #+
  # scale_x_log10()
ptroph


t.test(log10(d)~PredictedTrophicMode,data=alexander_df)


# put figure together ----------------------------------------------------------


setwd("~/eeggo/Figures")
pdf("mag_panels.pdf",width=12,height=8)
ggarrange(ggarrange(ggarrange(pcomp1,
                              pcomp2,
                              ncol=2,
                              widths=c(5,1),
                              labels=c("(a)","(b)"),
                              vjust=1,
                              hjust=0),
                    ptroph,
                    nrow=2,
                    heights=c(3,1),
                    labels=c("","(c)")),
          ptax_all,
          ncol=2,
          labels=c("","(d)"),
          vjust=1,
          widths=c(3,2))
dev.off()

