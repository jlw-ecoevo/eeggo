# JLW - 2021
# Plot performance data for gRodon v2.0.0 predictions on eukaryotes

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(MASS)
library(reshape2)
library(psych)

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

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}

# Load Training Data -----------------------------------------------------------

setwd("~/eeggo/Data")
load("growth_MMETSP.RData")
load("growth_RefSeq.RData")
growth_df$Source <- "GenBank"
growth_df_mmetsp$Source <- "MMESTP"
growth_df <- rbind(growth_df,growth_df_mmetsp)

# Training Data Characteristics ------------------------------------------------

# Check Normality
qqnorm(growth_df$Doubling.Time %>% log10())
qqline(growth_df$Doubling.Time %>% log10())

# Test for outliers
rosnerTest(growth_df$Doubling.Time %>% log10())

setwd("~/eeggo/Figures")
pdf("training_data.pdf",width=8,height=4)
ggplot(growth_df,aes(x=Doubling.Time,fill=Source)) +
  geom_density(alpha=0.5) +
  scale_x_log10() +
  xlab("Minimum Doubling Time (Hours)") +
  # geom_vline(xintercept = 40,lty=2) +
  theme_pubclean() +
  scale_fill_manual(values=c("black","lightgray"))
dev.off()

# RefSeq data faster growth rates than MMETSP data
t.test(log10(Doubling.Time)~Source,data=growth_df)

# Train Models ---------------- ------------------------------------------------

train_df <- growth_df

#Test if CUB correlated w/ growth
cor.test(log(train_df$Doubling.Time),train_df$CUBHE)

#linear model w/ temperature
m_milc <- lm(Doubling.Time~CUBHE+OGT,data=train_df)

#Box-Cox transformation
bc_milc <- MASS::boxcox(Doubling.Time~CUBHE+OGT,data=train_df)
lambda_milc <- bc_milc$x[which.max(bc_milc$y)]
# re-run with transformation
mnew_milc <-
  lm(boxcoxTransform(Doubling.Time, lambda_milc) ~ CUBHE+OGT,data=train_df)

#No Temperature
#transformation
bc_milc_notemp <- MASS::boxcox(Doubling.Time~CUBHE,data=train_df)
lambda_milc_notemp <- bc_milc$x[which.max(bc_milc$y)]
# re-run with transformation
mnew_milc_notemp <-
  lm(boxcoxTransform(Doubling.Time, lambda_milc) ~ CUBHE,data=train_df)

#Only Temperature
#transformation
bc_temp <- MASS::boxcox(Doubling.Time~OGT,data=train_df)
lambda_temp <- bc_temp$x[which.max(bc_temp$y)]
# re-run with transformation
mnew_temp <-
  lm(boxcoxTransform(Doubling.Time, lambda_temp) ~ OGT,data=train_df)

#look at residuals
ggqqplot(m_milc$residuals)
ggqqplot(mnew_milc$residuals)
ggqqplot(mnew_milc_notemp$residuals)
ggqqplot(mnew_temp$residuals)

summary(mnew_milc)
summary(mnew_milc_notemp)
summary(mnew_temp)

#DF of actual and fitted values
train_df$dGR <- boxcoxTransform(mnew_milc$fitted.values,
                                lambda_milc,
                                back_transform = TRUE)
train_df$dGRnt <- boxcoxTransform(mnew_milc_notemp$fitted.values,
                                  lambda_milc_notemp,
                                  back_transform = TRUE)




p1la <- ggplot(train_df,aes(x=Doubling.Time,y=dGR,fill=OGT)) +
  geom_point(alpha=1,size=2,pch=21) +
  # xlim(0,125) +
  # ylim(0,125) +
  scale_x_log10(limits=c(2,250)) +
  scale_y_log10(limits=c(2,250)) +
  theme_pubclean() +
  #geom_smooth(color="darkgrey") +
  xlab("Empirical Minimal Doubling Time (Hours)") +
  ylab("Predicted Minimal Doubling Time (Hours)") +
  geom_abline(slope = 1,intercept = 0,lty=2) +
  geom_vline(xintercept = 40,lty=2,color="black") +
  # ggtitle("Model with CUB and OGT") +
  #geom_smooth(method="gam",color="gray") +
  labs(fill="Optimal Growth\nTemperature (C)") +
  # theme(legend.position = c(0.8,0.8)) +
  # scale_fill_gradient(high="white",low="black")
  scale_fill_gradient(low="white",high="black")
p1la


p3a <- ggplot(train_df,aes(x=Doubling.Time,y=CUBHE,fill=OGT)) +
  geom_point(alpha=1,size=2,pch=21) +
  #scale_x_log10(limits=c(2.5,80)) +
  xlim(2,250) +
  theme_pubclean() +
  #geom_smooth(method="gam",color="gray") +
  geom_smooth(method="loess",fill="black",color="white") +
  # geom_smooth(method="lm",aes(color=Doubling.Time>40),fill="black") +
  xlab("Empirical Minimal Doubling Time (Hours)") +
  ylab("Codon Usage Bias (Ribosomal Proteins)") +
  geom_vline(xintercept = 40,lty=2,color="black") +
  ggtitle("") +
  theme(legend.position = "none") +
  # scale_fill_gradient(high="white",low="black")
  scale_fill_gradient(low="white",high="black")
p3a


p4a <- ggplot() +theme_minimal()

# Load gRodon v1 and growthpred predictions on euks ----------------------------

setwd("~/eeggo/Data")
load("CodonStatistics_euk.RData")
assembly_df$Accession <- assembly_df$Assembly %>%
  gsub(pattern="[.].*",replace="")
gp_rates <- read.delim("growthpred_predicted_euk.tbl",sep="\t",head=F)
gp_rates$V2 <- NULL
names(gp_rates) <- c("Accession","d.growthpred")
gp_rates$Accession <- gp_rates$Accession %>%
  gsub(pattern="[.].*",replace="")
assembly_df <- merge.easy(assembly_df,gp_rates,key="Accession")

spp_id <- read.csv("mmetsp_spp.csv")
spp_id$Species.Name <- spp_id$Species
spp_id$Species.Name[!grepl("sp[.]",spp_id$Species.Name)] <-
  spp_id$Species.Name[!grepl("sp[.]",spp_id$Species.Name)] %>%
  strsplit(.,split=" ") %>%
  lapply(.,"[",1:2) %>%
  lapply(.,paste,collapse=" ") %>%
  unlist() %>%
  gsub(pattern=",",replace="")
spp_id <- spp_id %>%
  subset(Accession %in% assembly_df$Accession) %>%
  subset(select=c(Accession,Species.Name))


x <- read.csv("euk_growth_ogt_Weissman_RefSeq.csv",stringsAsFactors = F)
spp_id_ncbi <- data.frame(Accession=character(),
                           Species.Name=character(),
                           Doubling.Time=numeric(),
                           OGT=numeric())
for(i in 1:nrow(x)){
  spp_id_ncbi <- rbind(spp_id_ncbi,
                        data.frame(Accession=x[i,"Accession"] %>%
                                     strsplit(split=";") %>%
                                     unlist() %>%
                                     gsub(pattern="[.].*",replace=""),
                                   Species.Name=x[i,"Species.Name"]))
}

mmode_df <- merge.easy(assembly_df,
                       rbind(spp_id,spp_id_ncbi),
                       key="Accession") %>%
  subset(select=c(Species.Name,
                  d,
                  d.madin,
                  d.growthpred)) %>%
  group_by(Species.Name) %>%
  summarise_all(mean)


train_df <- merge.easy(train_df,mmode_df,key="Species.Name")

train_melt <- melt(train_df,id=c("Species.Name",
                                 "Doubling.Time",
                                 "OGT","CUBHE",
                                 "Source",
                                 "dGRnt",
                                 "d.madin"))
train_melt$MSE <- 
  (boxcoxTransform(train_melt$Doubling.Time, lambda_milc) - 
     boxcoxTransform(train_melt$value, lambda_milc))^2


p1pve <- ggplot(train_melt,aes(x=Doubling.Time,
                               y=value,
                               fill=variable,
                               shape=variable)) +
  geom_point(pch=21,size=2) +
  theme_pubclean() +
  scale_x_log10(limits=c(1,250)) +
  scale_y_log10(limits=c(0.02,100)) +
  geom_abline(intercept=0,slope=1,lty=2) +
  geom_vline(xintercept=40,lty=2) +
  xlab("Empirical Minimal Doubling Time (Hours)") +
  ylab("Predicted Minimal Doubling Time (Hours)") +
  scale_fill_manual(values=c("black","lightgray","white"),
                    labels=c("gRodon-Eukaryotic","gRodon-Prokaryotic","growthpred")) +
  labs(fill="") +
  theme(legend.position = c(0.1,0.87),
        legend.key=element_blank(),
        legend.background=element_blank())

p2pve <- ggplot(train_melt,aes(y=MSE,group=variable,x=variable)) +
  geom_boxplot() +
  theme_pubclean() +
  scale_x_discrete(labels=c("gRodon-Eukaryotic","gRodon-Prokaryotic","growthpred")) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60,hjust=1))

setwd("~/eeggo/Figures")
pdf("gRodon-me_performance_panels.pdf",width=10,height=10)
ggarrange(ggarrange(p4a,
                    ggarrange(p1la,
                              p3a,
                              nrow=1,
                              labels = c("(a)","(b)"),
                              hjust=0,
                              vjust=-1.75),
                    p4a,
                    nrow=3,
                    heights = c(1,9,1)),
          ggarrange(p1pve,
                    p2pve,
                    ncol=2,
                    widths=c(3,1),
                    labels=c("(c)","(d)"),
                    hjust=0,
                    vjust=-1),
          nrow=2)
dev.off()
