
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

set.seed(123)


pq <- function(m,q_thresh){
  p <- ggplot(m,aes(x=Tpot,y=TPM,color=TPM>quantile(TPM,q_thresh,na.rm=T))) +
    geom_point() +
    theme_pubclean() +
    scale_color_manual(values=c("black","red")) +
    theme(legend.position = "none") +
    geom_vline(xintercept=mean(m$Tpot[m$TPM>=quantile(m$TPM,q_thresh,na.rm=T)],na.rm=T),
               color="red",lty=2) +
    xlab("Potential Temperature (C)")
  return(p)
}


plotQuantiles <- function(m1,m2,m3,m4,m5,m6,m7,m8,m9,q_thresh){
  p1 <- pq(m1,q_thresh)
  p2 <- pq(m2,q_thresh)
  p3 <- pq(m3,q_thresh)
  p4 <- pq(m4,q_thresh)
  p5 <- pq(m5,q_thresh)
  p6 <- pq(m6,q_thresh)
  p7 <- pq(m7,q_thresh)
  p8 <- pq(m8,q_thresh)
  p9 <- pq(m9,q_thresh)
  p <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=3,ncol=3)
  return(p)
}

setwd("~/eeggo/Data/")
x <- read.csv("TOPAZ_MAG_abundance.csv")

rand_mags <- sample(unique(x$new_mag_name),9)
m1 <- x %>% subset(new_mag_name==rand_mags[1])
m2 <- x %>% subset(new_mag_name==rand_mags[2])
m3 <- x %>% subset(new_mag_name==rand_mags[3])
m4 <- x %>% subset(new_mag_name==rand_mags[4])
m5 <- x %>% subset(new_mag_name==rand_mags[5])
m6 <- x %>% subset(new_mag_name==rand_mags[6])
m7 <- x %>% subset(new_mag_name==rand_mags[7])
m8 <- x %>% subset(new_mag_name==rand_mags[8])
m9 <- x %>% subset(new_mag_name==rand_mags[9])

setwd("~/eeggo/Figures/")
pdf("example_quantiles_q99.pdf",height = 8,width=8)
plotQuantiles(m1,m2,m3,m4,m5,m6,m7,m8,m9,0.99)
dev.off()
