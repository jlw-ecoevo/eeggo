# JLW - 2021
# Supplemental figure phylogenetic patterns in growth rate

# Load Packages + Helper Functions ---------------------------------------------

library(dplyr)
library(ape)
library(ggtree)
library(phytools)
library(rcartocolor)
library(data.table)
library(ggpubr)
library(ggpointdensity)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,
                                    by=.EACHI,allow.cartesian=TRUE)),
                       stringsAsFactors=F))
}

# Load tree and growth data ----------------------------------------------------

setwd("~/eeggo/Data")
tree <- read.tree("busco_tree_pretrim.newick")
delmont_tax <- read.table("delmont_groups.tbl") %>% 
  subset(select=c(V1,V2))
names(delmont_tax) <- c("id","groups")
alexander_tax <- read.csv("TableS02_EukaryoticMAG.csv") %>% 
  subset(select=c(X,groups))
names(alexander_tax)[1] <- "id"
taxonomy <- rbind(delmont_tax,alexander_tax)
names(taxonomy)[2] <- "tax"
tree$tip.label <- tree$tip.label %>% gsub(pattern="[.].*",replace="")

### Growth Data
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

mmetsp_df <- mmetsp_df %>% 
  subset(nHE>10) %>%
  subset(select=c(Assembly,d))
alexander_df <- alexander_df %>% 
  subset(nHE>10) %>%
  subset(select=c(Assembly,d)) 
delmont_df <- delmont_df %>% 
  subset(nHE>10) %>%
  subset(select=c(Assembly,d)) 
growth_df <- rbind(mmetsp_df,alexander_df,delmont_df) %>% unique()
growth_df$ID <- gsub("[.].*","",growth_df$Assembly)
rownames(growth_df) <- growth_df$ID

gdf <- growth_df %>% subset(select=c(ID,d))
tree.mag <- drop.tip(tree,which(!tree$tip.label %in% growth_df$ID))

# Plot Tree --------------------------------------------------------------------

nColor <- 10
p <- ggtree(tree.mag, layout="rectangular") 
setwd("~/eeggo/Figures")
pdf("guidetree.pdf",height=20,width=10)
p %<+% taxonomy + geom_tippoint(aes(color=tax),alpha=1,size=1) + 
  scale_color_manual(values=carto_pal(nColor, "Safe")) +
  geom_tiplab(size=1)
dev.off()

t_off <- 0.06
h_just <- 0.1
pb <- ggtree(tree.mag, layout="rectangular") +
  geom_strip("TOPAZ_NPS1_E023","TARA_ARC_108_MAG_00245",
             barsize=2, color="black",
             label="Ochrophyta",angle=60,hjust=h_just,offset.text = t_off) +
  geom_strip("TOPAZ_SPS1_E038","MMETSP1464",
             barsize=2, color="black",
             label="Haptohpyta",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("TOPAZ_SPD1_E008","MMETSP0941",
             barsize=2, color="black",
             label="Chlorophyta",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("TOPAZ_NAD1_E011","TOPAZ_NAM1_E004",
             barsize=2, color="black",
             label="Fungi",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("TOPAZ_SAD1_E027","MMETSP0108",
             barsize=2, color="black",
             label="Cryptophyta",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("MMETSP1052","MMETSP0040",
             barsize=2, color="black",
             label="Rhizaria",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("TOPAZ_SAM1_E007","MMETSP0413",
             barsize=2, color="black",
             label="Amoebozoa",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("MMETSP1371","MMETSP0469",
             barsize=2, color="black",
             label="Dinophtya",angle=60,hjust=h_just,offset.text = t_off)+
  geom_strip("TOPAZ_NAS1_E015","TARA_ARC_108_MAG_00325",
             barsize=2, color="black",
             label="Dinophtya",angle=60,hjust=h_just,offset.text = t_off)
gdf2 <- gdf %>% 
  subset(select=c(d))
p1 <- gheatmap(pb,
         gdf2,
         low="white",
         high="black",
         colnames=F,
         legend_title= "Min. Doubling Time (h)",
         width = 0.5) +
  theme(legend.position = c(0.1,0.9)) 


# Plot distances (phylogenetic and pairwise difference in growth potential) ----

tree_dist <- cophenetic.phylo(tree.mag) %>% as.data.frame()
tree_dist$ID2 <- rownames(tree_dist) 
dist_df <- tree_dist %>% 
  melt(id.vars="ID2") %>%
  subset(ID2!=variable)
dist_df$variable <- as.character(dist_df$variable)
dist_df <- dist_df %>% 
  subset(variable %in% growth_df$ID) %>%
  subset(ID2 %in% growth_df$ID) %>% 
  subset(value>1e-6)

dist_df$d_diff <- NA
for(i in 1:nrow(dist_df)){
  dist_df$d_diff[i] <- 
    abs(growth_df[growth_df$ID==dist_df$variable[i],"d"] - 
          growth_df[dist_df$ID2[i],"d"])
}

dist_df$d_diff <- unlist(dist_df$d_diff)
summary(lm(log10(d_diff)~log10(value),data=dist_df))
summary(lm(log10(d_diff)~log10(value),data=dist_df))$coef[2,4]

p2 <- ggplot(dist_df,aes(x=value,y=d_diff)) + 
  geom_pointdensity() + 
  scale_color_continuous(low="black",high="lightgray") +
  geom_smooth(method="loess",color="red",se=F,lwd=2) +
  scale_x_log10() +
  theme_pubclean() + 
  geom_hline(yintercept = 6,lty=2,color="blue",lwd=1) +
  geom_hline(yintercept = 1,lty=2,color="blue",lwd=1) +
  geom_vline(xintercept = 0.1,lty=2,color="blue",lwd=1) +
  xlab("Patristic Distance") +
  ylab(expression("Difference in Min. Doubling Times (|" * d[1] - d[2] * "|)")) +
  theme(legend.position = "right",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))

setwd("~/eeggo/Figures")
png("tree_panels.png",height=1000,width=1500)
ggarrange(p1,p2,ncol=2,widths=c(2,3),labels=c("(a)","(b)"))
dev.off()
