
library(reshape2)
library(heatmap)
library(dplyr)

setwd("~/eukgrowth/Data")
x <- read.delim("busco_presences.txt",sep="\t",head=F)
x$presence <- 1

xw <- reshape(x,direction = "wide",idvar="V1",timevar = "V2")
xw[is.na(xw)] <- 0
heatmap(as.matrix(xw[,-1]),scale="column")

cutoff_org <- 255/2
plot(density(rowSums(xw[,-1])))
abline(v=cutoff_org)


cutoff <- .8*nrow(xw[rowSums(xw[,-1])>cutoff_org,-1])
hist(colSums(xw[rowSums(xw[,-1])>cutoff_org,-1]),breaks=50)
abline(v=cutoff)


names(xw[rowSums(xw[,-1])>cutoff_org,-1])[colSums(xw[rowSums(xw[,-1])>cutoff_org,-1]) > cutoff]
xw[rowSums(xw[,-1])>cutoff_org,1] %>% 
  gsub(pattern = "_.*", replace = "") %>%
  gsub(pattern="MMETSP.*",replace="MMETSP") %>%
  table() 


write(xw[rowSums(xw[,-1])>cutoff_org,1],
      file="organisms_for_tree_from_busco.txt")

write(names(xw[rowSums(xw[,-1])>cutoff_org,-1])[colSums(xw[rowSums(xw[,-1])>cutoff_org,-1]) > cutoff] %>%
        gsub(pattern="presence.",replace=""),
      file="orthologs_for_tree_from_busco.txt")
