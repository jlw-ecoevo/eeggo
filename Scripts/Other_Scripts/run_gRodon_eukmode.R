library(gRodon)
library(dplyr)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("Need filename, output folder, and temperature", call.=FALSE)
}

path_assembly <- args[1]
temperature <- args[3] %>% as.numeric()
out_folder <- args[2]

doAssembly <- function(path_assembly,out_folder,temperature){
  
  
  print(path_assembly)
  genes <- readDNAStringSet(path_assembly)
  ribo_names <- readLines(gsub("[.]cds",".riboprot",path_assembly))
  highly_expressed <- names(genes) %in% c(ribo_names,paste0(ribo_names,"-mRNA"))
  
  if(sum(highly_expressed)>1){
    
    pg <- predictGrowth(genes, 
                          highly_expressed, 
                          mode = "eukaryote",
                          temperature = temperature)
    pgnt <- predictGrowth(genes,
                          highly_expressed,
                          mode = "eukaryote")
    names(pgnt) <- paste0(names(pgnt),".nt")
    out <- c(list(assembly = path_assembly,
                  nHE = sum(highly_expressed)),
             pg,
             pgnt) 
    
    save(out,file=paste0(out_folder,basename(path_assembly),".gRodon.rda"))
    
    return(out)
    
  } else {
    return(NULL)
  }
}


tryDoAssembly <- function(path_assembly,out_folder,temperature){
  try(doAssembly(path_assembly,out_folder,temperature))
}

tryDoAssembly(path_assembly,out_folder,temperature)
