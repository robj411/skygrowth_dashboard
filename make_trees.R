# load libraries
require(sarscov2)
require(ape)
require(phangorn)
require(ggplot2)
require(ggtree)
require(lubridate)
require(limSolve)
require(skygrowth)
require(treedater)
require(dplyr)
require(stringr)

#source('prep_functions.R')


args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  region <- 'Peru'
}else{
  region <- args[1]
}


if(!(dir.exists( region ))) dir.create( region )
setwd(region)



tree_file <- paste0("trees.Rds")
nwk_file <- paste0("startTrees.nwk")
fastafn <- paste0("algn3.fasta")
end_files <- c(tree_file,nwk_file)
if(sum(file.exists(end_files))==length(end_files)) return()

n_startingtrees = 1

# Make the starting trees
if(!file.exists(nwk_file)|!file.exists(tree_file)){
  if (inherits(fastafn, "phylo")) {
    labs <- fastafn$tip.label
  }else{ 
    d <- read.dna(fastafn, format = "fasta")
    labs <- dimnames(d)[[1]]
  }
  spltlabs <- sapply(labs,function(x)strsplit(x,'\\|')[[1]])
  epi <- sum(!grepl('^EPI_ISL_[0-9]{6}$',spltlabs[2,]))
  dt <- sum(!str_detect(spltlabs[3,],'^([0-9]{4}[:punct:][0-9]{2}[:punct:][0-9]{2})|([0-9]{2}[:punct:][0-9]{2}[:punct:][0-9]{4})$'))
  if(epi | dt){
    stop('Sequence names must be pipe-separated strings with GISAID ID (e.g. EPI_ISL_420899) in position 2 and date (e.g. 2020-03-11) in position 3.' )
  }
  
  if (inherits(fastafn, "phylo")) {
    tr <- fastafn
  }else{ 
    tr = .mltr(fastafn)
  }
  saveRDS(tr,'pretree.Rds')
  max_length <- decimal_date(today())-decimal_date(as.Date("2019-10-01"))
  
  treeok <- F
  while(!treeok){
    
    # Note requires IQTREE, ncpu > 1 will not work on windows
    tds = make_starting_trees(tr, treeoutfn = nwk_file, ntres = n_startingtrees, ncpu = 6)
    
    lnths <- tds[[1]]$edge.length[ match( 1:Ntip(tds[[1]]) , tds[[1]]$edge[,2] ) ]
    if(sum(lnths>max_length)==0){
      treeok <- T
    }else{
      tr <- ape::drop.tip(tr,tr$tip.label[lnths>max_length])
    }
  }
  
  # save trees for later
  saveRDS(tds,tree_file)
  rm(tr)
  rm(tds)
}
