# load packages
library(dplyr)
library(sarscov2)

## A simple name for the region, state, or country:
args = commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)==0){
  region <- 'Munich'; population <- 1500000
  region <- 'Sweden'; population <- 10230000
  region <- 'Switzerland'; population <- 8570000
}else{
  region <- args[1]
  population <- as.numeric(args[2])
}


setwd(region)

lineages <- list()

files <- list.files('./skygrowth3/')
gtdsfiles <- files[sapply(files,function(x)grepl('\\-gtds',x))]
sgfiles <- files[sapply(files,function(x)grepl('\\-sg.',x))]
ids <- sapply(gtdsfiles,function(txt)gsub('.rds','',gsub('skygrowth3-gtds','',txt)))

## store information
trees <- list()
sequences <- list()
orig_trees <- readRDS(paste0('trees.Rds'))
i <- which.max(sapply(orig_trees,function(y)y[[7]]))

md <- read.csv( '../datasets/gisaid.tsv',sep='\t', stringsAs=FALSE ) 
locs <- c('country','division','location')
loc <- locs[which(sapply(locs,function(x)region%in%md[[x]]) )[1]]
md$Location <- md[[loc]]
md$Date <- md$date
md$Sequence <- md$strain # sapply(md$seqName,function(y)strsplit(y,'\\|')[[1]][1])

lineages <- readRDS(paste0('skygrowth3/',sgfiles[i]))
x <- readRDS(paste0('skygrowth3/',gtdsfiles[i]))
maxind <- which.max(sapply(x,function(y)y[[7]]))
trees <- x[[maxind]]
sequences <- as.data.frame(cbind(sapply(x[[1]]$tip.label,function(y)strsplit(y,'\\|')[[1]][1]),ids[i])) 
colnames(sequences) <- c('Sequence','Lineage')
sequences <- left_join(sequences,md[,colnames(md)%in%c('Sequence','Location','Date')],by='Sequence')
rownames(sequences) <- NULL
trees$data <- as.data.frame(t(sapply(x[[1]]$tip.label,function(y)strsplit(gsub('hCoV-19/','',y),'\\|')[[1]][1:3])),stringsAsFactors=F) # metadata[match(sequences[[i]]$Sequence,metadata$Sequence),] 
rownames(trees$data) <- NULL
trees$data$Location <- sapply(trees$data[,1],function(x)strsplit(x,'/')[[1]][1])
colnames(trees$data) <- c('Sequence name','GISAID','Date','Location')

## store items
parms <- list()
parms$sequences <- sequences#do.call(rbind,sequences)
parms$tree <- trees
parms$lineages <- lineages


pldf <- compute_timports( orig_trees )
parms$imports <- pldf

saveRDS(parms,paste0(region,'.Rds'))
saveRDS(parms,paste0('../input_files/',region,'.Rds'))

