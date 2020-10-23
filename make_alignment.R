# load packages
require(sarscov2)
require(ape)
require(phangorn)
require(ggplot2)
require(ggtree)
require(lubridate)
require(limSolve)
library( treedater )


## A simple name for the region, state, or country:
args = commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)==0){
  region <- 'Munich'
  region <- 'Sweden'
  region <- 'Switzerland'
}else{
  region <- args[1]
}


if(!(dir.exists( region ))) dir.create( region )
setwd(region)
print(region)
  
  
fastafn <- paste0("algn3.fasta")
  
# Load the edited gisaid metadata 
md <- read.csv( '../datasets/gisaid.tsv',sep='\t', stringsAs=FALSE ) 
if(!region%in%c(unique(md$country),unique(md$location),unique(md$division))) 
  stop(paste0('\n"',region,'" is not among cities, counties, countries, regions and states.\n\n'))
seqnm <- paste0('hCoV-19/',md$strain)
md$seqName <- sapply(1:length(seqnm),function(x)paste0(c(seqnm[x],md$gisaid_epi_isl[x],md$date[x],md$region[x]),collapse='|'))
md$sampleDate <- md$date

## path to the GISAID alignment: 
gisaid_file <- '../datasets/gisaid.fasta'
## A path to the file containing small genetic distance pairs
distfn = '../datasets/tn93.txt'
  
# define some parameters 
## How many sequences to include within the region?
n_region = 1000
## How many to include from the international reservoir?
n_reservoir = 50 # (actual number to be included will be greater since it also includes close distance matches)
n_startingtrees = 1
  
  
# Set dedup = TRUE in region sampler and exog sampler to remove duplicate sequences
regiontips = region_sampler1( md, n = n_region  , inclusion_rules = list( c('country', paste0('^',region,'$')),
                                                                              c('location', paste0('^',region,'$')),
                                                                              c('division', paste0('^',region,'$'))), dedup = FALSE, time_stratify=F)
# Sample a set of closely related external sequences, 
md2 <- read.csv(distfn,sep=',',stringsAsFactors = F)
exogtips = exog_sampler2( md, n=n_reservoir, smallGDpairs=md2, region_sample=regiontips, exclusion_rules = list( c('division', paste0('^',region,'$')),
                                                                                                                     c('location', paste0('^',region,'$')),
                                                                                                                     c('country', paste0('^',region,'$')) ), dedup= FALSE )
# You could select based on a different geographic criterion, current fields in metadata are  "Continent" "Country" "RegionOrState" "CityOrCounty"
rm(md2)
#~ Here is a subset of the metadata that corresponds to the sample. Check this carefully- is the sample including everything you want and not including things you dont want?
mdregion = md[ match( regiontips, md$seqName ) , ]
    


prep_tip_labels_seijr <- function (algnfn, outfn, regiontips, exogtips, metadata) {
  md = metadata
  d = read.dna(algnfn, "fasta")
  rownames(d) <- sapply(rownames(d),function(x)gsub(' ','',x))
  s = intersect(c(regiontips, exogtips), rownames(d))
  .md <- md[match(s, md$seqName), ]
  sts <- setNames(lubridate::decimal_date(lubridate::ymd(as.character(.md$sampleDate))), 
                  .md$seqName)
  dd = d[s, ]
  rm(d)
  nms = rownames(dd)
  demes <- setNames(rep("exog", length(nms)), nms)
  demes[nms %in% regiontips] <- "Il"
  nms = paste(sep = "|", nms, sts[nms], paste0("_", demes[nms]))
  rownames(dd) <- nms
  write.dna(dd, file = outfn, format = "fasta")
  invisible(dd)
}

prep_tip_labels_seijr(algnfn=gisaid_file, outfn=fastafn, regiontips=regiontips, exogtips=exogtips, metadata=md  )

rm(md)
  
  
  
