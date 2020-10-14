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

source('alignment_to_tree.R')


## A simple name for the region, state, or country:
args = commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)==0){
  region <- 'Peru'
}else{
  region <- args[1]
}


if(!(dir.exists( region ))) dir.create( region )
setwd(region)

get_region_tree()
