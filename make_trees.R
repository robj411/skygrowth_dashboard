# load packages
require(sarscov2)
require(ape)
require(phangorn)
require(ggplot2)
require(ggtree)
require(lubridate)
require(limSolve)
library(skygrowth)
library( treedater )
library(dplyr,quietly = T)

source('alignment_to_tree.R')


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


if(!(dir.exists( region ))) dir.create( region )
setwd(region)

get_region_tree(region,population=population)
