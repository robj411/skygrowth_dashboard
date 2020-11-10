# load packages
require(ape)
library(dplyr)
require(ggplot2)
require(ggtree)
require(limSolve)
require(lubridate)
require(phangorn)
library(sarscov2)
require(skygrowth)
require(stringr)
require(treedater)

source('prep_functions.R')

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  region <- 'Peru'
  ## How many sequences to include within the region?
  n_region = 1000
  ## How many to include from the international reservoir?
  n_reservoir = 50 # (actual number to be included will be greater since it also includes close distance matches)
  missing_threshold <- 0.9
}else{
  region <- args[1]
  ## How many sequences to include within the region?
  n_region = as.numeric(args[2])
  ## How many to include from the international reservoir?
  n_reservoir = as.numeric(args[3]) # (actual number to be included will be greater since it also includes close distance matches)
  missing_threshold <- as.numeric(args[4])
}


if(!(dir.exists( region ))) dir.create( region )
setwd(region)
print(region)

get_alignment_function(region,n_region,n_reservoir,missing_threshold)

make_trees_function()

setwd('..')

run_skygrowth(region)

setwd(region)

process_for_shiny_function(region)

  
  
