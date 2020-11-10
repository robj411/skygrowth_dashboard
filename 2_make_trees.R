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

source('prep_functions.R')


args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  region <- 'Peru'
}else{
  region <- args[1]
}


if(!(dir.exists( region ))) dir.create( region )
setwd(region)

make_trees_function()

