# load libraries
require(sarscov2)
require(ape)
require(phangorn)
require(ggplot2)
require(ggtree)
require(lubridate)
require(limSolve)
library(skygrowth)
library( treedater )
require(stringr)

source('prep_functions.R')


args = commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)==0){
  region <- 'Peru'
}else{
  region <- args[1]
}



run_skygrowth(region)

