# load libraries
library(dplyr)
library(sarscov2)

source('prep_functions.R')

args = commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)==0){
  region <- 'Peru'
}else{
  region <- args[1]
}


setwd(region)

process_for_shiny_function(region)
