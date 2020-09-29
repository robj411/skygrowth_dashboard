
get_region_tree <- function(region, population=1100000,rerun=F){
  
  
  fasta_file <- paste0("fasta.Rds")
  tree_file <- paste0("trees.Rds")
  nwk_file <- paste0("startTrees.nwk")
  fastafn <- paste0("algn3.fasta")
  end_files <- c(fasta_file,tree_file,nwk_file)
  if(sum(file.exists(end_files))==length(end_files)&rerun==F) return()
  
  # Load the edited gisaid metadata 
  md <- read.csv( '../datasets/gisaid_meta.csv', stringsAs=FALSE ) 
  if(!region%in%c(unique(md$Country),unique(md$RegionOrState),unique(md$CityOrCounty))) 
    stop(paste0('\n"',region,'" is not among cities, counties, countries, regions and states.\n\n'))
  
  # define some parameters 
  ## When do internal SEIR dynamics initiate: 
  startTime = 2020.15
  ## What is the population size of the region of interest?
  popSize = population
  ## How many sequences to include within the region?
  n_region = 1000
  ## How many to include from the international reservoir?
  n_reservoir = 50 # (actual number to be included will be greater since it also includes close distance matches)
  ## How many starting trees? A BEAST run will be carried out for each 
  n_startingtrees = 1 
  
  
  ################################################################################################
  #   We strongly recommend checking the tree and alignment produced by this step and removing   #
  #                          outlying sequences caused by misalignment                           #
  ################################################################################################
  
  # Make the starting trees
  if(!file.exists(nwk_file)|!file.exists(tree_file)|!file.exists(fasta_file)){

     if (inherits(fastafn, "phylo")) {
       tr <- fastafn
     }else{ tr = .mltr(fastafn)}
     saveRDS(tr,fasta_file)
     

    # Note requires IQTREE, ncpu > 1 will not work on windows
    tds = make_starting_trees(tr, treeoutfn = nwk_file, ntres = n_startingtrees, ncpu = 6)
    # save trees for later
    saveRDS(tds,tree_file)
    rm(tr)
    rm(tds)
  }
  return()
  
}
  
  
  
