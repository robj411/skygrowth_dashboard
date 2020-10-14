
get_region_tree <- function(rerun=F){
  
  
  tree_file <- paste0("trees.Rds")
  nwk_file <- paste0("startTrees.nwk")
  fastafn <- paste0("algn3.fasta")
  end_files <- c(tree_file,nwk_file)
  if(sum(file.exists(end_files))==length(end_files)&rerun==F) return()
  
  n_startingtrees = 1
  
  # Make the starting trees
  if(!file.exists(nwk_file)|!file.exists(tree_file)){
    if (inherits(fastafn, "phylo")) {
      tr <- fastafn
    }else{ 
      tr = .mltr(fastafn)
    }
    
    # Note requires IQTREE, ncpu > 1 will not work on windows
    tds = make_starting_trees(tr, treeoutfn = nwk_file, ntres = n_startingtrees, ncpu = 6)
    # save trees for later
    saveRDS(tds,tree_file)
    rm(tr)
    rm(tds)
  }
  return()
  
}
  
  
  
