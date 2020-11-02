
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
      labs <- fastafn$tip.label
    }else{ 
      d <- read.dna(fastafn, format = "fasta")
      labs <- dimnames(d)[[1]]
    }
    spltlabs <- sapply(labs,function(x)strsplit(x,'\\|')[[1]])
    epi <- sum(!grepl('^EPI_ISL_[0-9]{6}$',spltlabs[2,]))
    dt <- sum(!str_detect(spltlabs[3,],'^([0-9]{4}[:punct:][0-9]{2}[:punct:][0-9]{2})|([0-9]{2}[:punct:][0-9]{2}[:punct:][0-9]{4})$'))
    if(epi | dt){
      stop('Sequence names must be pipe-separated strings with GISAID ID (e.g. EPI_ISL_420899) in position 2 and date (e.g. 2020-03-11) in position 3.' )
    }
    
    if (inherits(fastafn, "phylo")) {
      tr <- fastafn
    }else{ 
      tr = .mltr(fastafn)
    }
    saveRDS(tr,'pretree.Rds')
    max_length <- decimal_date(today())-decimal_date(as.Date("2019-10-01"))
    
    treeok <- F
    while(!treeok){
      
      # Note requires IQTREE, ncpu > 1 will not work on windows
      tds = make_starting_trees(tr, treeoutfn = nwk_file, ntres = n_startingtrees, ncpu = 6)
      
      lnths <- tds[[1]]$edge.length[ match( 1:Ntip(tds[[1]]) , tds[[1]]$edge[,2] ) ]
      if(sum(lnths>max_length)==0){
        treeok <- T
      }else{
        tr <- ape::drop.tip(tr,tr$tip.label[lnths>max_length])
      }
    }
      
    # save trees for later
    saveRDS(tds,tree_file)
    rm(tr)
    rm(tds)
  }
  return()
  
}


run_skygrowth <- function(region){
  # library(devtools)
  # install_github('https://github.com/emvolz/treedater')
  
  path = paste0(region,"/trees.Rds") # directory with tree file
  mdfn = 'datasets/gisaid.tsv' # path to gisaid metadata 
  
  use_gtds = TRUE 
  minsize = 50
  maxsize = 1000 # will downsample to this value 
  
  nthreads = 1
  ncpu = 16
  
  ntrees = 20
  SGSTART <- as.Date('2020-01-15' )
  MHSTEPS <- 5e5
  RES = 60
  CLOCKRATE <- 0.0008254
  TAU0 <- 1e-5
  
  tres = readRDS( path )[[1]]
  
  md = read.csv(mdfn, sep='\t', stringsAs=FALSE, header=TRUE)
  md$sequence_name <- md$gisaid_epi_isl  # sapply(md$seqName,function(x)strsplit(as.character(x),'\\|')[[1]][2])
  tips <- sapply(tres$tip,function(y)strsplit(as.character(y),'\\|')[[1]][2])
  
  md <- md[md$sequence_name%in%tips,]
  md$sample_time <- decimal_date( as.Date( md$date ))
  
  tr <- tres
  nm <- str_extract(tr$tip,'EPI_ISL_[0-9]{6}') # sapply(tr$tip,function(x)strsplit(as.character(x),'\\|')[[1]][2])
  tr2sts <- setNames( md[ match(nm,md$sequence_name) , 'sample_time'] , tr$tip.label )
  
  #fix seq id, remove missing sts
  sts = tr2sts #[[k]]
  tr = tres #[[k]]
  todrop = names( sts ) [ is.na( sts ) ]
  if ( length( todrop ) > (Ntip(tr)-3))
    todrop <- todrop[1:(Ntip(tr)-3)]
  tres <- drop.tip(tr, todrop  )
  
  .cutie_treedater <- function (tr, ntres = 10, threads = 1, ncpu = 5, mrl = c(CLOCKRATE, CLOCKRATE+1E-4)) 
  {
    
    sts <- sapply(strsplit(tr$tip.label, "\\|"), function(x) {
      as.numeric(tail(x, 2)[1])
    })
    names(sts) <- tr$tip.label
    trpl <- tr
    tr <- di2multi(tr, tol = 1e-05)
    tres <- lapply(1:ntres, function(i) {
      tr = unroot(multi2di(tr))
      tr$edge.length <- pmax(1/29000/5, tr$edge.length)
      tr
    })
    tds <- parallel::mclapply(tres, function(tr) {
      dater(unroot(tr)
            , tr$sts  # sts[tr$tip.label]
            , s = 29000
            , meanRateLimits = mrl
            , ncpu = ncpu
            , searchRoot = ncpu  - 1, omega0 = 0.0012, 
            numStartConditions = 0
      )
    }, mc.cores = threads)
    tds
  }
  
  .lineage_skygrowth <- function(tr, sts)
  {
    # labels 
    tr3 = tr
    tr3$tip.label <- paste(sep='|', tr3$tip.label, sts[tr3$tip.label], 'Il' )
    
    # time trees 
    st0 = Sys.time()
    tds = .cutie_treedater(tr3, ntres = ntrees, threads = nthreads, ncpu = ncpu) 
    st1 <- Sys.time()
    print( paste( 'treedater', st1 - st0 ))
    
    
    # smooth node times 
    if ( use_gtds  & length(tds)>1 ){
      a = capture.output({
        gtds = parallel::mclapply( tds, function(td) gibbs_jitter( td, returnTrees=2 ) , mc.cores = ncpu*nthreads )
      })
      gtds1 <- do.call( c, gtds )
      sg0 = skygrowth1( gtds1, tau0 = TAU0, res = RES, ncpu = ncpu ,  tstart = decimal_date(SGSTART), mhsteps = MHSTEPS)
    } else{
      sg0 = skygrowth1( tds, tau0 = TAU0, res = RES, ncpu = ncpu ,  tstart = decimal_date(SGSTART), mhsteps = MHSTEPS)
      gtds1 = tds
    }
    
    # skygrowth 
    sg0$GRmat  <- sg0$GRmat[ , sample(1:ncol(sg0$GRmat), size = 50, replace=FALSE) ] # shrink it
    sg0$Rmat <- NULL 
    sg0$Ntip = Ntip( tr3 ) 
    
    list( sg = sg0, gtds = gtds1, tds = tds  )
  }
  
  print( Sys.time() )
  #~ st0 = system.time( { sg0 = .lineage_skygrowth(tres[[1]])  } ) #71 sec 
  
  if(!(dir.exists( paste0(region,'/skygrowth3') ))) dir.create( paste0(region,'/skygrowth3') )
  
    tr = tres
    sts <- tr2sts 
    sg_filename <- paste0(region,'/skygrowth3/skygrowth3-sg.rds' )
    if(!file.exists(sg_filename)){
      # downsample as needed 
      if ( Ntip( tr ) > maxsize ){
        m = Ntip(tr) - maxsize 
        drop = sample( tr$tip.label, size = m, replace=FALSE)
        tr = drop.tip( tr,  drop )
        sts = sts[ tr$tip.label ]
      }
      st0 = Sys.time() 
      o = .lineage_skygrowth(tr, sts )
      print(o)
      sg = o$sg
      gtds = o$gtds
      tds = o$tds
      st1 = Sys.time() 
      print( paste( 'skygrowth', st1 - st0 ))
      
      saveRDS(gtds, file = paste0(region,'/skygrowth3/skygrowth3-gtds.rds' ) )
      saveRDS(tds, file = paste0(region,'/skygrowth3/skygrowth3-tds.rds' ) )
      saveRDS(sg, file = sg_filename )
      print( Sys.time() )
    }else{
      sg <- readRDS(sg_filename)
    }
    sgs <-   sg
  
  saveRDS(sgs, file=paste0(region,'/skygrowth3/skygrowth3-sgs.rds' ) )
  
}
  
  
