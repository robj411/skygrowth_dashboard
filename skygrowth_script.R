
run_skygrowth <- function(region){
  # library(devtools)
  # install_github('https://github.com/emvolz/treedater')
  
  path = paste0(region,"/startTrees.nwk") # directory with tree files
  fsta <- readRDS(paste0(region,"/fasta.Rds"))
  mdfn = 'datasets/gisaid_meta.csv' # path to csv 
  
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
  
  tres = read.tree( path )
  
  md = read.csv(mdfn,  stringsAs=FALSE, header=TRUE)
  md$sequence_name <- sapply(md$seqName,function(x)strsplit(as.character(x),'\\|')[[1]][2])
  tips <- sapply(tres$tip,function(y)strsplit(as.character(y),'\\|')[[1]][2])
  
  md <- md[md$sequence_name%in%tips,]
  md$sample_time <- decimal_date( as.Date( md$date_submitted ))
  
    tr <- tres
    nm <- sapply(tr$tip,function(x)strsplit(as.character(x),'\\|')[[1]][2])
    tr2sts <- setNames( md[ match(nm,md$sequence_name) , 'sample_time'] , tr$tip.label )

  #fix seq id, remove missing sts
    sts = tr2sts #[[k]]
    tr = tres #[[k]]
    todrop = names( sts ) [ is.na( sts ) ]
    if ( length( todrop ) > (Ntip(tr)-3))
      todrop <- todrop[1:(Ntip(tr)-3)]
    
    tres <- drop.tip(tr, todrop  )
  
  last_sample_time <- max(unlist(tr2sts),na.rm=T)
  
  .cutie_treedater <- function (tr, ntres = 10, threads = 1, ncpu = 5, mrl = c(CLOCKRATE, CLOCKRATE+1E-5)) 
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
            , sts[tr$tip.label]
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
    tds = .cutie_treedater(fsta, ntres = ntrees, threads = nthreads, ncpu = ncpu) 
    st1 <- Sys.time()
    print( paste( 'treedater', st1 - st0 ))
    
    
    # smooth node times 
    if ( use_gtds  & length(tds)>1 ){
      a = capture.output({
        gtds = parallel::mclapply( tds, function(td) gibbs_jitter( td, returnTrees=2 )[[2]] , mc.cores = ncpu*nthreads )
      })
      sg0 = skygrowth1( gtds, tau0 = TAU0, res = RES, ncpu = ncpu ,  tstart = decimal_date(SGSTART), mhsteps = MHSTEPS)
    } else{
      sg0 = skygrowth1( tds, tau0 = TAU0, res = RES, ncpu = ncpu ,  tstart = decimal_date(SGSTART), mhsteps = MHSTEPS)
      gtds = tds
    }
    
    # skygrowth 
    sg0$GRmat  <- sg0$GRmat[ , sample(1:ncol(sg0$GRmat), size = 50, replace=FALSE) ] # shrink it
    sg0$Rmat <- NULL 
    sg0$Ntip = Ntip( tr3 ) 
    
    list( sg = sg0, gtds = gtds, tds = tds  )
  }
  
  print( Sys.time() )
  #~ st0 = system.time( { sg0 = .lineage_skygrowth(tres[[1]])  } ) #71 sec 
  
  if(!(dir.exists( paste0(region,'/skygrowth3') ))) dir.create( paste0(region,'/skygrowth3') )
  
    tr = tres
    sts <- tr2sts 
    ln <- c()
    sg_filename <- paste0(region,'/skygrowth3/skygrowth3-sg', ln, '.rds' )
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
      
      saveRDS(gtds, file = paste0(region,'/skygrowth3/skygrowth3-gtds', ln, '.rds' ) )
      saveRDS(tds, file = paste0(region,'/skygrowth3/skygrowth3-tds', ln, '.rds' ) )
      saveRDS(sg, file = sg_filename )
      print( Sys.time() )
    }else{
      sg <- readRDS(sg_filename)
    }
    sgs <-   sg
  
  saveRDS(sgs, file=paste0(region,'/skygrowth3/skygrowth3-sgs', '.rds' ) )
  
}
  
  
