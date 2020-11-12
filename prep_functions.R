
get_alignment_function <- function(region,n_region,n_reservoir,missing_threshold){
  
  ## destination path for extracted alignment
  fastafn <- paste0("algn3.fasta")
  
  ## Load gisaid metadata from www.gisaid.org
  md <- read.csv( '../datasets/gisaid.tsv',sep='\t', stringsAs=FALSE ) 
  if(!region%in%c(unique(md$country),unique(md$location),unique(md$division))) 
    stop(paste0('\n"',region,'" is not among cities, counties, countries, regions and states.\n\n'))
  seqnm <- paste0('hCoV-19/',md$strain)
  md$seqName <- sapply(1:length(seqnm),function(x)
    gsub(' ','',paste0(c(seqnm[x],md$gisaid_epi_isl[x],md$date[x],md$region[x]),collapse='|')))
  md$sampleDate <- md$date
  
  ## path to the GISAID alignment
  gisaid_file <- '../datasets/gisaid.fasta'
  
  # Set dedup = TRUE in region sampler and exog sampler to remove duplicate sequences
  regiontips = region_sampler1( md, n = n_region  , inclusion_rules = list( c('country', paste0('^',region,'$')),
                                                                            c('location', paste0('^',region,'$')),
                                                                            c('division', paste0('^',region,'$'))), dedup = FALSE, time_stratify=F)
  
  ## A path to the file containing small genetic distance pairs
  #md2 <- read.csv('../datasets/tn93.txt',sep=',',stringsAsFactors = F)
  # distances from pairsnp
  nms <- readRDS('../datasets/names.Rds')
  dists <- readRDS('../datasets/snp_dist.Rds')
  md2 <- data.frame(ID1=as.character(nms[dists[,1]]),ID2=as.character(nms[dists[,2]]),Distance=dists[,3])
  
  # Sample a set of closely related external sequences, 
  exogtips = exog_sampler2( md, n=n_reservoir, smallGDpairs=md2, region_sample=regiontips, exclusion_rules = list( c('division', paste0('^',region,'$')),
                                                                                                                   c('location', paste0('^',region,'$')),
                                                                                                                   c('country', paste0('^',region,'$')) ), dedup= FALSE )
  
  rm(md2)
  
  prep_tip_labels_seijr <- function (algnfn, outfn, regiontips, exogtips, metadata) {
    md = metadata
    # read in global alignment
    d = read.dna(algnfn, "fasta")
    rownames(d) <- sapply(rownames(d),function(x)gsub(' ','',x))
    # extract tips
    s = intersect(c(regiontips, exogtips), rownames(d))
    .md <- md[match(s, md$seqName), ]
    sts <- setNames(lubridate::decimal_date(lubridate::ymd(as.character(.md$sampleDate))), 
                    .md$seqName)
    dd = d[s, ]
    # omit sequences with less than missing_threshold coverage
    nonmissing <- sapply(1:nrow(dd),function(x)sum(base.freq(dd[x,],all=T)[1:4])>missing_threshold)
    dd <- dd[nonmissing,]
    rm(d)
    # label names
    nms = rownames(dd)
    demes <- setNames(rep("exog", length(nms)), nms)
    demes[nms %in% regiontips2] <- "Il"
    nms = paste(sep = "|", nms, sts[nms], paste0("_", demes[nms]))
    rownames(dd) <- nms
    # write alignment out
    write.dna(dd, file = outfn, format = "fasta")
    invisible(dd)
  }
  
  prep_tip_labels_seijr(algnfn=gisaid_file, outfn=fastafn, regiontips=regiontips, exogtips=exogtips, metadata=md  )
  
  rm(md)
}

make_trees_function <- function(){
  
  tree_file <- paste0("trees.Rds")
  nwk_file <- paste0("startTrees.nwk")
  fastafn <- paste0("algn3.fasta")
  end_files <- c(tree_file,nwk_file)
  if(sum(file.exists(end_files))==length(end_files)) return()
  
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
}

run_skygrowth <- function(region){
  # library(devtools)
  # install_github('https://github.com/emvolz/treedater')
  
  path = paste0(region,"/trees.Rds") # directory with tree file
  mdfn = 'datasets/gisaid.tsv' # path to gisaid metadata 
  
  use_gtds = TRUE 
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
  
  
process_for_shiny_function <- function(region){
  
  lineages <- list()
  
  files <- list.files('./skygrowth3/')
  gtdsfiles <- files[sapply(files,function(x)grepl('\\-gtds',x))]
  sgfiles <- files[sapply(files,function(x)grepl('\\-sg.',x))]
  ids <- sapply(gtdsfiles,function(txt)gsub('.rds','',gsub('skygrowth3-gtds','',txt)))
  
  ## store information
  trees <- list()
  sequences <- list()
  orig_trees <- readRDS(paste0('trees.Rds'))
  i <- which.max(sapply(orig_trees,function(y)y[[7]]))
  
  md <- read.csv( '../datasets/gisaid.tsv',sep='\t', stringsAs=FALSE ) 
  locs <- c('country','division','location')
  loc <- locs[which(sapply(locs,function(x)region%in%md[[x]]) )[1]]
  md$Location <- md[[loc]]
  md$Date <- md$date
  md$Sequence <- paste0('hCoV-19/',md$strain) # sapply(md$seqName,function(y)strsplit(y,'\\|')[[1]][1])
  
  lineages <- readRDS(paste0('skygrowth3/',sgfiles[i]))
  x <- readRDS(paste0('skygrowth3/',gtdsfiles[i]))
  maxind <- which.max(sapply(x,function(y)y[[7]]))
  trees <- x[[maxind]]
  trees$tip.label <- unname(trees$tip.label)
  trees$intree$tip.label <- unname(trees$intree$tip.label)
  trees$lnd.mean.rate.prior <- NULL
  sequences <- as.data.frame(cbind(sapply(x[[1]]$tip.label,function(y)strsplit(y,'\\|')[[1]][1]),ids[i])) 
  colnames(sequences) <- c('Sequence','Lineage')
  sequences <- left_join(sequences,md[,colnames(md)%in%c('Sequence','Location','Date')],by='Sequence')
  rownames(sequences) <- NULL
  trees$data <- as.data.frame(t(sapply(x[[1]]$tip.label,function(y)strsplit(gsub('hCoV-19/','',y),'\\|')[[1]][1:3])),stringsAsFactors=F) # metadata[match(sequences[[i]]$Sequence,metadata$Sequence),] 
  rownames(trees$data) <- NULL
  trees$data$Location <- sapply(trees$data[,1],function(x)strsplit(x,'/')[[1]][1])
  colnames(trees$data) <- c('Sequence name','GISAID','Date','Location')
  
  ## store items
  parms <- list()
  parms$sequences <- sequences#do.call(rbind,sequences)
  parms$tree <- trees
  parms$lineages <- lineages
  
  
  pldf <- compute_timports( orig_trees )
  parms$imports <- pldf[[1]]
  
  saveRDS(parms,paste0(region,'.Rds'))
  saveRDS(parms,paste0('../input_files/',region,'.Rds'))
  
}