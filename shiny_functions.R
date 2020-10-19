## plotting functions ####################################


quick_annotated_treeplot <- function( td , annotation='Location', maxdate = NULL){ #date_decimal( max(tr$sts))
  tr = td 
  len <- length(unique(tr$data[[annotation]]))
  if(annotation=='d614g') {
    cols <- c('navyblue','turquoise','darkorange')
    names(cols) <- c('D','G','X')
    leg.dir <- "vertical"
    shapes <- rep(19,3)
  }else{
    if(len<10){
      cols <- rainbow(len)
      names(cols) <- unique(tr$data[[annotation]])
      shapes <- rep(19,len)
    }else{
      pchs <- c(15:18,25,5,4)
      if(len<60)  pchs <- c(15:18,25,5)
      if(len<50)  pchs <- c(15:18,25)
      if(len<40)  pchs <- c(15:18)
      reps <- ceiling(len/length(pchs))
      cols <- rep(rainbow(reps),times=length(pchs))[1:len]
      names(cols) <- unique(tr$data[[annotation]])
      shapes <- rep(pchs,each=reps)[1:len]
      names(shapes) <- unique(tr$data[[annotation]])
    }
    leg.dir <- "horizontal"
  }
  colScale <- scale_colour_manual(name = annotation,values = cols)
  shapeScale <- scale_shape_manual(name = annotation,values = shapes) 
  class( tr ) = 'phylo'
  maxdate <- max((tr$data$Date)) # date_decimal( as.Date(max(tr$data$Date))) # max(decimal_date(tr$data$Date))
  #tipdeme <-  grepl( tr$tip.label, pat = region ) 
  tipdata <- data.frame( 
    taxa = tr$tip.label, 
    anno =  tr$data[[annotation]],
    text = paste0(tr$data$Country,' / ',tr$data$Date)
  )
  tipdata$size <- .75
  #tipdata$size[ !tipdata$d614g ] <- 0
  #tipdata$d614g[ !tipdata$d614g ] <- NA
  btr = ggtree(tr, mrsd= maxdate, ladderize=TRUE, as.Date=T,aes(text=text))  + theme_tree2() 
  btr <- btr %<+% tipdata 
  metat <- btr$data 
  
  #btr = btr + geom_tippoint( aes(color = anno, pch=anno, size = size, text=text), na.rm=TRUE, show.legend=TRUE, size =1.5) 
  btr = btr + geom_point( data=na.omit(metat),aes(color = anno, pch=anno)) 
  
  btr = btr + colScale + shapeScale +
    theme(legend.position='top', 
          legend.justification='left',
          legend.title=element_text(size=14), 
          legend.text=element_text(size=12),
          #legend.direction=leg.dir,
          axis.text=element_text(size=12)) +
    scale_x_date(date_labels = "%Y/%m/%d") + theme(legend.title=element_blank())
  
  ## add some space under legend otherwise it covers tree
  hgt <- length(tr$Ti)
  height = 3*hgt + 1000 + ceiling(len/7)*20
  btrplotly <- ggplotly(btr,height=height, tooltip = "text") %>%
    config(displayModeBar = F)   %>%
    style(hoverinfo = "none", traces = 1:2) %>%
    layout(legend = list(orientation = "h",xanchor='center',
                         yanchor='bottom',
                         x=0.4,y=1))
  # y=1 at the top and y=0 at the bottom. then the legend takes up a space ceiling(len/7)/hgt
  #leg_height <- ceiling(len/7)/height
  #newy <- 1 + 20*leg_height + (1-btrplotly$x$layout$legend$y) #+ 4*hgt_unit*floor(len/7)
  #btrplotly <- btrplotly %>% 
  #  layout(legend = list(orientation = "h",xanchor='center',
  #                       yanchor='bottom',
  #                       x=0.4,y=newy))
  btrplotly
}


gg_sarscov2skygrowth <- function(x, metric='growth', log_size=F,date_limits = c( as.Date( '2020-03-01'), NA ) ,ci,col,... ){
  #stopifnot( inherits( x, 'sarscov2skygrowth' ))
  y = x[[metric]]
  taxis = as.Date( y$time )
  
  if ( is.na( date_limits[2]) )
    date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
  #qs <- c( .5, .025, .975 )
  
  pldf <- data.frame( Date = taxis , reported=FALSE )
  pldf$out = y$pc50
  pldf$`2.5%` = y$pc2.5
  pldf$`97.5%` = y$pc97.5
  
  pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
  pl = ggplot( pldf) + 
    geom_path( aes(x = Date, y = out ), lwd=0.8,col=col) 
  if(log_size) pl <- pl  +  scale_y_continuous(limits=c(1e-3,NA), trans='log',label=scientific_format(digits=2)) 
  
  if(ci==1) {
    pl <- pl +  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`) , alpha = .1 ,col=col,lwd=0) 
  }
  
  y_lab <- 'Effective growth rate'
  if(metric=='R') y_lab <- 'Effective reproduction number'
  if(metric=='Ne') y_lab <- 'Effective population size'
  
  if(metric%in%c('R','growth')) pl <- pl + geom_hline( aes(yintercept = 1 ), colour = 'red' )
  pl <- pl + xlab('')  +theme(axis.text=element_text(size=12),axis.title=element_text(size=14) ,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank())+ 
    ylab ( y_lab) 
  pl
}


plot_sample_distribution <- function (trees) {
  nex <- trees
  algn3 = data.frame(seq_id = nex$tip.label) %>% 
    separate(seq_id, c("Seq_ID", "GISAID", "Date", "Date_dec", "Region"), "[|]") %>% 
    mutate(Region = ifelse(Region == "_Il", "Local", "Global")) %>% 
    mutate(Date = as.Date(date_decimal(as.numeric(Date_dec)), "%Y-%m-%d"))
  pl = ggplot(algn3, aes(x = Date, color = Region, fill = Region)) + 
    geom_histogram(aes(y = ..density..), position = "identity", 
                   alpha = 0.5, bins = as.numeric(max(algn3$Date) - min(algn3$Date)) + 1) + geom_density(alpha = 0.4) + 
    theme_minimal() + theme(legend.title=element_text(size=14), 
                  legend.text=element_text(size=12),
                  axis.text=element_text(size=14),
                  axis.title=element_text(size=14)) + scale_color_manual(values = c("#999999", "#E69F00")) + 
    scale_fill_manual(values = c("#999999", "#E69F00")) + labs(x = "", y = "Sampling density")
  return(pl)
}
