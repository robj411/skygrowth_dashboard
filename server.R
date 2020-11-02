library( lubridate )
library( shiny )
library(shinyjs)
library(grid)
#library(ape)
library(DT)
library(scales)
library(ggtree)
library( ggplot2 )
library(Hmisc)
library(plotly)
library(dplyr)
library(RColorBrewer)
library(plotrix)
library(sarscov2)
library(svDialogs)
library(tidyr)
source('shiny_functions.R')

#scale_fill_discrete <- function(...) {
#  scale_fill_brewer(..., palette="Accent")
#}


shiny::shinyServer(function(input, output, session) {
  
  preprepared <- sort(c('Sweden'))
  updateSelectInput(session, inputId='ti_region', label = 'Select input', c('Choose from file',preprepared),selected='Sweden')
  
  preprepared_files <- lapply(preprepared,function(region)
    readRDS(paste0('input_files/',region,'.Rds'))
    )
  names(preprepared_files) <- preprepared
  
  observeEvent({
    input$ti_region
  },{
    if(input$ti_region%in%preprepared){
      region <<- input$ti_region
      parms <<- preprepared_files[[region]]
      output$ti_regionname <- renderText({""})
    }else{
      path <- rstudioapi::selectFile(caption = "Select RDS File",
                                     filter = "RDS Files (*Rds)",
                                     path='.',
                                     existing = TRUE)
      parms <<- readRDS(path)
      if(!all(c('sequences','tree','lineages','imports')%in%names(parms))){
        print('break')
        break
      }
      splts <- strsplit(path,'/')[[1]]
      fname <- splts[length(splts)]
      default <- strsplit(fname,'\\.')[[1]][1]
      regions <- sort(table(parms$sequences$Location),decreasing=T)
      region <<- names(regions)[1]
      output$ti_regionname <- renderText({paste0("Location: ", region)})
      if(!region%in%preprepared){
        preprepared_files[[region]] <<- parms
        preprepared <<- sort(c(region,preprepared))
      }
      #region <<- dlg_input(message = "Enter location name", default = default, gui = .GUI)
      
    }
    updateSelectInput(session, inputId='ti_region', label = 'Select input', c('Choose from file',preprepared),selected=region)
  })
  
  
  
  
  ## initialise checkboxes
  #combined_labs <- c('Combined "D"','Combined "G"','Combined "D" and "G"')
  #updateCheckboxGroupInput(session, inputId='ti_DGX', label='Genotypes', choices = parms$geno_labs)
  #updateCheckboxGroupInput(session, inputId='ti_groups', label='', choices = combined_labs)
  #updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Lineages', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
  #updateSelectInput(session, inputId='ti_filename_tree', label = 'Lineage', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
  #updateSelectInput(session, inputId='ti_filename_region', label = 'Lineage', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
  #updateSelectInput(session, inputId='ti_region', label = 'Region', choices = unique(parms$sequences_geog$region), selected = 'london')
  
  ## initialise reactive values
  #re.filename <- reactiveVal( parms$filename )
  #re.filename_tree <- reactiveVal( parms$filename_tree )
  #re.hist <- reactiveVal( 'Combined' )
  
  
  ## observe input events
  ## geography ########################################################
  ## plot hist for geography
  observeEvent({input$ti_region},
               {
    output$hist_by_location <- renderPlot({
      plot_sample_distribution(parms$tree)
    })#, 
    #height = 3*hgt+500)
  })
  
  ## skygrowth ########################################################
  ## plot lineage or region?
  
  prep_plot <- function(metric='growth',ci=1,log_size=F){
    log_size <- metric=='Ne' & log_size
    #if(!inherits(sgs[[1]], "sarscov2skygrowth")) print(lin[[2]])
    p0 <- gg_sarscov2skygrowth(parms$lineages,metric=metric,log_size=log_size,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,'navyblue')
    return(p0)
  }
  
  ## plot
  observeEvent({
    input$ti_log_size
    input$ti_region
  },{
    #if(!is.na(lin[[2]][1])){
    ci=1
      output$GR <- renderPlot({
        suppressWarnings(plot( prep_plot(metric='growth',ci=ci)))
      })
      output$R <- renderPlot({
        suppressWarnings(plot( prep_plot(metric='R',ci=ci)))
      })
      log_size <- input$ti_log_size
      output$Ne <- renderPlot({
        suppressWarnings(plot( prep_plot(metric='Ne',ci=ci,log_size=log_size)))
      })
    #}
  })
  
  ## tree ########################################################
  
  
  ## plot tree with height depending on number of sequences
  observeEvent({input$ti_region},{
    anno <- 'Location'
    #hgt <- length(parms$tree[[l_ind]]$Ti)
    #output$tree <- renderPlot({
      output$tree <- renderPlotly({
        #quick_region_treeplot(parms$tree,'Il')
        quick_annotated_treeplot(parms$tree,annotation = anno) #%>% layout(height = 3*hgt+500)
      })#, 
    #height = 3*hgt+500)
  })
  
  ## imports #################################################
  observeEvent({input$ti_region},{
    output$imports <- renderPlot({
      par(mar=c(5,5,2,2))
      plot_importations(parms$imports,cex.axis=1.5,cex.lab=1.5) 
    })
  })
  
  output 
  
}
)
