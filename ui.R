shinyUI(
  bootstrapPage(
    theme = shinythemes::shinytheme("spacelab"),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$link(rel = "shortcut icon", type = "image/x-icon",
                href = "favicon.ico")
    ),
    # include the js code
    shinyjs::useShinyjs(),
    includeScript("scripts.js"),
    column(12,
           HTML(
             "
             <meta name='keywords' content='infectious,disease,epidemiology,genome,coronavirus,SARS-CoV-2,phylodynamic,phylogenetics,lineage,mutation,Spike,UK'>
             <meta name='author' content='Rob Johnson'>
             <h1>Skygrowth dashboard</h1>
             "
           ),
      tabsetPanel(
        tabPanel("About",
          HTML(
             "
             <p>Samples from lineages are summarised over time and geography in the UK according to the 
             details of the constituent sequences. Scalable phylodynamic methods are used to estimate 
             time-scaled phylogenies and estimate effective population size through time.  
             </p>
             <p>Sequences and alignments from  
             <a href='https://gisaid.org/' target='_blank'>Gisaid</a>.</p>
             "
             )
        ),
        tabPanel("Methods",
                 HTML(
                   "
             <p>Importations are estimated with <a href='https://github.com/emvolz-phylodynamics/sarscov2Rutils' target='_blank'>sarscov2Rutils</a>. 
             </p>
             <p>Skygrowth curves show aspects of regional outbreaks constructed using <a href='https://github.com/emvolz-phylodynamics/sarscov2Rutils' target='_blank'>sarscov2Rutils</a>. 
             Graphs show the time-dependent effective reproduction number (the average number of secondary cases per primary case over time), the effective population size (the number 
             of individuals weighted by their contribution to producing the next generation) and the effective growth rate (the growth rate of the effective population size). </p>
             <p>Phylogenetic trees are constructed by <a href='https://cran.r-project.org/web/packages/treedater/index.html' target='_blank'>treedater</a>.</p>
             "
                 )
        ),
        tabPanel("References",
                 HTML(
                   "
                   <p>A Rambaut, EC Holmes, Á O’Toole et al. 
             <a href='https://doi.org/10.1038/s41564-020-0770-5' target='_blank'>A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology</a>
              Nat Microbiol (2020).</p>
             <p>EM Volz, SDW Frost 
             <a href='https://doi.org/10.1093/ve/vex025' target='_blank'>Scalable relaxed clock phylogenetic dating</a>
              Virus Evolution (2017).</p>
              <p>EM Volz, X Didelot
             <a href='https://doi.org/10.1093/sysbio/syy007' target='_blank'> Modeling the growth and decline of pathogen effective population size provides insight into epidemic dynamics and drivers of antimicrobial resistance</a>
              Systematic Biology (2018).</p>
              <p>EM Volz, V Hill, JT McCrone et al. 
             <a href='https://www.medrxiv.org/content/10.1101/2020.07.31.20166082v1' target='_blank'>Evaluating the effects of SARS-CoV-2 Spike mutation D614G on transmissibility and pathogenicity</a>
               (2020).</p>
             "
                 )
        )
        ),
      tags$hr(style="background: #cccccc; height: 4px;"),
      column(4,
        selectInput('ti_region', label = 'Select input', c('Sweden'),selected='Sweden')
        ),
      column(4,
             textOutput('ti_regionname')
      ),
      column(4),
      tags$hr(style="background: #cccccc; height: 4px;")
    ),
    column(12, id = "plot",
           tabsetPanel(
             tabPanel("Samples over time",
                      column(9, 
                             fluidRow( plotOutput( 'hist_by_location', width = "100%", height = "600px"), align="right")
                      )
             ),
             tabPanel("Importations",
                      #tableOutput("estimated_r_output")
                      fluidRow( plotOutput( 'imports', width = "100%", height = "600px"), align="right")
                      
             ),
             tabPanel("Skygrowth curves",
                      column(3, id = "menu",
                             div(id = "1.2",
                                 div(id = "incidence_data_type_error_box", class = "error_box",
                                     #selectInput('ti_curve', label = '', c('By lineage','By region'),selected='By lineage'),
                                     #checkboxInput('ti_ci', label = "Credible intervals", 1),
                                     checkboxInput(inputId='ti_log_size', label='Effective size on log scale', value=F)
                                     #, checkboxGroupInput(inputId='ti_DGX', label='Genotypes', choices = NULL, selected = NULL)
                                     #, checkboxGroupInput(inputId='ti_groups', label='', choices = NULL, selected = NULL)
                                     #, checkboxGroupInput(inputId='ti_filename', label='Lineages', choices = NULL, selected = NULL)
                                     
                                 )
                             )
                             
                      ),
                      #         downloadButton("save_plot", "Save Image"),
                      column(9, fluidRow( plotOutput( 'GR', width = "100%", height = "300px"))
                             , fluidRow( plotOutput( 'R', width = "100%", height = "300px" ) )
                             , fluidRow( plotOutput( 'Ne', width = "100%", height = "300px" ) ))#,
                      #column(3,plotOutput('legend'))
             ),
             tabPanel("Phylogenetic tree",
                      #downloadButton("save_incidence",
                      #               "Save Table"),
                      #selectInput('ti_filename_tree', label = 'Lineage', NULL),
                      #selectInput('ti_tree_colour', label = 'Colour by', choices=c('Genotype','Country','Location'),selected='Genotype'),
                      fluidRow( plotly::plotlyOutput( 'tree', width = "100%", height = "auto"), align="right")
                      #fluidRow( plotOutput( 'tree', width = "100%", height = "auto"), align="right")
             )
           )
    )
  )
)
