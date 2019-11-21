#dependencies
library(shiny)
library(shinyBS)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)
library(DT)
#sqlite libraries
library(DBI)
library(RSQLite)

source("utility/readConfig.R")
#get database from config file
databases <- config::get("SQLiteDB")

jsCode <- "shinyjs.progressToggle = function(params){if(params == 1){console.log();document.querySelectorAll(\"[id^='spinner-']\").forEach(el => {el.style = 'display:none'});} else{console.log();document.querySelectorAll(\"[id^='spinner-']\").forEach(el => {el.style = ''});}}"

#create handle to sqlite file
dbhandle <- dbConnect(RSQLite::SQLite(), databases[1])

#get the sample ids in the database
res <- dbGetQuery(dbhandle, "SELECT sample_id from GENOM_SAMPLE")
options(shiny.sanitize.errors = TRUE)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsCode, functions="progressToggle"),
  theme = shinytheme("united"),
  # Application title
  titlePanel("TarPan Viewer"),
  img(src = "logo.jpg", align = "right"),
  tags$head(tags$link(rel = "shortcut icon", href = "favicon.ico")),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput("database", width = '100%',
                           "Database:", selected = as.character(databases[1]),
                           c(unique(as.character(databases)))),
      checkboxInput("noprogress", "Hide progress indicator", value=0),
      h4("Data Filtering Options"),
      checkboxInput("type", "Show all chromosomes", value = 1, width = NULL),
      checkboxInput("hideXY", "Hide X/Y chromosomes", value= 1, width = NULL),
      conditionalPanel(
        condition = "input.type != 1",
      sliderInput("chromnum",
                  "Chromosome Number",
                  min = 1,
                  max = 24,
                  value = 1,
                  step = 1
                  )
      ),
      checkboxInput("chromnorm", "Use specific chromosomes for normalization", 
                    value = FALSE, width = NULL),
      conditionalPanel(
        condition = "input.chromnorm == 1",
        fluidRow(
          column(4,
            checkboxInput("chrom1", "1", value = FALSE, width = NULL),
            checkboxInput("chrom2", "2", value = FALSE, width = NULL),
            checkboxInput("chrom3", "3", value = FALSE, width = NULL),
            checkboxInput("chrom4", "4", value = FALSE, width = NULL),
            checkboxInput("chrom5", "5", value = FALSE, width = NULL),
            checkboxInput("chrom6", "6", value = FALSE, width = NULL),
            checkboxInput("chrom7", "7", value = FALSE, width = NULL),
            checkboxInput("chrom8", "8", value = FALSE, width = NULL)
          ),
          column(4,
            checkboxInput("chrom9", "9", value = FALSE, width = NULL),
            checkboxInput("chrom10", "10", value = FALSE, width = NULL),
            checkboxInput("chrom11", "11", value = FALSE, width = NULL),
            checkboxInput("chrom12", "12", value = FALSE, width = NULL),
            checkboxInput("chrom13", "13", value = FALSE, width = NULL),
            checkboxInput("chrom14", "14", value = FALSE, width = NULL),
            checkboxInput("chrom15", "15", value = FALSE, width = NULL),
            checkboxInput("chrom16", "16", value = FALSE, width = NULL)
          ),
          column(4,
            checkboxInput("chrom17", "17", value = FALSE, width = NULL),
            checkboxInput("chrom18", "18", value = FALSE, width = NULL),
            checkboxInput("chrom19", "19", value = FALSE, width = NULL),
            checkboxInput("chrom20", "20", value = FALSE, width = NULL),
            checkboxInput("chrom21", "21", value = FALSE, width = NULL),
            checkboxInput("chrom22", "22", value = FALSE, width = NULL)
          )
        )
      ),
      checkboxInput("show_blacklist", "Show blacklist genes", 
                    value = FALSE, width=NULL),
      conditionalPanel(
        "$('li.active a').first().html() == 'Mutations' | $('li.active a').first().html() == 'RCircos'",
        h4("Mutation Options"),
        checkboxInput("all_mutation_information", 
                      "Show all mutation information", 
                      value = FALSE, width = NULL),
        checkboxInput("show_failed_mutations", 
                      "Show mutations that failed filter", 
                      value = FALSE, width = NULL),
        checkboxInput("hide_blank_mutations", 
                      "Hide blank mutations", 
                      value = FALSE, width = NULL),
        conditionalPanel(
          "$('li.active a').first().html() == 'RCircos'",
          checkboxInput("hide_circos_muts", "Hide Mutations", 
                       value = FALSE, width = NULL)
        )
      ),
      conditionalPanel(
        "$('li.active a').first().html() == 'Structural Variants' | $('li.active a').first().html() == 'RCircos'",
        h4("Structural Variant Options"),
        checkboxInput("show_failed_structvars", "Show SVs that failed filter", value = FALSE, width = NULL),
        conditionalPanel(
          "$('li.active a').first().html() == 'RCircos'",
          checkboxInput("show_inter_svs", "Show Inter-chromosomal SVs",value = FALSE, width = NULL)
        )
      ),
      conditionalPanel(
        "$('li.active a').first().html() == 'Copy Number Plot'",
        h4("Plotting Options"),
        conditionalPanel(
          condition = "input.type != 1",
          checkboxInput("continuous", "Continuous output", 
                        value = FALSE, width=NULL)
        ),
        conditionalPanel(
          condition = "input.type != 1 & input.continuous == 1",
          checkboxInput("showlast", "Indicate first gene position")
        ),
        checkboxInput("pandq", "Show p and q boundaries", 
                      value = FALSE, width = NULL),
        checkboxInput("nolim", "Change max limit", 
                      value = FALSE, width = NULL),
        conditionalPanel(
          condition = "input.nolim == 1",
          textInput("maxlim", label = "Maximum mean ratio", value = "0")
        ),
        
        checkboxInput("normline", "Show normal line", 
                      value = FALSE, width = NULL),
        conditionalPanel(
          condition = "input.type != 1",
          checkboxInput("showsnps", "Show SNPs for these regions", 
                        value = FALSE, width = NULL)
        )
      ),
      conditionalPanel(
        condition = "0 == 1",
        textInput("currentchrom", label = "", value = -1),
        textInput("currentregion", label = "", value = "ALL")
      ),
      conditionalPanel(
        "$('li.active a').first().html() == 'Copy Number Plot'",
        downloadButton('downloadCopyNumberPlot', 'Download CN Plot'),
        conditionalPanel(
          condition = "input.showsnps == 1",
          br(),
          downloadButton('downloadSNPPlot', 'Download SNP Plot')
        )
      ),
      conditionalPanel(
        "$('li.active a').first().html() == 'RCircos'",
        downloadButton('downloadRCircosPlot', 'Download RCircos')
      ),
      conditionalPanel(
        "$('li.active a').first().html() == 'Auto CN'",
        downloadButton('downloadDataCN', 'Download')
      ),
      conditionalPanel(
        "$('li.active a').first().html() == 'Auto CN Groups'",
        downloadButton('downloadDataCNG', 'Download')
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      width = 10,
      selectInput("sample", width = '100%',
                           "Sample:", selected = as.character(res$sample_id[1]),
                           c(unique(as.character(res$sample_id)))),
      tabsetPanel(type = "tabs", 
                  tabPanel("Copy Number Plot", plotOutput("plot") %>% withSpinner(color="#e41a1c"), 
                           plotOutput("plot2")),
                  tabPanel("Mutations", DT::dataTableOutput("mutations_table")  %>% withSpinner(color="#e41a1c")),
                  tabPanel("Structural Variants", 
                           DT::dataTableOutput("manta_table") %>% withSpinner(color="#e41a1c")),
                  tabPanel("QC Metrics", DT::dataTableOutput("metrics_table") %>% withSpinner(color="#e41a1c")),
                  tabPanel("Auto CN", DT::dataTableOutput("AutoCN") %>% withSpinner(color="#e41a1c")),
                  tabPanel("Auto CN Groups", DT::dataTableOutput("AutoCNGroups") %>% withSpinner(color="#e41a1c")),
                  # tabPanel("Debug", verbatimTextOutput("debug")),
                  tabPanel("RCircos", plotOutput("rcircos") %>% withSpinner(color="#e41a1c"))
      )
    )
  )
))
