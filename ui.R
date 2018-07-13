library(shiny)
library(shinyjs)

appCSS <- "
#loading-content {
position: absolute;
background: white !important;
opacity: 0.8;
z-index: 1000000;
left: 0;
right: 0;
top: 0;
bottom: 0;
font-size: 50px;
text-align: center;
color: #black;
}
#loading-gif { 
opacity: 0.8; 
display: block;
margin-left: auto;
margin-right: auto;
vertical-align: middle;
z-index: 1000000;
}"
fluidPage(
  useShinyjs(),
  #extendShinyjs(text=jsCode),
  #tags$head(tags$script(src="resizing.js")),
  inlineCSS(appCSS),
  hidden(div(id="loading-content",
             img(id="loading-gif", src="screen-loading.gif"),
             p(id="loadingText", "WORKING"),
             p("...")
  )),
navbarPage("FunMappOne",id = "page_id",
                tabPanel("Input",
                    sidebarLayout(
                      sidebarPanel(
                        wellPanel(
                        fluidRow(
                          column(6,
                                 radioButtons("fileType","FileType",
                                              choices = c('Genes' = "GenesOnly", 'Genes and Modifications' = "genesFC"),
                                              selected = "GenesOnly") ),
                          column(6,radioButtons("disp", "Display",choices = c(Head = "head", All = "all"),selected = "head"))
                          ),
                        #textInput("nSample", label = h3("Number of Samples"), value = "Enter text..."),
                        fileInput("file1", "Choose CSV File",
                                  multiple = FALSE,
                                  accept = c("text/csv/xlsx",
                                             "text/comma-separated-values,text/plain/excel",
                                             ".xlsx"))
                        
                        # Horizontal line ----
                         #tags$hr(),
                         
                        
                        # fluidRow(column(3,
                        #                 # Input: Checkbox if file has header ----
                        #                 checkboxInput("header", "Header", TRUE)),
                        #          column(3,      
                        #                 # Input: Select separator ----
                        #                 radioButtons("sep", "Separator",
                        #                              choices = c(Comma = ",",
                        #                                          Semicolon = ";",
                        #                                          Tab = "\t"),
                        #                              selected = ";")),
                        #          column(3,
                        #                 # Input: Select quotes ----
                        #                 radioButtons("quote", "Quote",
                        #                              choices = c(None = "",
                        #                                          "Double Quote" = '"',
                        #                                          "Single Quote" = "'"),
                        #                              selected = '"')),
                        #          column(3,
                        #                 # Input: Select number of rows to display ----
                        #                 radioButtons("disp", "Display",
                        #                              choices = c(Head = "head",
                        #                                          All = "all"),
                        #                              selected = "head"))
                        # )
                        ),
                        
                        
                        # # Horizontal line ----
                        # tags$hr(),
                        
                        
                        wellPanel(
                        fluidRow(
                          column(4,
                                 radioButtons("organism","Organisms",
                                              choices = c(human = "Human", mouse = "Mouse"),selected = "Human") 
                          ) ,
                          column(4,radioButtons("idtype","GeneID",
                                                choices = c(symbols = "SYMBOL", ensemble = "ENSEMBL",entrez = "ENTREZID"),
                                                selected = "SYMBOL")
                          ),
                          column(4,
                                 # radioButtons("continuous","Continuous or Discrete",
                                 #              choices = c(Continuous = "continuous",
                                 #                          Discrete = "discrete"),
                                 #              selected = "continuous")
                                 radioButtons("continuous","Plot modification",
                                              choices = c(value = "continuous",
                                                          sign = "discrete"),
                                              selected = "continuous")
                          ) 
                          
                        ),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        fluidRow(
                           
                              column(4,radioButtons("aggregation","Aggregation Function",
                                           choices = c(min = "min", max = "max",mean = "mean",
                                                       median = "median"),
                                           selected = "mean")  
                              ),
                              column(4,radioButtons("pcorrection","Correction Method",
                                                    choices = c(fdr = "fdr", bonferroni = "bonferroni",none = "none"),
                                                    selected = "fdr")
                                     ),
                             column(4,selectInput(inputId = "pvalueTh", label = "P-value threshold:",choices = list(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),selected = 0.05)

                            )
                            
                        ),
                        
                        tags$hr(),
                        
                        fluidRow(
                          column(3,radioButtons("EnrichType","Select Enrichment",
                                                choices = c(KEGG = "KEGG", REACTOME="REACTOME",GO = "GO"),
                                                selected = "KEGG")
                          ),
                          column(3,radioButtons("GOType","Select GO",
                                                choices = c(BP = "BP", CC="CC",MF = "MF"),
                                                selected = "BP")
                          ),
                          column(3,
                            #uiOutput("valueType")
                            radioButtons("MapValueType","Choose Values Type",
                                         choices = c(Pvalue = "PVAL", GenesModifications="FC",GenesModifications_PValue  = "FCPV"),
                                         selected = "FC")
                          ),
                          column(3,actionButton("computePathways","Generate Map"))

                        )
                       
                      )),
                        
                      mainPanel(
                        wellPanel(
                          #fluidRow(tableOutput("contents")),
                          fluidRow(dataTableOutput("contents")),
                          tags$hr(),
                          
                          fluidRow(textOutput("updatedTable")),
                          fluidRow(dataTableOutput("colSums"))
                          ),
                        wellPanel(
                          fluidRow(textOutput("updatedPat")),
                          fluidRow(dataTableOutput("colSumsPat"))
                        )
                      )      
                    )# end sidebar layout
           ),
           tabPanel("Plot Maps",value="PlotMaps",
              sidebarLayout(
                sidebarPanel(
                  wellPanel(
                    fluidRow(                      
                      selectInput(inputId = "level", label = "Browse hierarchy: choose a level",choices = list(1,2,3))
                    ),
                    fluidRow(
                      column(4,uiOutput("chose_lev1")),
                      column(4,uiOutput("chose_lev2")),
                      column(3,uiOutput("chose_lev3"))
                    ),
                    fluidRow(
                      uiOutput("selectColumn")
                    ),
                    fluidRow(
                      actionButton("do", "Plot Map"),
                      downloadButton('downloadData')
                    )
                  ),
                  wellPanel(
                    fluidRow(
                      column(4,uiOutput("nClust")),
                      column(4,selectInput("ClusterMethod","Select aggregation method",list("ward","complete","single"),selected = "complete")),
                      column(4,selectInput("Distance","Select distance",list("euclidean","jaccard","jaccard+euclidean"),selected = "jaccard"))
                    ),
                    fluidRow(
                      column(6,actionButton("doCluster", "Cluster Samples")),
                      column(6,actionButton("resetCluster","Reset Cluster"))
                    )
                  )
                  
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Heatmap", plotOutput(outputId="heatmap",height = 900)),
                    tabPanel("Clustering",plotOutput(outputId="hclust_plot",height = 900))
                  )
                  #fluidRow(plotOutput(outputId="heatmap",height = 900))
                )
              )
           )
)
)
# 
# #Define UI for application that plots random distributions
# shinyUI(fluidPage(
# 
#   # Application title
#   titlePanel("Kegg Hierarchy Visualization"),
# 
#   # Sidebar with a slider input for number of observations
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file1", "Choose CSV File",
#                 multiple = FALSE,
#                 accept = c("text/csv",
#                            "text/comma-separated-values,text/plain",
#                            ".csv")),
#       
#       # Horizontal line ----
#       tags$hr(),
#       fluidRow(column(4,
#                       # Input: Checkbox if file has header ----
#                       checkboxInput("header", "Header", TRUE)),
#                 column(4,      
#                       # Input: Select separator ----
#                       radioButtons("sep", "Separator",
#                                    choices = c(Comma = ",",
#                                                Semicolon = ";",
#                                                Tab = "\t"),
#                                    selected = ",")),
#                 column(4,
#                       # Input: Select quotes ----
#                       radioButtons("quote", "Quote",
#                                    choices = c(None = "",
#                                                "Double Quote" = '"',
#                                                "Single Quote" = "'"),
#                                    selected = '"'))
#       ),
#      
#       
#       # Horizontal line ----
#       tags$hr(),
#       
#       # Input: Select number of rows to display ----
#       radioButtons("disp", "Display",
#                    choices = c(Head = "head",
#                                All = "all"),
#                    selected = "head"),
#       
#       selectInput(inputId = "level", label = "Choose a level:",choices = list(1,2,3)),
#       uiOutput("chose_lev1"),
#       uiOutput("chose_lev2"),
#       uiOutput("chose_lev3"),
#       actionButton("do", "Click Me")
#       
#     ),
# 
#     # Show a plot of the generated distribution
#     mainPanel(
#       fluidRow(wellPanel(
#         plotOutput(outputId="heatmap",height = 900)
#         #plotlyOutput("heatmap")
#       )),
#       fluidRow(wellPanel(
#         tableOutput("contents")
#       )))
#     
#   )
# ))
