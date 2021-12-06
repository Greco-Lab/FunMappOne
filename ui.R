library(shiny)
library(shinyjs)
library(shinyBS)


jsCode <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}

shinyjs.addCustomTooltip = function(iVars) {
  var name = iVars[0]
  var msg = iVars[1]
  var tab_li = $('a[data-value=' + name + ']').parent();
  tab_li.attr('data-toggle', 'tooltip')
  tab_li.attr('title', msg)
}

shinyjs.removeCustomTooltip = function(name) {
  var tab_li = $('a[data-value=' + name + ']').parent();
  tab_li.removeAttr('data-toggle');
  tab_li.removeAttr('title');
}
"

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
}
.nav li a.disabled {
  background-color: #aaa !important;
  color: #333 !important;
  cursor: not-allowed !important;
  border-color: #aaa !important;
}"
fluidPage(
  useShinyjs(),
  extendShinyjs(text=jsCode,functions = c("disableTab","enableTab","addCustomTooltip","removeCustomTooltip")),
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
                          tags$h5("1. Input gene lists"),
                          wellPanel(
                            fluidRow(
                              column(4,
                                     radioButtons("organism","1) Organisms",
                                                  choices = c(human = "Human", mouse = "Mouse",rat = "Rat"),selected = "Human"),
                                     shinyBS::bsTooltip(id = "organism",title = "Note: select organism and gene ID before uploading the file",placement = "bottom")

                              ) ,
                              column(4,radioButtons("idtype","2) GeneID",
                                                    choices = c(symbols = "SYMBOL", ensemble = "ENSEMBL",entrez = "ENTREZID"),
                                                    selected = "SYMBOL"),
                                     shinyBS::bsTooltip(id = "idtype",title = "Note: select organism and gene ID before uploading the file",placement = "bottom")


                              ),
                              column(4,
                                     fileInput("file1", "3) Choose Excel File",
                                               multiple = FALSE,
                                               accept = c("text/csv/xlsx",
                                                          "text/comma-separated-values,text/plain/excel",
                                                          ".xlsx")))

                            )
                            # column(4,
                            # radioButtons("continuous","Plot modification",
                            #              choices = c(value = "continuous",
                            #                          sign = "discrete"),
                            #              selected = "continuous")
                            # )



                            # fluidRow(
                            #   column(6,
                            #          radioButtons("fileType","FileType",
                            #                       choices = c('Genes' = "GenesOnly", 'Genes and Modifications' = "genesFC"),
                            #                       selected = "GenesOnly") ),
                            #   column(6,radioButtons("disp", "Display",choices = c(Head = "head", All = "all"),selected = "head"))
                            #   ),
                            #   fileInput("file1", "Choose Excel File",
                            #             multiple = FALSE,
                            #             accept = c("text/csv/xlsx",
                            #                        "text/comma-separated-values,text/plain/excel",
                            #                        ".xlsx"))
                          ),
                          tags$h5("2. Functional annotation parameters"),
                          wellPanel(
                            fluidRow(
                              column(6,radioButtons("EnrichType","Select Functional Annotation",
                                                    choices = c(KEGG = "KEGG", REACTOME="REACTOME",GO = "GO"),
                                                    selected = "KEGG")
                              ),
                              column(6,radioButtons("GOType","Select GO",
                                                    choices = c(BP = "BP", CC="CC",MF = "MF"),
                                                    selected = "BP")
                              )),
                            fluidRow(column(6,
                                            radioButtons("MapValueType","Choose Values Type",
                                                         choices = c(Pvalue = "PVAL", GenesModifications="FC",GenesModifications_PValue  = "FCPV"),
                                                         selected = "FC")
                            ),
                            column(6,selectInput(inputId = "pvalueTh", label = "P-value threshold:",choices = list(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),selected = 0.05)


                            )),
                            fluidRow( #TOBECHANGED
                              column(6,checkboxInput("only_annotated", "Annotated genes only", value = TRUE)),
                              column(6,checkboxInput("only_significant", "Only significant", value = TRUE)),
                              column(6, sliderInput("min_intersection", "Minumum number of genes in the intersection:",
                                                    min = 0, max = 100, value = 0))
                            )

                          ),
                          tags$h5("3. Display parameters"),
                          wellPanel(
                            # Horizontal line ----
                            fluidRow(

                              column(6,radioButtons("aggregation","Aggregation Function",
                                                    choices = c(min = "min", max = "max",mean = "mean",
                                                                median = "median"),
                                                    selected = "mean")
                              ),
                              column(6,radioButtons("pcorrection","Correction Method",
                                                    #choices = c(gSCS = "analytical", fdr = "fdr", bonferroni = "bonferroni"),
                                                    choices = c(gSCS = "analytical", fdr = "fdr", bonferroni = "bonferroni", Nominal = "none"),     #TOBECHANGED

                                                    selected = "analytical"),
                                     shinyBS::bsTooltip(id = "pcorrection",
                                                        title = "Default is g:SCS. Check g:Profiler web page for more info",
                                                        placement = "top")

                              )),


                            fluidRow(column(6,radioButtons("continuous","Plot modification",
                                                           choices = c(value = "continuous",
                                                                       sign = "discrete"),
                                                           selected = "continuous")),
                                     column(6,actionButton("computePathways","Generate Map")))


                          )

                        ),

                        mainPanel(
                          wellPanel(
                            #fluidRow(tableOutput("contents")),
                            fluidRow(DT::dataTableOutput("contents")),
                            tags$hr(),

                            fluidRow(textOutput("updatedTable")),
                            fluidRow(DT::dataTableOutput("colSums"))
                          ),
                          wellPanel(
                            fluidRow(textOutput("updatedPat")),
                            fluidRow(DT::dataTableOutput("colSumsPat"))
                          )
                        )
                      )# end sidebar layout
             ),
             tabPanel("Plot Maps",value="PlotMaps",
                      sidebarLayout(
                        sidebarPanel(
                          wellPanel(
                            tags$h5("1. Data Selection"),
                            fluidRow(
                              selectInput(inputId = "level", label = "Browse hierarchy: choose a level",choices = list(1,2,3))
                            ),
                            fluidRow(
                              column(4,uiOutput("chose_lev1")),
                              shinyBS::bsTooltip(id = "chose_lev1",title = "Note: remove ALL from the list for specific selection.",placement = "top"),
                              column(4,uiOutput("chose_lev2")),
                              shinyBS::bsTooltip(id = "chose_lev2",title = "Note: remove ALL from the list for specific selection.",placement = "top"),
                              column(3,uiOutput("chose_lev3")),
                              shinyBS::bsTooltip(id = "chose_lev3",title = "Note: remove ALL from the list for specific selection.",placement = "top")

                            ),
                            fluidRow(
                              uiOutput("selectColumn"),
                              shinyBS::bsTooltip(id = "selectColumn",title = "Note: remove ALL from the list for specific selection.",placement = "top")

                            )
                          ),
                          wellPanel(
                            tags$h5("2. Plot section"),

                            fluidRow(
                              column(4,checkboxInput("doGrouping", "Show categories", value = TRUE)),
                              column(4,checkboxInput("aspectRatio", "Keep aspect ratio", value = FALSE)),
                              column(4,actionButton("do", "Plot Map")),
                              shinyBS::bsTooltip(id = "do",title ="NOTE: press the Plot Mat button every time you update the map!",placement = "bottom")

                            )
                          ),
                          wellPanel(
                            tags$h5("3. Download Selection"),
                            fluidRow(
                              column(4,textInput(inputId ="img_width", value = 15,label = "Width")), #width
                              column(4,textInput(inputId ="img_height", value = 30,label = "Height")),
                              column(4,downloadButton('downloadData')),
                              shinyBS::bsTooltip(id = "downloadData",title ="NOTE: when downloading, specify image size in inches ",placement = "bottom")
                            )
                          ),
                          wellPanel(
                            tags$h5("4. Clustering Selection"),
                            fluidRow(
                              column(4,uiOutput("nClust")),
                              column(4,selectInput("ClusterMethod","Select aggregation method",list("ward","complete","single"),selected = "complete")),
                              column(4,selectInput("Distance","Select distance",list("euclidean","jaccard","jaccard+euclidean"),selected = "jaccard"))
                            ),
                            fluidRow(
                              column(6,actionButton("doCluster", "Cluster Samples")),
                              column(6,actionButton("resetCluster","Reset Cluster"))
                            )
                          ),
                          tags$h5("5. Download enriched term lists as excel file"),
                          downloadButton("downloadEnrichedPathwayTables", "Download")

                        ),
                        mainPanel(
                          tags$h5("Use scrollbars to navigate and see the whole map"),

                          tabsetPanel(

                            tabPanel("Heatmap",fluidRow(column(12,align="left",shinycssloaders::withSpinner(plotOutput(outputId="heatmap"), type=6)))),
                            tabPanel("Enrichment Table",
                                     fluidRow(column(12,
                                        uiOutput("selectExperiment")
                                     )),
                                     fluidRow(
                                       column(12,
                                              DT::dataTableOutput("PAT_table")
                                       ))
                            ),
                            tabPanel(title="Heatmap Genes", id="hmGenes", value="hmGenes",
                                     fluidRow(column(4,
                                                     selectInput(inputId="levelGene", label="Choose a hierarchy level", choices=list(1,2,3))
                                     ),column(4,
                                              uiOutput("choosePath")
                                     ),column(4,
                                              selectInput(inputId="selScoreType", label="Show Values", choices=list("logFC"="lfc", "P-Value"="pval", "Combined"="comb"), selected="lfc")
                                     )),fluidRow(column(4,
                                                        #shinyBS::bsButton("doGeneHeatMap", label="Plot", style="danger", icon=icon("exclamation-circle"))
                                                        actionButton("doGeneHeatMap", "Plot")
                                     )),fluidRow(column(12,align="center",
                                                        shinycssloaders::withSpinner(plotOutput(outputId="heatmapGenes"), type=6)
                                     ))
                            ),
                            tabPanel("Clustering",plotOutput(outputId="hclust_plot", width = "100%"))
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
