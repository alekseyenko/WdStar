library(dplyr)
library(shiny)
library(glue)
library(DT)
library(futile.logger)
library(ggplot2)

data_controls_box_height <- 210

# UI
analysisPanelUI <- function(id) {
  ns = NS(id)

  tagList(
    fluidRow(h1("Please wait for plots to appear, notification is at bottom of screen")),
    fluidRow(div(id="xdata_controls_panel",
       box(width = 4,
           background = "olive",
           height = data_controls_box_height,
           selectInput(ns("distanceMethod"), "Distance",
                 distance_choices, selected = "jsd")),
       box(width = 4, background="olive",
           height = data_controls_box_height,
           uiOutput(ns("strataFactor")),
           uiOutput(ns("mainFactorDropdown"))))),
        fluidRow(
          box(title = "Ordination Plot", width = 6, plotOutput(ns("plot"))),
          box(title = "Plot", width = 6, plotOutput(ns("newplot")))))
}


#' Analysis Panel module server-side processing
#'
#' @param input,output,session standard \code{shiny} boilerplate
#'
#' @return list with following components
#' \describe{
#'   \item{globalData}{reactive character indicating x variable selection}
#' }
analysisPanel <- function(input, output, session, 
                          globalData) {

  ns <- session$ns
  globalData <- globalData
  sampleMainFactor <- reactiveVal(NULL)
  selectedFactorChoices <- reactiveVal(NULL)

  observeEvent(input$mainFactor, {
  })
  

  availableFactors <- reactive({
    av <- globalData()[["factor_metadata"]]
    af <- av %>% 
            dplyr::filter(ready == TRUE) %>% 
            pull(name)
    af
  })

  # output$mainFactor ----
  output$mainFactorDropdown = renderUI({
    af <- availableFactors()
    if(length(af) > 0) {
      selectInput(ns("mainFactor"), 'Main factor', af )
    } else {
      selectInput(ns("mainFactor"), 'Main factor', 
                  c("Missing") )
    }
  })

  # output$strataFactor ----
  output$strataFactor = renderUI({
    af <- append("None", availableFactors())
    if(length(af) > 0) {
      selectInput(ns('strataFactor'), 'Strata factor', af )
    } else {
      selectInput(ns('strataFactor'), 'Strata factor', 
                  c("Missing"))
    }
  })


  # physeqDataFactor ----
  physeqDataFactor <- reactive({
     req(input$mainFactor)

     #physeq <- physeqData()

     sd <- globalData()[["sample"]]
     q <- sd[ , input$mainFactor]
     
     # save main factor
     sampleMainFactor(q)
     #sampleMainFactor(sample_data(physeq)[[input$mainFactor]])

     #soilrep
     physeq
  })

  # distanceMatrices ----
  distanceMatrices <- reactive({
    req(physeqDataFactor())

    physeq <- physeqDataFactor()
    dist_matrices = distance(physeq, method=c(input$distanceMethod))
  })

  # output$plot ----
  output$plot <- renderPlot({
    req(physeqDataFactor(), input$mainFactor)
    withProgress(message = 'Creating ordination plot', value = 0, {

     physeq <- physeqDataFactor()
     dist_matrices = distance(physeq, method=c(input$distanceMethod))


     ordination_list = ordinate(physeq, method="MDS",
            distance = distanceMatrices())

     p <- plot_ordination(physeq = physeq,
         ordination = ordination_list,
         type = "samples",
         color = input$mainFactor ) +
         ggplot2::theme_minimal()

    flog.info("analysis : plot D")
     p <- p + ggplot2::scale_color_discrete(labels = paste(levels(sampleMainFactor()),
                                          table(sampleMainFactor())))
                  })

     return(p)
  })


  # display new plot ----
  output$newplot <- renderPlot({
     req(physeqDataFactor(), input$mainFactor)

     withProgress(message = 'Creating plot', value = 0, {
     physeq <- physeqDataFactor()
     dist_matrices = distance(physeq, method=c(input$distanceMethod))


     phy <- physeq
     jsd.dist = dist_matrices
     jsd.pco = ade4::dudi.pco(ade4::cailliez(jsd.dist), scannf = FALSE, nf = 2)

     if(input$strataFactor != "None") {
       location.wca = wca(jsd.pco, sample_data(phy)[[input$mainFactor]],
                          scannf = FALSE, nf = 2)
       s.class(location.wca$li,
               sample_data(phy)[[input$strataFactor]])

       #ade4::s.class(jsd.pco$li, interaction(sample_data(phy)[[input$mainFactor]],
      #                                       #sample_data(phy)[[input$strataFactor]]
       #                                      ))
     } else {
       ade4::s.class(jsd.pco$li, sample_data(phy)[[input$mainFactor]])
     }
   })

  })

  outin <- reactiveValues(inputs = NULL)
  return(outin)
}
