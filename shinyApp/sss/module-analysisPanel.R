library(dplyr)
library(shiny)
library(glue)
library(ggplot2)
library(DT)
library(futile.logger)

# UI
analysisPanelUI <- function(id) {
  ns = NS(id)

  tagList(uiOutput(ns("analysis")))
}

# Server
analysisPanel <- function(input, output, session, data) {

  data <- data

  sampleMainFactor <- reactiveVal(NULL)
  selectedFactorChoices <- reactiveVal(NULL)

  data_controls_box_height <- 210

  numTests <- reactive({ 0 })

  availableFactors <- reactive({
    md <- data()[["factor_metadata"]]
    af <- md %>% dplyr::filter(ready == TRUE) %>% pull(name)
  })

  observeEvent(input$runTestsBtn1, {
    show(id = "running-tests", anim = TRUE, animType = "fade")
    runTestsAction()
    hide(id = "running-tests", anim = TRUE, animType = "fade")
    updateTabsetPanel(session, inputId="sidebar", selected="tests")
  })

  # output$mainFactor ----
  output$mainFactor = renderUI({
    af <- availableFactors()
    if(length(af) > 0) {
      ns <- session$ns
      selectInput(ns('mainFactor'), 'Main factor', af )
    }
  })

  # output$strataFactor ----
  output$strataFactor = renderUI({
    af <- append("None", availableFactors())
    if(length(af) > 0) {
      ns <- session$ns
      selectInput(ns('strataFactor'), 'Strata factor', af )
    }
  })

  # output$physeqData ----
  physeqData <- reactive({
     req(data())

     sdata <- data()[["sample"]]
     asv.t <- data()[["genomic"]]

     phyloseq( otu_table(asv.t, taxa_are_rows = F), sample_data(sdata))
  })

  # physeqDataFactor ----
  physeqDataFactor <- reactive({
     req(input$mainFactor)

     physeq <- physeqData()
     sampleMainFactor(sample_data(physeq)[[input$mainFactor]])

     physeq
  })

  # distanceMatrices ----
  distanceMatrices <- reactive({
    req(physeqDataFactor())

    physeq <- physeqDataFactor()
    dist_matrices = distance(physeq, method=c(input$distanceMethod))
  })

  observeEvent(input$clearTestsBtn, {
     testsDT(NULL)
  })

  # ADD TEST BUTTON ----
  observeEvent(input$addTestBtn, {
     req(input$mainFactor)

     flog.info("ADD TEST BUTTON")
     aRow <- data_frame(test_id = numTests() + 1,
                        strata_factor = input$strataFactor,
                        main_factor = input$mainFactor,
                        distance = input$distanceMethod,
                        status = "not run")

     newTable <- bind_rows(testsDT(), aRow)
     testsDT(newTable)
  })


  # analysis tab ----
  output$analysis <- renderUI({
    flog.info("117: analysis render ui")

    d <- data()
    print(names(d))

    sd <- d[["sample"]]
    gd <- d[["genomic"]]

    if(!is.null(sd) && !is.null(gd)) {
      ns <- session$ns
      flog.info("120: analysis render ui: all data available")
      mainPanel(
      fluidRow(div(id="analysis_panel",
        column(4, selectInput(ns("distanceMethod"), "Distance",
                              distance_choices, selected = "jsd")),
        column(4, uiOutput(ns("mainFactor"))),
        column(4, uiOutput(ns("strataFactor")))),

      fluidRow(
        column(2, selectInput(ns("plotType"), "Plot Type",
                                    c("basic","newplot","both"))),
        column(10, uiOutput(ns("plot")))))
      )
    } else {
      flog.info("analysis render ui: no data available")
      div(h1("No data loaded."))
    }
  })

  output$plot <- renderUI({
      ns <- session$ns

      if (input$plotType == "both") {
        fluidRow(column(6, plotOutput(ns("origplot"))),
                 column(6, plotOutput(ns("newplot"))))
      } else if (input$plotType == "basic") {
        plotOutput(ns("origplot"))
      } else {
        plotOutput(ns("newplot"))
      }
  })

  # output$origplot ----
  output$origplot <- renderPlot({
    flog.info(glue("analysis : origplot {input$mainFactor}"))

    req(physeqDataFactor(), input$mainFactor)

    physeq <- physeqDataFactor()
    # MDS - calls ape::pcoa method
    ordination_list = phyloseq::ordinate(physeq,
                                        method="MDS",
            distance = distanceMatrices())

    p <- plot_ordination(physeq = physeq,
         ordination = ordination_list,
         type = "samples",
         color = input$mainFactor ) +
         theme_minimal() +
         # update the legend
         scale_color_discrete(labels =
                paste(levels(sampleMainFactor()), table(sampleMainFactor()))) +
         # add larger points
         geom_point(size = 7, alpha = 0.75)

    return(p)
  })


  # display new plot ----
  output$newplot <- renderPlot({
     req(physeqDataFactor(), input$mainFactor)

     physeq <- physeqDataFactor()

     flog.info(glue::glue("newplot: {input$mainFactor} {input$strataFactor} {input$distanceMethod}"))

     jsd.dist = phyloseq::distance(physeq, method=c(input$distanceMethod))
     # cailliez computes the smallest positive constant
     # that makes Euclidean a distance matrix
     # and applies it
     #
     # dudi.pco - principal coordinates analysis of
     # Euclidean distance matrix, returns objects of
     # class pco and dudi
     jsd.pco = ade4::dudi.pco(ade4::cailliez(jsd.dist), scannf = FALSE, nf = 2)

     if(input$strataFactor != "None") {
       # within class analysis
       #  Performs a particular case of an
       #  Orthogonal Principal Component Analysis with respect to
       #  Instrumental Variables (orthopcaiv),
       #  in which there is only a single factor as covariable
       w = ade4::wca(jsd.pco, sample_data(physeq)[[input$strataFactor]],
                          scannf = FALSE, nf = 2)
       li <- w$li
       colnames(li) <- c("A1", "A2")
     } else {
       li <- jsd.pco$li
     }

     scale_labels <- paste(levels(sampleMainFactor()), table(sampleMainFactor()))
     yyy <- sample_data(physeq)[[input$mainFactor]]

     zzz <- ggplot(li) +
            theme_minimal() +
            geom_point(aes(x = A1, y = A2, col = yyy), size = 7, alpha=0.75)  +
            scale_color_discrete(name = input$mainFactor, labels = scale_labels)

     #ade4::s.class(li, sample_data(physeq)[[input$mainFactor]])
     zzz
  })

  outin <- reactiveValues(inputs = NULL)
  return(outin)
}
