library(dplyr)
library(shiny)
library(glue)
library(DT)
library(futile.logger)

# UI
analysisPanelUI <- function(id) {
  ns = NS(id)

  tagList(uiOutput(ns("analysis")))
}

# Server
analysisPanel <- function(input, output, session, computed_details, sampleData, 
                       allDataAvailable, availableFactors, numTests,
                       testsData, testsDT) {

  numTests <- numTests
  allDataAvailable <- allDataAvailable
  computed_details <- computed_details
  testsData <- testsData
  sampleMainFactor <- reactiveVal(NULL)

  selectedFactorChoices <- reactiveVal(NULL)

  data_controls_box_height <- 210

  data_controls_panel <- reactive({
    ns <- session$ns
     fluidRow(div(id="xdata_controls_panel",
       box(width = 4,
           background = "olive",
           height = data_controls_box_height,
           selectInput(ns("distanceMethod"), "Distance",
                              distance_choices, selected = "jsd")),
       box(width = 4, background="olive",
           height = data_controls_box_height,
           uiOutput(ns("strataFactor")),
           uiOutput(ns("mainFactor"))),
       box(id = ns("add_test_box"), width = 4, background= "teal",
           height = data_controls_box_height,
           actionButton(ns("addTestBtn"), "Add Test"),
           numericInput(ns("numPermutations"), "Number of Permutations",
                            value = 999, min = 1, max = 100000),
           if(numTests() > 0) actionButton(ns("runTestsBtn1"), "Run Tests"),
           if(numTests() > 0) actionButton(ns("clearTestsBtn"), "Clear Tests"),
           if(numTests() > 0) textOutput(ns("testTable2"))
           )))})


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
      selectInput('mainFactor', 'Main factor', af )
    }
  })

  # output$strataFactor ----
  output$strataFactor = renderUI({
    af <- append("None", availableFactors())
    if(length(af) > 0) {
      selectInput('strataFactor', 'Strata factor', af )
    }
  })


  # physeqDataFactor ----
  physeqDataFactor <- reactive({
     flog.info(glue("78: physeqDataFactor {input$mainFactor}"))
     req(input$mainFactor)

     flog.info(glue("81: physeqDataFactor {input$mainFactor}"))
     physeq <- physeqData()

     flog.info(glue("84: physeqDataFactor {input$mainFactor}"))
     # save main factor
     sampleMainFactor(sample_data(physeq)[[input$mainFactor]])
     flog.info(glue("87: physeqDataFactor {input$mainFactor}"))

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
    if(allDataAvailable()) {
      ns <- session$ns
      flog.info("120: analysis render ui: all data available")
      div(data_controls_panel(), 
        fluidRow(
          box(title = "Ordination Plot", width = 6, plotOutput(ns("plot"))),
          box(title = "Plot", width = 6, plotOutput(ns("newplot")))))
    } else {
      flog.info("analysis render ui: no data available")
      div(h1("No data loaded."))
    }
  })

  # output$plot ----
  output$plot <- renderPlot({
    flog.info(glue("analysis : plot {input$mainFactor}"))
    flog.info(glue("analysis : plot {is.null(physeqDataFactor())}"))
    req(physeqDataFactor(), input$mainFactor)
    flog.info("analysis : plot B")

     physeq <- physeqDataFactor()
     dist_matrices = distance(physeq, method=c(input$distanceMethod))

     #print(glue::glue("plot_ordination: {input$mainFactor} {input$distanceMethod}"))

    flog.info("analysis : plot C")
     ordination_list = ordinate(physeq, method="MDS",
            distance = distanceMatrices())

     p <- plot_ordination(physeq = physeq,
         ordination = ordination_list,
         type = "samples",
         color = input$mainFactor ) +
         theme_minimal()

    flog.info("analysis : plot D")
     p <- p + scale_color_discrete(labels = paste(levels(sampleMainFactor()),
                                          table(sampleMainFactor())))

     return(p)
  })


  # display new plot ----
  output$newplot <- renderPlot({
     flog.info("analysis render ui new plot A checking requirements")
     req(physeqDataFactor(), input$mainFactor)

     flog.info("analysis render ui new plot B requiements met")
     physeq <- physeqDataFactor()
     dist_matrices = distance(physeq, method=c(input$distanceMethod))

     #print(glue::glue("newplot: {input$mainFactor} {input$strataFactor} {input$distanceMethod}"))

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

  outin <- reactiveValues(inputs = NULL)
  return(outin)
}
