library(dplyr)
library(shiny)
library(glue)
library(DT)
library(futile.logger)

# UI
testsPanelUI <- function(id) {
  ns = NS(id)

  tagList(uiOutput(ns("tests")))
}

# Server
testsPanel <- function(input, output, session, computed_details, sdata, 
                       allDataAvailable, numTests,
                       testsData, testsDT) {

  numTests <- numTests
  allDataAvailable <- allDataAvailable
  computed_details <- computed_details
  testsData <- testsData

  selectedFactorChoices <- reactiveVal(NULL)
  my_sample_data <- sdata

  output$tests <- renderUI({
    flog.info("40: model-testsPanel::tests - renderUI - TESTS")
    flog.info(glue("=== num Tests {allDataAvailable()} #tests: {numTests()}"))
    if(allDataAvailable()) {
      if ( numTests() > 0 ) {
        ns <- session$ns
        fluidRow(id="test_results_panel",
           box(width = 12, background="teal",
             actionButton(ns("runTestsBtn"), "Run Tests"),
             actionButton(ns("deleteRowsBtn"), "Delete Selected Rows"),
             actionButton(ns("clearTestsBtn1"), "Clear All"),
             DTOutput(ns("testTable"))))
      } else {
         div(h1("Data loaded. No tests"))
      }
    } else {
      div(h1("No data loaded."))
    }
  })

  # computeTresFtn ----
  computeTresFtn <- function(strata, mainf, dist) {

    physeq <- physeqData()
    dist_matrices = distance(physeq, method = dist)

    tres = WdS.test(dist_matrices,
                   f = sample_data(physeq)[[mainf]],
                   nrep = input$numPermutations,
                   strata = sample_data(physeq)[[strata]]
                   )
  }

  # computeTres ----
  computeTres <- reactive({

    physeq <- physeqDataFactor()
    dist_matrices = distance(physeq, method=c(input$distanceMethod))

    tres = WdS.test(dist_matrices,
                   values$sampleMainFactor,
                   nrep = input$numPermutations,
                   sample_data(physeq)[[input$strataFactor]])

    if(!is.list(tres)) {
      message(str_c("tres results invalid", class(tres)))
      message(tres)
      return(NULL)
    }
    tres
  })

  # distTestOutput ----
  distTestOutput <- reactive({

    tres = computeTres()
    dist_matrices = distance(physeqDataFactor(), method=c(input$distanceMethod))

    tres$d = dist.cohen.d(dist_matrices, values$sampleMainFactor)
    #flog.info("distTestOutput")
    df <- as_data_frame(tres$d) %>% rename(dist_cohen_d = value)
  })

  # statisticalOutput ----
  statisticalOutput <- eventReactive(input$runPermutations, {
    st <- computeTres()
    if(is.list(st)) {
      df <- as_data_frame(st)
    }
  })


  sayecho <- function(tid, strata, mainf, dist) {
    status <- "SUCCESS"
    err <- NULL
    x <- tryCatch({
            computeTresFtn(strata, mainf, dist)
    }, warning = function(w) {
            flog.warn(glue::glue("warning: {w}"))
            NULL
    }, error = function(e) {
            err = e
            if("message" %in% names(e)) {
               flog.error(e$message)
               flog.error(glue::glue("error : {e$message} {e$call}"))
            }
            err
            NULL
    })

  if(!is.null(x)) {
    # p.value, statistic, nrep
    data_frame(
               test_id = tid,
               strata_factor = strata,
               main_factor = mainf,
               p.value = x[["p.value"]],
               nrep = x[["nrep"]],
               statistic = x[["statistic"]],
               status = "SUCCESS")
  } else {
    flog.error(glue::glue("error : {err}"))
    data_frame(
               test_id = tid,
               strata_factor = strata,
               main_factor = mainf,
               p.value = 0,
               nrep = 0,
               statistic = 0,
               status = str_c("ERROR"))
  }
  }

  # RUN TESTS BUTTON ----
  runTestsAction <- reactive({
    nt <- numTests()
    n <- nt

    flog.info(str_c("Running tests: ", nt, " ", nrow(testsDT())))
    flog.info(testsDT())

    w_df = map(1:n, ~
                sayecho(testsDT()$test_id[.],
                        testsDT()$strata_factor[.],
                        testsDT()$main_factor[.],
                        testsDT()$distance[.]))  %>%
    bind_rows()

    a_df <- testsDT() %>%
            select(test_id, distance) %>%
            left_join(w_df, by="test_id")

    flog.info(str_c("results: ", nrow(a_df)))
    flog.info(a_df)
    if(nrow(a_df) != nrow(w_df)){
      browser()
    }

    testsDT(a_df)
  })


  observeEvent(input$runTestsBtn, {
    show(id = "running-tests", anim = TRUE, animType = "fade")
    runTestsAction()
    hide(id = "running-tests", anim = TRUE, animType = "fade")
  })

  # DELETE ROWS BUTTON ----
  observeEvent(input$deleteRowsBtn, {
    if(!is.null(input$testTable_rows_selected)) {
      nn = as.numeric(input$testTable_rows_selected)
      dt <- testsDT()
      testsDT(dt[-nn, ])
    }
  })

  observeEvent(input$clearTestsBtn, {
      testsDT(NULL)
  })

  observeEvent(input$clearTestsBtn1, {
      testsDT(NULL)
  })

  output$testTable1 <- renderDT({
    req(testsData())

    if(nrow(testsData()) > 0) {
      #datatable(testsData(), options = list(dom='t')) %>%
      t <- testsDT() %>%
          mutate(test=str_c(distance, strata_factor, main_factor, sep="/")) %>%
          select(test)
      datatable(t, options = list(dom='t'))
    }
  })

  output$testTable <- renderDT({
    req(testsData())

    flog.info("74: testsPanel::output$testTable ==> TEST TABLE")
    if(nrow(testsData()) > 0) {
            flog.info("HERE IS THE TESTDATA")
            flog.info(testsData())
    datatable(testsData(), options = list(dom='t')) %>%
        formatStyle('status',
                    color = styleEqual(c("SUCCESS", "ERROR"), c("lightgreen","red")),
                    backgroundColor = styleEqual(c("SUCCESS", "ERROR"), c("darkgreen","pink")))
    }
  })

  outin <- reactiveValues(inputs = NULL)
  return(outin)
}
