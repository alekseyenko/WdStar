library(dplyr)
library(shiny)
library(glue)
library(DT)
library(shinywdstar)
library(futile.logger)

# UI
testsPanelUI <- function(id) {
  ns = NS(id)

  tagList(uiOutput(ns("tests")))
}

# Server
testsPanel <- function(input, output, session, data) {

  data <- data
  testsData <- reactiveVal(data.frame())
  numTests <- reactive({ nrow(testsData()) })
  dataLoaded <- reactive({ !is.null(data()[["sample"]]) &&
          !is.null(data()[["genomic"]]) })

  availableFactors <- reactive({ 
    md <- data()[["factor_metadata"]] 
    af <- md %>% dplyr::filter(ready == TRUE) %>% pull(name)
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


  data_controls_box_height <- 210

  data_controls_panel <- reactive({
    ns <- session$ns
     fluidRow(div(id="t2_xdata_controls_panel",
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


  output$tests <- renderUI({
    flog.info("40: model-testsPanel::tests - renderUI - TESTS")
    if(dataLoaded()) {
        ns <- session$ns
        tagList( data_controls_panel(),
                fluidRow(id = "test_results_panel",
           box(width = 12, background="teal",
             actionButton(ns("runTestsBtn"), "Run Tests"),
             actionButton(ns("deleteRowsBtn"), "Delete Selected Rows"),
             actionButton(ns("clearTestsBtn1"), "Clear All"),
             DTOutput(ns("testTable")))))
    } else {
      div(h1("No data loaded."))
    }
  })


  # output$physeqData ----
  physeqData <- reactive({
     sdata <- data()[["sample"]] 
     asv.t <- data()[["genomic"]] 

     phyloseq( otu_table(asv.t, taxa_are_rows = F), sample_data(sdata))
  })

  # computeTresFtn ----
  computeTresFtn <- function(physeq, strata, mainf, dist, nrep = 999) {
    dist_matrices  <- tryCatch({
      dist_matrices = distance(physeq, method = dist)
    }, error = function(e) {
      browser()
      err = e
      print(e$message)
    })

    #nrep = input$numPermutations,

    dstrata = NULL
    if(strata != "None") {
      dstrata = sample_data(physeq)[[strata]]
    }
    f = sample_data(physeq)[[mainf]]
    tres = shinywdstar::WdS.test(dist_matrices,
                   f = f, nrep = nrep, strata = dstrata)
    tres
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

  sayecho <- function(physeq, tid, strata, mainf, dist) {
    status <- "SUCCESS"
    err <- NULL

    nrep = 999
    x <- tryCatch({
      computeTresFtn(physeq, strata, mainf, dist, nrep = nrep)
    }, warning = function(w) {
      flog.warn(glue::glue("warning: {w}"))
      NULL
    }, error = function(e) {
      err = e
      print(e$message)
      browser()
      if("message" %in% names(e)) {
        print(e$message)
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

    flog.info(glue("Running tests: {numTests()}"))

    print("before running tests")
    t_df <- testsData()
    physeq <- physeqData()
    print(t_df)
    w_df = map(1:nrow(t_df), function(x) {
                       print(glue("RUNNING TESTS {x}"))
                       sayecho(physeq, 
                               t_df$test_id[x], 
                               t_df$strata_factor[x], 
                               t_df$main_factor[x],
                               t_df$distance[x]) }) %>% 
      bind_rows()

    a_df <- testsData() %>%
            select(test_id, distance) %>%
            left_join(w_df, by="test_id")

    flog.info(str_c("results: ", nrow(a_df)))
    flog.info(a_df)
    if(nrow(a_df) != nrow(w_df)){
      browser()
    }

    testsData(a_df)
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
      dt <- testsData()
      testsData(dt[-nn, ])
    }
  })

  observeEvent(input$clearTestsBtn, {
      testsData(data.frame())
  })

  observeEvent(input$clearTestsBtn1, {
      testsDate(data.frame())
  })

  output$testTable1 <- renderDT({
    req(testsData())

    if(nrow(testsData()) > 0) {
      #datatable(testsData(), options = list(dom='t')) %>%
      t <- testsData() %>%
          mutate(test=str_c(distance, strata_factor, main_factor, sep="/")) %>%
          select(test)
      datatable(t, options = list(dom='t'))
    }
  })

  output$testTable <- renderDT({
    req(testsData())

    if(nrow(testsData()) > 0) {
      datatable(testsData(), options = list(dom='t')) %>%
        formatStyle('status',
                    color = styleEqual(c("SUCCESS", "ERROR"), c("lightgreen","red")),
                    backgroundColor = styleEqual(c("SUCCESS", "ERROR"), c("darkgreen","pink")))
    }
  })


  observeEvent(input$runTestsBtn1, {
    show(id = "running-tests", anim = TRUE, animType = "fade")
    runTestsAction()
    hide(id = "running-tests", anim = TRUE, animType = "fade")
    updateTabsetPanel(session, inputId="sidebar", selected="tests")
  })


  observeEvent(input$clearTestsBtn, {
     testsData(data.frame())
  })

  # ADD TEST BUTTON ----
  observeEvent(input$addTestBtn, {
     req(input$mainFactor)
     
     aRow <- data_frame(test_id = numTests() + 1,
                        strata_factor = input$strataFactor,
                        main_factor = input$mainFactor,
                        distance = input$distanceMethod,
                        status = "not run")

     newTable <- bind_rows(testsData(), aRow)
     testsData(newTable)
  })


  outin <- reactiveValues(inputs = NULL)
  return(outin)
}
