library(dplyr)
library(shiny)
library(glue)
library(readr)
library(DT)
library(shinywdstar)
library(futile.logger)
library(skimr)

source("module-assignFactorLevels.R")

# UI
factorPanelUI <- function(id) {
  ns = NS(id)

  tagList(
    fluidRow(
      column(4, uiOutput(ns("factorDropdown")),
                uiOutput(ns("availableFactorDropdown"))),
      column(4, uiOutput(ns("variableDetailsDND"))),
      column(4, uiOutput(ns("editVariableDetails")))
      ),
    fluidRow(column(12,
      DT::dataTableOutput(ns("factorsDataTable")))))
}

# Server
factorPanel <- function(input, output, session, data) {

  data <- data

  #fp <- callModule(assignFactorLevels, "fp", data)
  #output$result <- renderPrint(str(input$dragvars))

  # Drop down containing potential factors
  output$factorDropdown <- renderUI({
    if(is.null(data()[["factor_metadata"]])) {
      return(p("No data to display"))
    }

    d <- data()
    af <- d[["factor_metadata"]] %>%
            dplyr::filter(ready == FALSE)

    if(is.null(af)) {
      return(p("No factor data available. NO DETAILS"))
    }
    ns <- session$ns
    #print(str_c("factors: ", str_c(af, collapse = ", ")))
    selectInput(ns("editFactor"), 'Variables', af %>% pull(name))
  })




  output$availableFactorDropdown <- renderUI({
    if(is.null(data()[["factor_metadata"]])) {
      return(NULL)
    }

    d <- data()
    af <- d[["factor_metadata"]] %>%
            dplyr::filter(ready == TRUE)

    if(is.null(af)) {
      return(p("No factor data available. NO DETAILS"))
    }
    ns <- session$ns
    #print(str_c("factors: ", str_c(af, collapse = ", ")))
    selectInput(ns("availableFactor"), 'Available Variables', af %>% pull(name))
  })

  output$factorsDataTable <- DT::renderDataTable({
    if(is.null(data()[["factor_metadata"]])) {
      return(NULL)
    }
     d <- data()[["factor_metadata"]]
     sd <- data()[["sample"]]
     #get the 
     skim_df <- skim(sd) %>% 
       dplyr::filter(stat == "hist") %>%
       select(name = variable, Histogram = formatted)

     d1 <- d %>% left_join(skim_df, by="name") 
     
     #d2 <- d1 %>%
     #        knitr::kable("html") %>%
     #        kable_styling("striped", full_width = FALSE)
     d1
  })

  valsForSelected <- reactive({

    sd <- data()[["sample"]]
    if(is.null(sd)) {
       print(" no sd vd returning")
       return(NULL)
    }

    print(glue("sample class: {ncol(sd)}"))
    if(ncol(sd) < 1) {
       print(glue("ERROR missing sample data for selected variable {input$editFactor}"))
       flog.error(glue("ERROR sample data is not in correct format {input$editFactor}"))
       return(p(glue("Missing data for {input$editFactor}")))
    }

    vals <- sd %>% pull(input$editFactor)

  })

  output$variableDetailsDND <- renderUI({

    print("VARIABLE DETAILS {input$editFactor}")
    vals <- valsForSelected()

    ns <- session$ns
    if (is.factor(vals)) {
      choices <- levels(vals)
      print(choices)
      dnd <- esquisse::dragulaInput(inputId = ns("dragvars"),
      #dnd <- shinywdstar::dragulaInput(inputId = ns("dragvars"),
         sourceLabel = "Levels",
         targetsLabels = c("Level 0", "Level 1"),
         targetsIds = c("level0", "level1"),
         choices = choices,
         badge = TRUE, width = "400px", 
         height = "100px", replace = FALSE) 
      return(dnd)
      #return(assignFactorLevelsUI(id = ns("fp")))
    }
    
    h1(glue("Non-factor {input$editFactor}"))
  })



  output$editVariableDetails <- renderUI({

    print(glue("edit VARIABLE DETAILS {input$editFactor}"))
    vals <- valsForSelected()

    ns <- session$ns
    if (is.factor(vals)) {
      return(tagList(
          textInput(ns("level0_label"), "Level 0 Label")
        , textInput(ns("level1_label"), "Level 1 Label")
        , actionButton(ns("applyFactorChanges"), "Apply")))
    }
    
  })

  outin <- reactiveValues(inputs = NULL)

  observeEvent(input$editFactor, {
    print(glue("editFactor changed: {input$editFactor}"))
    outin$inputs <- reactiveValuesToList(input)
  })

  observeEvent(input$applyFactorChanges, {
    print("FACTOR CHANGES")
    flog.info(glue("apply factor LABELS"))
    print(glue("apply factor: Labels 0= [{input$level0_label}] 1=[{input$level1_label}]"))

    if ("source" %in% names(input$dragvars)) {
      print("NOT APPLYING CHANGES: ALL VARIABLES MUST BE ACCOUNTED FOR ")
      print(glue("SOURCE {input$dragvars$source}"))
      #need to show an alert here
      return(NULL)
    }

    target <- input$dragvars$target
    print(str_c("LEVEL0:", str_c(target[["level0"]], collapse=", ")))
    print(str_c("LEVEL1:", str_c(target[["level1"]], collapse=", ")))
    print("NOT APPLYING CHANGES: NOT IMPLEMENTED")
    
    #outin$inputs <- reactiveValuesToList(input)
  })
  
  observeEvent(input$dragvars, {
    print(glue("dragvars changed: {input$dragvars}"))
    outin$inputs <- reactiveValuesToList(input)
  })

  return(outin)
}
