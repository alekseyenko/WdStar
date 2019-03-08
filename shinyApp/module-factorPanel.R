library(dplyr)
library(shiny)
library(glue)
library(readr)
library(DT)
library(shinywdstar)
library(futile.logger)

source("module-assignFactorLevels.R")

# UI
factorPanelUI <- function(id) {
  ns = NS(id)

  tagList(
    fluidRow(
      column(4, uiOutput(ns("factorDropdown"))),
      column(4, assignFactorLevelsUI(id = ns("fp")))),
      fluidRow(column(12,
        DT::dataTableOutput(ns("factorsDataTable")))))
}

# Server
factorPanel <- function(input, output, session, computed_details, sdata, gdata) {

  computed_details <- computed_details

  selectedFactorChoices <- reactiveVal(NULL)

  my_sample_data <- sdata
  genomic_data <- gdata

  fp <- callModule(assignFactorLevels, "fp", 
                   selectedFactorChoices, 
                   my_sample_data, genomic_data)

  output$result <- renderPrint(str(input$dragvars))

  # Drop down containing potential factors
  output$factorDropdown <- renderUI({
     af <- computed_details()
     if(is.null(af)) {
       return(p("No factor data available. NO DETAILS"))
     }
     ns <- session$ns
     #print(str_c("factors: ", str_c(af, collapse = ", ")))
     selectInput(ns("editFactor"), 'Factor', af %>% pull(name))
  })

  # Drop down containing potential factors
  output$factorDragAndDrop <- renderUI({
    af <- computed_details()
    if(is.null(af)) {
      return(p("NO DETAILS"))
    }
    ns <- session$ns

    choices <- af %>% pull(name)
    #sapply(my_sample_data(), is.factor)
    #sapply(my_sample_data(), function(x) {length(levels(x))} )

    dnd <- 
       esquisse::dragulaInput(inputId = ns("dragvars"),
                                 sourceLabel = "SSSLevels",
         targetsLabels = c("Level 0", "Level 1"),
         targetsIds = c("level0", "level1"),
         choices = choices,
         badge = TRUE, width = "400px", 
         height = "100px", replace = FALSE) 

       dnd
  })

  output$factorsDataTable <- DT::renderDataTable({
     computed_details()
  })

  outin <- reactiveValues(inputs = NULL)

  observeEvent(input$editFactor, {
    print(glue("editFactor changed: {input$editFactor}"))
    if(!input$editFactor %in% colnames(my_sample_data())) {
       print("ERROR sample data is not in the correct format")
       return(NULL)
    }

    #####################################
    # set the tab based on the type of factor selected

    vals <- my_sample_data() %>% pull(input$editFactor)
    uvals <- unique(vals)
    if(is.factor(vals)) {
      prev <- selectedFactorChoices()
      tmp  <- data_frame(name=levels(vals))
      print(tmp)
      selectedFactorChoices(tmp)
    }
    outin$inputs <- reactiveValuesToList(input)
  })

  observeEvent(input$dragvars, {
    print(glue("dragvars changed: {input$dragvars}"))
    outin$inputs <- reactiveValuesToList(input)
  })

  return(outin)
}
