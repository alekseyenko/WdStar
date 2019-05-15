library(dplyr)
library(shiny)
library(glue)
library(readr)
library(DT)
library(futile.logger)

#source("module-assignFactorLevels.R")

# UI
factorPanelUI <- function(id) {
  ns = NS(id)

  tagList(
    fluidRow(
      column(3, uiOutput(ns("factorDropdown"))),
      column(9, uiOutput(ns("factorDragAndDrop")))),
      fluidRow(column(12,
        DT::dataTableOutput(ns("factorsDataTable")))))
}

# Server
factorPanel <- function(input, output, session, globalData) {

  globalData <- globalData
  outin <- reactiveValues(inputs = NULL)

  selectedFactorChoices <- reactiveVal(NULL)

  output$result <- renderPrint(str(input$dragvars))

  # Drop down containing potential factors
  output$factorDropdown <- renderUI({
     af <- globalData()[["computed_details"]]
     if(is.null(af)) {
       return(p("No factor data available. NO DETAILS"))
     }
     ns <- session$ns
     selectInput(ns("editFactor"), 'Factor', af %>% pull(name))
  })

  # Drop down containing potential factors
  output$factorDragAndDrop <- renderUI({

    if(is.null(input$editFactor)) {
      return(NULL)
    }

    af <- globalData()[["computed_details"]]
    if(is.null(af)) {
      return(p("a1) NO DETAILS"))
    }
    ns <- session$ns

    print(paste0("SELECTED: ", input$editFactor))

    fmd <- globalData()[["factor_metadata"]]
    fmdthis <- dplyr::filter(fmd, name == input$editFactor)
    fchoices <- pull(fmdthis, unique_values)
    print(fchoices)
    them <- str_split(fchoices, ", ") %>% unlist()
    choices <- them

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
     globalData()[["factor_metadata"]]
  })


  observeEvent(input$editFactor, {
    mysd <- globalData()[["sample"]]
    if(is.null(mysd)) {
      return(NULL)
    }
    print(glue("editFactor changed: [{input$editFactor}]"))
    if(!input$editFactor %in% colnames(mysd)) {
       print("ERROR sample data is not in the correct format")
       futile.logger::flog.warn("ERROR sample data not in correct format")
       browser()
       return(NULL)
    }

    #####################################
    # set the tab based on the type of factor selected
    vals <- mysd %>% pull(input$editFactor)
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
    # set the data in factors metadata
    fmd <- globalData()[["factor_metadata"]]
    #fmda <- fmd %>% 
    outin$inputs <- reactiveValuesToList(input)
  })

  return(outin)
}
