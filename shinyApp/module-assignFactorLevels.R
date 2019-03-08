library(dplyr)
library(shiny)
library(glue)
library(readr)
library(DT)
library(shinywdstar)
library(futile.logger)

# UI
assignFactorLevelsUI <- function(id) {
  ns = NS(id)
  uiOutput(ns("factorDragAndDrop"))
}

# Server
assignFactorLevels <- function(input, output, session, computed_details, sdata, gdata) {

  computed_details <- computed_details
  sample_data <- sdata
  genomic_data <- gdata

  # Drop down containing potential factors
  output$factorDragAndDrop <- renderUI({

    af <- computed_details()
    if(is.null(af)) {
      return(p("NO DETAILS"))
    }
    ns <- session$ns

    choices <- af %>% pull(name)
    #sapply(sample_data, is.factor)
    #sapply(sample_data, function(x) {length(levels(x))} )

    dnd <- 
       esquisse::dragulaInput(inputId = ns("dragvars"),
                                 sourceLabel = "Levels",
         targetsLabels = c("Level 0", "Level 1"),
         targetsIds = c("level0", "level1"),
         choices = choices,
         badge = TRUE, width = "400px", 
         height = "100px", replace = FALSE) 

       dnd
  })


  output$dnd_id <- renderPrint({ 
          ns <- session$ns
          ns("factorDragAndDrop")
  })

  outin <- reactiveValues(inputs = NULL)

  observeEvent(input$dragvars, {
    print(glue("dragvars changed: {input$dragvars}"))
    outin$inputs <- reactiveValuesToList(input)

    ns <- session$ns
    outin$inputs$dnd_id <- ns("factorDragAndDrop")
  })
  return(outin)
}
