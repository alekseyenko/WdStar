library(dplyr)
library(shiny)
library(glue)
library(readr)
library(DT)
library(shinywdstar)
library(futile.logger)
library(DiagrammeR)

# UI
statusPanelUI <- function(id) {
  ns = NS(id)
  tagList(
      #grVizOutput('diagram'),
      uiOutput(ns("statusPanel")))
}

# Server
statusPanel <- function(input, output, session, data) {

  data <- data

    spec <- "
   digraph {
       graph [overlap = true, rankdir = LR]
       node [shape = box, fontname = Helvetica, color = blue]
       edge [color = gray]
       Load; Preprocess; Analysis; Test
       Load->Preprocess                   
       Preprocess->Analysis [color = red]
       Analysis->Test
    }"
  
  # Drop down containing potential factors
  output$statusPanel <- renderUI({

    sd <- data()[["sample"]]
    gd <- data()[["genomic"]]
    if(is.null(sd) && is.null(gd)) {
      return(div("No data loaded."))
    }
    sample_msg <- ifelse(is.null(sd),"No sample data", glue("Samples: {nrow(sd)}"))
    g_msg <- ifelse(is.null(gd),"No genomic data", glue("Genomic: {nrow(gd)}"))

      div(str_c(sample_msg, " ", g_msg))
#      ,  grVizOutput("processGraph")
#      , DiagrammeROutput('processGraph', height = "50px", width = "auto")
  })

  output$diagram <- renderGrViz({
      grViz(spec)
      #grViz("digraph test {A; B; A-> B }")
  })
  return(NULL)
}
