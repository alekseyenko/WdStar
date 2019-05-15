library(dplyr)
library(shiny)
library(glue)
library(DT)
library(futile.logger)

#' load panel to allow user to load data
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
loadPanelUI <- function(id) {
  ns = NS(id)

  tagList(uiOutput(ns("page")))
}

#' loadPanel module server-side processing
#'
#' @param input,output,session standard \code{shiny} boilerplate
#'
#' @return list with following components
#' \describe{
#'   \item{genomic}{reactive character indicating x variable selection}
#'   \item{sample}{reactive character indicating y variable selection}
#' }
loadPanel <- function(input, output, session) {

  ns = session$ns

  data <- reactiveVal(list(genomic = NULL, 
                           sample = NULL,
                           orig_factor_data = NULL,
                           factor_metadata = NULL,
                           orig_factor_metadata = NULL,
                           computed_details = NULL
                           ))

  physeqToSample <- function(physeq) {
    # -- convert the physeq to the same format as sdata
    dd <- physeq@sam_data@.Data
    names(dd) <- physeq@sam_data@names
    dd_df <- data.frame(dd)

    rownames(dd_df) <- physeq@sam_data@row.names
    dd_df
  }

  physeqToRawCount <- function(physeq) {
    t(physeq@otu_table)
  }

  generate_metadata_from_sample_data <- function(sd) {
    if(ncol(sd) > 0) {

      md <- map_df(colnames(sd), ~ handle_factor(., dplyr::pull(sd, .)))
    } else {
      NULL
    }
  }

  #this doesn't need to be reactive, it could be a list
  loadDemoData <- reactive({
    sdata <- physeqToSample(physeq)
    rcdata <- physeqToRawCount(physeq)

    orig <- data()
    orig[["sample"]] <- sdata
    orig[["genomic"]] <- rcdata

    md <- generate_metadata_from_sample_data(sdata)
    orig[["factor_metadata"]] <- md
    orig[["orig_factor_metadata"]] <- md

    avfactors <- md %>% 
                dplyr::filter(type == "factor") %>% 
                select(name, ready, description)
        
    orig[["computed_details"]] <- avfactors 
    data(orig)
  })


  btnStyle <- "color: #fff; background-color: #337ab7; border-color: #2e6da4;"

  loadBtn <- function(id,nm) {
    actionButton(id, nm, icon("paper-plane"), 
      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  }

  output$page = renderUI({
    flog.info("67:loadPage")
    div(
       p(
         if(allDataAvailable()) {
           loadBtn(ns("clearDataBtn"), "Clear Data")
         } else { 
           tagList(actionButton(ns("loadDataBtn"), "Load Phylo Demo Data"))
       }))
  })

  # data
  factorForLoaded <- reactive({
    d <- data()
    if (!is.null(d[["sample"]])) {
      NULL 
    } 
  })


  # allDataAvailable ----
  allDataAvailable <- reactive({
    d <- data()
    if (!is.null(d[["sample"]]) && !is.null(d[["genomic"]])) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })


  # sampleData ----
  sampleDataFromFile <- reactive({
    infile = input$sampleDataFile
    if(is.null(infile)) return(NULL)

    sdata = read.csv(infile$datapath)
    rownames(sdata) <- as.character(sdata$id_full)
    return(sdata)
  })

  sampleData <- reactive({
     if (is.null(data()[["sample"]])) {
       orig <- data()
       sd <- sampleDataFromFile()
       orig[["sample"]] <- sd
       sd
     } else {
       data()[["sample"]]
     }
  })

  # LOAD DEMO DATA BUTTON ----
  ### Load the data then TO DO: change to the analysis tab
  observeEvent(input$loadDataBtn, {
    loadDemoData()
  })

  # CLEAR DEMO DATA BUTTON ----
  observeEvent(input$clearDataBtn, {
    orig <- data()
    orig[["sample"]] <- NULL
    orig[["genomic"]] <- NULL

    orig[["factor_metadata"]] <- NULL
    orig[["orig_factor_metadata"]] <- NULL
    orig[["computed_details"]] <- NULL

    data(orig)
  })

  return(data)
}
