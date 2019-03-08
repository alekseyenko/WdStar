library(dplyr)
library(shiny)
library(glue)
library(DT)
library(futile.logger)

# UI
loadPanelUI <- function(id) {
  ns = NS(id)

  tagList(includeMarkdown("load.md"),
          uiOutput(ns("load")))
}

# Server
loadPanel <- function(input, output, session) {

  data <- reactiveVal(list(genomic = NULL, 
                           sample = NULL, 
                           orig_sample = NULL, 
                           factor_metadata = NULL,
                           orig_factor_metadata = NULL))

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

  loadDemoData2 <- reactive({
    flog.info("24:loadDemoData2")

    orig <- data()
    sdata = read.csv("data/sample_data.csv")
    rownames(sdata) <- as.character(sdata$id_full)

    asv = read.table("data/raw_count.csv", row.names = 1, sep = ",", check.names = FALSE, header = TRUE)
    rcdata = t(asv)

    orig <- data()
    orig[["sample"]] <- sdata
    orig[["orig_sample"]] <- sdata
    orig[["genomic"]] <- rcdata

    data(orig)
    init_factors()
  })

  loadDemoData1 <- reactive({
    flog.info("41:loadDemoData1")
  })


  ## process factor data ----
  init_factors <- reactive({
    flog.info("61:init_factors ")
    d <- data()
    sd <- d[["sample"]]
    d[["orig_factor_data"]] <- sd

    # compute factor meta dta
    if(ncol(sd) > 0) {
      md <- map_df(colnames(sd), ~ handle_factor(., dplyr::pull(sd, .))) 
    
      d[["factor_metadata"]] <- md
      d[["orig_factor_metadata"]] <- md

    } else {
      d[["factor_metadata"]] <- NULL
      d[["orig_factor_metadata"]] <- NULL
    }
    data(d)
    flog.info("85: EXIT init_factors")
  })

  loadDemoData <- reactive({
    flog.info("46:loadDemoData")

    orig <- data()
    sdata <- physeqToSample(physeq)
    orig[["sample"]] <- sdata
    orig[["orig_sample"]] <- sdata
    orig[["genomic"]] <- physeqToRawCount(physeq)

    data(orig)
    init_factors()
  })


  btnStyle <- "color: #fff; background-color: #337ab7; border-color: #2e6da4;"

  loadBtn <- function(id,nm) {
    actionButton(id, nm, icon("paper-plane"), 
      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  }

  ## load ----
  output$load = renderUI({
    flog.info("67:load")
    ns <- session$ns
    div(
         p(
           #menuItem('Raw Count/CSV', tabName = 'data'),
         if(dataLoaded()) {
            loadBtn(ns("clearDataBtn"), "Clear Data")
         } else { 
            tagList(p(
            loadBtn(ns("loadDemoDataBtn"), "Load Demo Data"),
            loadBtn(ns("loadDemoDataBtn2"), "Load Demo Data 2")))
         }))})


  ## dataLoaded ----
  dataLoaded <- reactive({
    return(!is.null(data()[["sample"]]) && !is.null(data()[["genomic"]]))
  })


  ## sampleDataFromFile ----
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

  ## loadDemoDataBtn ---- 
  observeEvent(input$loadDemoDataBtn, {
    flog.info("159:loadDemoDataBtn")
    print("164: LOAD Demo Data Btn")
    loadDemoData()
  })

  observeEvent(input$loadDemoDataBtn2, {
    flog.info("163:loadDemoDataBtn2")
    loadDemoData2()
  })

  ## clearDataBtn ----
  observeEvent(input$clearDataBtn, {
    flog.info("169: clearDataBtn")

    orig <- data()
    orig[["orig_sample"]] <- NULL
    orig[["sample"]] <- NULL
    orig[["genomic"]] <- NULL

    data(orig)
  })

  return(data)
}
