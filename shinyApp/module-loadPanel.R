library(dplyr)
library(shiny)
library(glue)
library(DT)
library(futile.logger)

# UI
loadPanelUI <- function(id) {
  ns = NS(id)

  tagList(h3("loadPanel"), 
          actionButton("loadDataBtn", "Load Demo Data Btn"),
          uiOutput(ns("loadMenu")),
          uiOutput(ns("load")))
}

# Server
loadPanel <- function(input, output, session) {

  flog.info("20: loadPanel")
  data <- reactiveVal(list(genomic=NULL, sample=NULL))

  loadDemoData2 <- reactive({
    flog.info("24:loadDemoData2")

    orig <- data()
    sdata = read.csv("data/sample_data.csv")
    rownames(sdata) <- as.character(sdata$id_full)

    asv = read.table("data/raw_count.csv", row.names = 1, sep = ",", check.names = FALSE, header = TRUE)
    rcdata = t(asv)

    orig <- data()
    orig[["sample"]] <- sdata
    orig[["genomic"]] <- rcdata

    data(orig)
  })

  loadDemoData1 <- reactive({
    flog.info("41:loadDemoData1")
  })


  loadDemoData <- reactive({
    flog.info("46:loadDemoData")

    sdata <- physeqToSample(physeq)
    rcdata <- physeqToRawCount(physeq)

    orig <- data()
    orig[["sample"]] <- sdata
    orig[["genomic"]] <- rcdata

    data(orig)
  })


  btnStyle <- "color: #fff; background-color: #337ab7; border-color: #2e6da4;"

  loadBtn <- function(id,nm) {
    actionButton(id, nm, icon("paper-plane"), 
      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  }

  output$loadMenu = renderUI({
    flog.info("67:loadMenu")
     div(
         p(
           #menuItem('Raw Count/CSV', tabName = 'data'),
         if(allDataAvailable()) {
            loadBtn("clearDataBtn", "Clear Data")
         } else { 
            tagList(p(
            loadBtn("xx", "xx"),
            loadBtn("loadDemoDataBtn", "Load Demo Data"),
            loadBtn("loadDemoDataBtn2", "Load Demo Data 2")))
         }))})



  # allDataAvailable ----
  allDataAvailable <- reactive({
    d <- data()
    if (!is.null(d[["sample"]]) && !is.null(d[["genomic"]])) {
      print("513: ALL AVAILABLE DATA")
      flog.info("513: ALL AVAILABLE DATA")

      #????
      sd <- sampleData()
      d[["orig_factor_data"]] <- sd
      if(ncol(sd) > 0) {
        md <- map_df(colnames(sd),
                       ~ handle_factor(., dplyr::pull(sd, .)))
        flog.info("allDataAvailable")
        flog.info(md)

        md <- md %>%
                mutate(type = map_chr(colnames(sd), ~ class(sd[[.]])))

        d[["factor_metadata"]] <- md
        d[["orig_factor_metadata"]] <- md

        d[["computed_details"]] <- data_frame(name = availableFactors(),
                                              ready = FALSE,
                                              description = "init")

        computed_details(d[["computed_details"]])
        flog.info(d[["computed_details"]])
      } else {
        values$factor_metadata = NULL
        values$orig_factor_metadata = NULL
        d[["factor_metadata"]] <- NULL
        d[["orig_factor_metadata"]] <- NULL
      }
      data(d)
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
  ### Load the data then change to the analysis tab
  observeEvent(input$loadDataBtn, {
    flog.info("150:loadDataBtn")

    #show(id = "updating-content", anim = TRUE, animType = "fade")
    loadDemoData()

    #hide(id = "updating-content", anim = TRUE, animType = "fade")
    #updateTabsetPanel(session, inputId = "sidebar", selected = "analysis")
  })

  observeEvent(input$loadDemoDataBtn, {
    flog.info("159:loadDemoDataBtn")
    print("LOAD")
  })

  observeEvent(input$xx, {
    flog.info("166:xx loadDemoDataBtn")
    print("LOAD")
  })

  observeEvent(input$loadDemoDataBtn2, {
    flog.info("163:loadDemoDataBtn2")
  })


  # CLEAR DEMO DATA BUTTON ----
  observeEvent(input$clearDataBtn, {
    flog.info("169: clearDataBtn")

    orig <- data()
    orig[["sample"]] <- NULL
    orig[["genomic"]] <- NULL

    data(orig)
  })

  observe({
          print("181: loadPanel HELLO")
          flog.info("181: loadPanel HELLO")
  })


  return(data)
}
