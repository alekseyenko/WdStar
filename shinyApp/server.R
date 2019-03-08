shinyAppServer <- function(input, output, session) {

  loadData <- callModule(loadPanel, "loadPanel") 
  
  #genomic_data <- loadData()[["genomic"]]

  # loaded data
  values <- reactiveValues(sampleMainFactor = NULL,
                           sdata = NULL,
                           rcdata = NULL,
                           computed_details = NULL,
                           factor_metadata = NULL,
                           orig_factor_metadata = NULL,
                           orig_factor_data = NULL)

  computed_details <- reactiveVal(NULL)
  sample_details <- reactiveVal(NULL)
  testsDT <- reactiveVal(data.frame())

  # sampleData ----
  sampleDataFromFile <- reactive({
    infile = input$sampleDataFile
    if(is.null(infile)) return(NULL)

    sdata = read.csv(infile$datapath)
    rownames(sdata) <- as.character(sdata$id_full)
    return(sdata)
  })

  sampleData <- reactive({
     if(is.null(values$sdata)) {
       sampleDataFromFile()
     } else {
       values$sdata
     }
 })


  cd <- read_csv("cd.csv")

  factorPanel <- callModule(factorPanel, 
       "factorPanel", computed_details, sampleData, cd)

  load_data <- function() {
     hide(id = "loading-content", anim = TRUE, animType = "fade")
     show("app-content")
  }

  load_data()
  print("AFTER LOAD DATA")

  # rawCountData ----
  rawCountDataFromFile <- reactive({
    infile = input$rawCountFile
    if(is.null(infile)) return(NULL)

    asv = read.table(infile$datapath, row.names = 1, sep = ",", check.names = FALSE, header = TRUE)
    return(t(asv))
  })

  rawCountData <- reactive({
     if(is.null(values$rcdata)) {
       rawCountDataFromFile()
     } else {
       values$rcdata
     }
  })

  # output$physeqData ----
  physeqData <- reactive({
     #req(sampleData(), rawCountData())

     if( !is.null(values$userPhyseq)) {
        #print("return USER physeqData")
        return(values$userPhyseq)
     }

     ##sdata <- sampleData()
     ##asv.t <- rawCountData()
     sdata <- values$sdata
     asv.t <- values$rcdata

     phyloseq( otu_table(asv.t, taxa_are_rows = F), sample_data(sdata))
  })




  # computeTresFtn ----
  computeTresFtn <- function(strata, mainf, dist) {

    physeq <- physeqData()
    dist_matrices = distance(physeq, method = dist)

    #print(glue::glue("computeTresFtn get dist_matrics"))
    # Tw2.test or WdS.test
    #tres = WT.test(dist_matrices,  ...)

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
    #print("distTestOutput")
    df <- as_data_frame(tres$d) %>% rename(dist_cohen_d = value)
  })

  sayecho <- function(tid, strata, mainf, dist) {
    status <- "SUCCESS"
    err <- NULL
    x <- tryCatch({
            computeTresFtn(strata, mainf, dist)
    }, warning = function(w) {
            print(glue::glue("warning: {w}"))
            NULL
    }, error = function(e) {
            err = e
            if("message" %in% names(e)) {
               print(e$message)
               print(glue::glue("error : {e$message} {e$call}"))
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
    print(glue::glue("error : {err}"))
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
    n <- numTests()
    print(str_c("Running tests: ", n, " ", nrow(testsDT())))

    w_df = map_df(1:n, ~ sayecho(testsDT()$test_id[.],
                           testsDT()$strata_factor[.],
                           testsDT()$main_factor[.],
                           testsDT()$distance[.])) 

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
  testsData <- reactive({
     testsDT()
  })

  numTests <- reactive({
    td <- testsData()
    if(is.null(td)) {
       0
    } else {
       nrow(td)
    }
  })

  output$numTests <- reactive({
    numTests()
  })

  output$testTable2 <- renderText({

    req(testsData())
    if(nrow(testsData()) > 0) {
      #datatable(testsData(), options = list(dom='t')) %>%
      t <- testsDT() %>%
          mutate(test=str_c(distance, strata_factor, main_factor, sep="/")) %>%
          select(test)
      str_c(t$test, glue_collapse = ", ")
    }
  })

  # HERE ----

  ### applyFactorNumeric
  applyFactorNumeric <- function(ft, fmethod) {
    flog.info(str_c("applyFactorNumeric: ", ft, " fmethod ", fmethod))
    flog.info(str_c("applyFactorNumeric: ", ft, " fmethod ", fmethod))
    orig_md <- orig_metadata()
    md <- metadata()

    selected_factor <- orig_md %>% pull(ft)
    if (fmethod == "median") {
       midpoint = median(selected_factor, na.rm = TRUE)
    } else {
       midpoint = mean(selected_factor, na.rm = TRUE)
    }
    md[ , ft ] <- orig_md[ , ft] > midpoint
    metadata(md)

    details <- allOriginalFactor()

    mp <- round(midpoint,2)
    new_descript <- str_c("0 = below or at ", fmethod,
        " 1 = above (", mp, ")",sep = "")

    n_df <- data_frame(name=ft,
                       new_ready=TRUE,
                       new_method_applied=fmethod,
                       des=new_descript)

    k_df <- details %>% left_join(n_df, by=c("name"="name"))

    z_df <- k_df %>%
            mutate(ready = ifelse(is.na(new_ready), ready, new_ready)) %>%
            mutate(description = ifelse(is.na(des), description, des)) %>%
            mutate(method_applied =
              ifelse(is.na(new_method_applied),
                    method_applied, new_method_applied)) %>%
            select(-new_ready, -des, -new_method_applied)

    details <- z_df

    cd <- computedDetails() %>%
            mutate(ready, ifelse(name==ft,TRUE,ready)) %>%
            mutate(method_applied, ifelse(name==ft,fmethod,method_applied))  %>%
            mutate(description, ifelse(name==ft,new_descript,description))


    k_df <- computedDetails() %>% left_join(n_df, by=c("name"="name"))

    z_df <- k_df %>%
      mutate(ready = ifelse(is.na(new_ready), ready, new_ready)) %>%
      mutate(description = ifelse(is.na(des), description, des)) %>%
      mutate(method_applied =
               ifelse(is.na(new_method_applied),
                      method_applied, new_method_applied)) %>%
      select(-new_ready, -des, -new_method_applied)

    cd <- z_df

    flog.info("UPDATING applyFactorNumeric computedDetails")
    computedDetails(cd)
    setDetails(details)
  }


  # availableFactors ----
  availableFactors <- reactive({
      sd <- sampleData()
      if(!is.null(sd)) {
        choices <- colnames(sd)
        # DIFFERENCE *****  ignore first column if it is is X
        # ignore first column, which is X
        #choices[1:length(choices)]
        choices
      }
  })

  # allDataAvailable ----
  allDataAvailable <- reactive({
    if (!is.null(sampleData()) && !is.null(rawCountData())) {
      #print("513: ALL AVAILABLE DATA")

      sd <- sampleData()
      values$orig_factor_data = sd
      if(ncol(sd) > 0) {
        md <- map_df(colnames(sd),
                       ~ handle_factor(., dplyr::pull(sd, .)))
        flog.info("allDataAvailable")
        flog.info(md)

        md <- md %>%
                mutate(type = map_chr(colnames(sd), ~ class(sd[[.]])))

        values$factor_metadata = md
        values$orig_factor_metadata = md

        values$computed_details <- data_frame(name = availableFactors(),
                                              ready = FALSE,
                                              description = "init")

        computed_details(values$computed_details)
        flog.info(values$computed_details)
      } else {
        values$factor_metadata = NULL
        values$orig_factor_metadata = NULL
      }
      return(TRUE)
    } else {
      return(FALSE)
    }
  })


  testsPanel <- callModule(testsPanel, "testsPanel", computed_details, 
       sampleData, allDataAvailable, numTests, testsData, testsDT)

  analysisPanel <- callModule(analysisPanel, "analysisPanel", computed_details, 
       sampleData, 
       allDataAvailable, availableFactors, numTests, testsData, testsDT)


  output$dataLoaded <- renderText({
    flog.info("350:server dataLoaded")
    if (allDataAvailable()) {
       flog.info("352:server dataLoaded:data loaded")
       "data loaded"
    } else {""}
  })

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

  output$intro <- renderText({
    markdown::renderMarkdown('intro.md')
  })

}
