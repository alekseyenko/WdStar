# UI
statusPanelUI <- function(id) {
  ns = NS(id)
  tagList(textOutput(ns("msg")))
}

# Server
statusPanel <- function(input, output, session, data) {

  ns <- session$ns
  data <- data

  output$msg <- renderText({

    sd <- data()[["sample"]]
    gd <- data()[["genomic"]]
    if(is.null(sd) && is.null(gd)) {
      return("No data loaded.")
    }

    sample_msg <- ifelse(is.null(sd),"No sample data", glue("Samples: {nrow(sd)}"))
    g_msg <- ifelse(is.null(gd),"No genomic data", glue("Genomic: {nrow(gd)}"))

    return(str_c(sample_msg, " ", g_msg))
  })

  return(NULL)
}
