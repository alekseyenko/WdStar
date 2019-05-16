# Shiny app server object ----
shinyAppUI <- 
  tagList(
    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css')),
    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles-dad.css')),
    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'dragula/dragula.min.css')),
    navbarPage("WdStar", id = "tabs", inverse = TRUE,
      header = div(id = "status-panel", statusPanelUI("statusPanel"))
      , tabPanel("Intro", mainPanel(includeMarkdown("intro.md")))
      , tabPanel("Load", loadPanelUI("loadPanel"))
#      , tabPanel("Preprocess", factorPanelUI("factorPanel")) 
      , tabPanel("Analysis", analysisPanelUI("analysisPanel")) 
#      , tabPanel("Test", testsPanelUI("testsPanel")) 
      , tabPanel("About", mainPanel(includeMarkdown("about.md")))
      ))
