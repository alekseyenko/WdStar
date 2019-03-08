distance_choices <- phyloseq::distanceMethodList

# input selection panel
input_selection_panel <-
  div(box(fileInput("rawCountFile", "Raw Count File")),
      box(fileInput("sampleDataFile", "Sample Data File")))

sidebar <-  dashboardSidebar(
  sidebarMenu(id="sidebar",
    menuItem("Introduction", tabName = "intro", icon = icon("home")),
    #menuItem('Load data', icon = icon("database"), tabName = 'data', uiOutput("loadMenu")),
    menuItem('Load', icon = icon("database"), tabName = 'load')
#    , menuItem('Factor', icon = icon("database"), tabName = 'factor'),
#    menuItem("Analysis", icon = icon("tasks"), tabName = "analysis" ),
#    menuItem('Tests', icon = icon("database"), tabName = 'tests')
    ))


# Shiny app server object ----
shinyAppUI <- dashboardPage(

  dashboardHeader(title = "Shiny WdStar"),
  dashboardSidebar(sidebar),
  dashboardBody(id = "dashboard",
    shinyjs::useShinyjs(),
    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css')),

    div(id = "loading-content", h1("LOADING...")),

    # some divs to display to show in-progress
    hidden(div(id = "updating-content", p("UPDATING...."))),
    hidden(div(id = "running-tests", h1("RUNNING TESTS...."))),

    tabItems(
      tabItem(tabName = "intro", htmlOutput("intro")),
      #tabItem(tabName = "data", input_selection_panel),
      tabItem(tabName = "load", loadPanelUI("loadPanel"))
#      , tabItem(tabName = "factor", factorPanelUI("factorPanel")),
#      tabItem(tabName = "tests", testsPanelUI("testsPanel")),
#      tabItem(tabName = "analysis", analysisPanelUI("analysisPanel"))
      )))
