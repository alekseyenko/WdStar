# Shiny app server object ----
shinyAppUI <- dashboardPage(

  dashboardHeader(title = "Shiny WdStar"),
  dashboardSidebar(sidebarMenu(id="sidebar",
    menuItem("Introduction", tabName = "intro", icon = icon("home")),
    menuItem('Load', icon = icon("database"), tabName = 'load'),
    menuItem('Preprocess', icon = icon("database"), tabName = 'factor')
    , menuItem('Analysis', icon = icon("database"), tabName = 'analysis')
    )),
  dashboardBody(id = "dashboard",
    shinyjs::useShinyjs(),
    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css')),

    div(id = "loading-content", h1("LOADING...")),

    # some divs to display to show in-progress
    hidden(div(id = "updating-content", p("UPDATING...."))),
    hidden(div(id = "running-tests", h1("RUNNING TESTS...."))),

    tabItems(
      tabItem(tabName = "intro", tagList(includeMarkdown("intro.md")))
      , tabItem(tabName = "load", loadPanelUI("loadPanel"))
      , tabItem(tabName = "factor", factorPanelUI("factorPanel"))
      , tabItem(tabName = "analysis", analysisPanelUI("analysisPanel"))
      )))
