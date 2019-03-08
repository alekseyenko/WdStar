shinyServer(function(input, output, session) {

  data <- callModule(loadPanel, "loadPanel") 

  statusPanel <- callModule(statusPanel, "statusPanel", data)

  factorPanel <- callModule(factorPanel, "factorPanel", data)
  analysisPanel <- callModule(analysisPanel, "analysisPanel", data)
  testsPanel <- callModule(testsPanel, "testsPanel", data)



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
  
  output$diagram <- renderGrViz({
      grViz(spec, quote = TRUE)
      #grViz("digraph test {A; B; A-> B }")
  })

})
