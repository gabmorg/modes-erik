# app.R
# script for loading Shiny app for the visualizations

library(shiny)
source("./helpers.R")

ui <-  fluidPage(
  titlePanel("iModEst: Visualizations explaining gene expression regulation across the genome"),
  mainPanel(
    fluidRow( 
      h3("Decomposing predictive accuracy of regulator types per gene across 21 cancers"),
      helpText("Figures 3 and 4 in the paper"),
      h2(textOutput("geneSpecified"))
      # h2("TRPM1"))
    ),
    
    fluidRow(
      textOutput("plotSel"),
      radioButtons("mirnaRadios", "choose a plot type:", 
                   choices = c( "Stacked" = "stack", "Grouped" = "dodge")),
      plotOutput("mirnaPlot", width = "100%", height = "400px", click = NULL,
                 dblclick = NULL, hover = NULL, hoverDelay = NULL,
                 hoverDelayType = NULL, brush = NULL, clickId = NULL,
                 hoverId = NULL, inline = FALSE)
    )
  )
)



server <- function(input, output) {
  
  plotTitle <- reactive({
    plotTitle <- input$geneSelection
    plotTitle
  })
  
  output$geneSpecified <- renderText({
    if(is.null(input$geneSelection)) {
      return("TRPM1")
    }
    else {
      return(plotTitle())
    }
  })

  # input$mirnaRadios
  # An observer is like a reactive expression in that it can read reactive
  # values and call reactive expressions, and will automatically re-execute when
  # those dependencies change. observers use eager evaluation; as soon as their
  # dependencies change, they schedule themselves to re-execute.
  # unlike reactive expressions, it doesn't yield a result and can't be used as
  # an input to other reactive expressions
  
 observeEvent(input$mirnaRadios,
              {output$mirnaPlot <-  renderPlot({
                plot(coefPlot("mirna", "TRPM1", "dodge"))},
                width = "auto",
                height = "auto",
                res = 72,
                env = parent.frame(),
                quoted = FALSE,
                execOnResize = FALSE
              )}
   
  )
 
 # reactive()
  
  # output$plotSel <- renderText({return(plotType())})
  #make reactive
  output$mirnaPlot <-  renderPlot({
    plot(coefPlot("mirna", "TRPM1", "dodge"))},
    width = "auto",
    height = "auto",
    res = 72,
    env = parent.frame(),
    quoted = FALSE,
    execOnResize = FALSE
  )
  
}

shinyApp(ui = ui, server = server)

# to call it from command line: runApp("myapp")


