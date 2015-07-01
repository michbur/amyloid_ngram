library(shiny)

load("properties.RData")


shinyUI(fluidPage(
  
  headerPanel("Properties comparision"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("checkProp", label = h3("Choose 1 or more properties"), 
                         choices = as.list(plot_id),
                         selected = c(1,2))
    ),
    
    mainPanel(
      plotOutput("plot", height = 800),
      tableOutput("value")
    )
  )))