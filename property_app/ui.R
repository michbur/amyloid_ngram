library(shiny)
library(reshape2)
load("properties.RData")

props_order <- c(1, 6, 10, 13, 16, 17, 18, 19, 20, #hydroph
                 3, 9, #polarity
                 2, 5, 8, 14, 21:26, #contact
                 4, 7, 11, #beta prob 
                 12, #stability
                 15) #size

shinyUI(fluidPage(
  
  headerPanel("Properties comparision"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("checkProp", label = h3("Choose 1 or more properties"), 
                         choices = as.list(plot_id)[props_order],
                         selected = 1)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Input table", tableOutput("value")),
        tabPanel("Barplot", plotOutput("plot", height = 800)),
        tabPanel("Corplot", plotOutput("corplot", height = 800)),
        tabPanel("PCAplot", plotOutput("pcaplot", height = 800))
      )
    )
  )))