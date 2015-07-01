library(shiny)
library(ggplot2)
library(Cairo)
load("properties.RData")

shinyServer(function(input, output) {
  
  output$value <- renderTable({ 
    tab <- t(plot_values[, as.numeric(input$checkProp)])
    rownames(tab) <- names(plot_id)[as.numeric(input$checkProp)]
    tab
    })
  
  output$plot <- renderPlot({
    ggplot(mplot_values[mplot_values[["id"]] %in% input$checkProp, ], aes(x = id, fill = id, y = value)) +
      geom_bar(stat="identity", position = "dodge") + 
      facet_wrap(~ aa, ncol = 5) +
      scale_x_discrete("Amino acid") + 
      scale_y_continuous("Normalized value") +
      scale_fill_discrete(name="Property name", 
                          labels = names(plot_id)[as.numeric(input$checkProp)]) +
      guides(fill=guide_legend(ncol=2)) + 
      theme(plot.background=element_rect(fill = "transparent",
                                         colour = "transparent"),
            panel.grid.major = element_line(colour="grey", linetype = "dashed", size = 0.5),
            panel.grid.major = element_line(colour="lightgrey", linetype = "dashed", size = 0.5),
            panel.background = element_rect(fill = "transparent",colour = "black"),
            legend.background = element_rect(fill = "NA"),
            legend.position = "bottom",
            strip.background = element_rect(fill = "NA", colour = "NA"))
  })
  
})


