library(shiny)
library(ggplot2)
load("properties.RData")
my_theme <- theme(plot.background=element_rect(fill = "transparent",
                                               colour = "transparent"),
                  panel.grid.major = element_line(colour="grey", linetype = "dashed", size = 0.5),
                  panel.grid.major = element_line(colour="lightgrey", linetype = "dashed", size = 0.5),
                  panel.background = element_rect(fill = "transparent",colour = "black"),
                  legend.background = element_rect(fill = "NA"),
                  legend.position = "bottom",
                  strip.background = element_rect(fill = "NA", colour = "NA"))

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
      guides(fill=guide_legend(ncol=2)) + my_theme
      
  })
  
  output$corplot <- renderPlot({
    nms <- names(plot_id)[as.numeric(input$checkProp)]
    br_nms <- sapply(nms, function(single_nm) {
      len <- nchar(single_nm)
      space_pos <-  which(strsplit(single_nm, "")[[1]] == " ")
      break_pos <- space_pos[which.min(abs(len/2 - space_pos))]
      paste0(substr(single_nm, 1, break_pos - 1), "\n", substr(single_nm, break_pos + 1, len))
    })
    
    corm <- melt(cor(plot_values[, as.numeric(input$checkProp)] - 0.5))
    corm[["Var1"]] <- factor(corm[["Var1"]], labels = br_nms)
    corm[["Var2"]] <- factor(corm[["Var2"]], labels = br_nms)
    
    ggplot(data = corm, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_fill_gradient2(low = "white", high = "red", mid = "blue", 
                           midpoint = 0, limit = c(-1,1), name="Correlation\ncoefficient") +
      my_theme
  })
  
  output$pcaplot <- renderPlot({
    pca_res <- princomp(plot_values[, as.numeric(input$checkProp)])
    ggdat <- cbind(data.frame(pca_res[["scores"]]), aa = rownames(pca_res[["scores"]]))
    ggplot(ggdat, aes(x = Comp.1, y = Comp.2, colour = Comp.3, label = aa)) +
      geom_point(size = 5) +
      scale_colour_gradient2(low = "red", high = "yellow", mid = "orange") +
      my_theme +
      geom_text(hjust=1.8, vjust=1.8)
  })
  
})


