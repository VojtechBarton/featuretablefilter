# featuretablefilter Shiny Dashboard
# Interactive filtering with real-time preview

library(shiny)
library(ggplot2)
library(DT)

# Source UI and Server components
source(file.path(dirname(__FILE__), "ui.R"))
source(file.path(dirname(__FILE__), "server.R"))

# Run the application
shinyApp(ui = ui, server = server)
