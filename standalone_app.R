# Standalone launcher for local development/testing of the microglia app on
# its own, outside the combined PolyLab_Tools shell.
#
# Copy this file alongside microglia_analysis_app.R and data_helpers.R into
# the same folder as your macro output CSVs (or set that folder as the R
# working directory), then open this file in RStudio and click "Run App"
# (or run shiny::runApp() from that folder).

source("microglia_analysis_app.R")

shinyApp(
  ui = microgliaUI("microglia"),
  server = function(input, output, session) microgliaServer("microglia")
)
