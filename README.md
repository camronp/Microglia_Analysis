# Microglia Analysis

Interactive Shiny app for exploring the ImageXpress microglia segmentation
and channel intensity quantification macro output (morphology + intensity
CSVs). Each visitor uploads their own CSVs through the browser - nothing is
read from disk at startup.

## Standalone use

Download `microglia_analysis_app.R`, `data_helpers.R`, and `standalone_app.R`
into any folder. Open `standalone_app.R` in RStudio and click **Run App**
(or run `shiny::runApp()` from that folder), then upload your CSVs once
it's running.

## Embedded use

`microglia_analysis_app.R` defines a Shiny module (`microgliaUI()` /
`microgliaServer()`) so it can be embedded as a tab in another app. It's
wired into [PolyLab_Tools](https://github.com/camronp/PolyLab_Tools) as
plain copied files (`microglia_app/`) - update this repo as usual, then
copy the changed file(s) over there to pick up the change.
