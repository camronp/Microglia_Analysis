# Microglia Analysis

Interactive Shiny app for exploring the ImageXpress microglia segmentation
and channel intensity quantification macro output (morphology + intensity
CSVs).

## Standalone use

Download `microglia_analysis_app.R`, `data_helpers.R`, and `standalone_app.R`
into the same folder as your macro output CSVs. Open `standalone_app.R` in
RStudio and click **Run App** (or run `shiny::runApp()` from that folder).

## Embedded use

`microglia_analysis_app.R` defines a Shiny module (`microgliaUI()` /
`microgliaServer()`) so it can be embedded as a tab in another app. It's
wired into [PolyLab_Tools](https://github.com/camronp/PolyLab_Tools) as a
git submodule - update this repo as usual, then bump PolyLab_Tools' submodule
reference to pick up the change there.
