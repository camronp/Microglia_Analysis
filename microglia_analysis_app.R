# Microglia Morphology Analysis interactive Shiny companion to analysis_v2.1.Rmd
#
# Run from the same folder as your morphology CSV(s) (the macro's analysis_csv/
# output) via shiny::runApp(). Data loading, filename/cell-name parsing, and
# statistics logic live in data_helpers.R (shared with analysis_v2.1.Rmd) so
# results stay consistent between the static report and this interactive view.
# The first load parses every source CSV and writes merged_morphology_data.csv
# next to them; later loads (including "Reload data" below) reuse that cache
# and are near-instant unless the source CSVs have changed.

required_pkgs <- c("shiny", "DT", "tidyverse", "rstatix", "ggpubr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       ".\nInstall with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


source("data_helpers.R")


df_raw <- load_morphology_data(".")

valid_vars      <- c("project", "well_position", "site", "wavelength", "image")
valid_vars      <- intersect(valid_vars, names(df_raw))
present_metrics <- intersect(all_cell_metrics, names(df_raw))

metric_choices <- setNames(
  lapply(names(metric_categories), function(cat) {
    ms <- intersect(metric_categories[[cat]], present_metrics)
    if (length(ms) == 0) return(NULL)
    setNames(ms, sapply(ms, metric_label))
  }),
  sapply(names(metric_categories), function(cat) str_to_title(gsub("_", " ", cat)))
)
metric_choices <- metric_choices[!sapply(metric_choices, is.null)]


ui <- fluidPage(
  titlePanel("Microglia Morphology"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      uiOutput("data_summary"),
      actionButton("reload_data", "Reload data", width = "100%",
                   title = "Rescan this folder for new/changed CSVs and rebuild the merged data cache"),
      br(), br(),
      selectInput("comparison_variable", "Compare by",
                  choices = setNames(valid_vars, sapply(valid_vars, pretty_var)),
                  selected = if ("well_position" %in% valid_vars) "well_position" else valid_vars[1]),
      fluidRow(
        column(6, actionButton("select_all", "Select all", width = "100%")),
        column(6, actionButton("select_none", "Select none", width = "100%"))
      ),
      br(),
      checkboxGroupInput("selected_levels", "Samples to plot", choices = NULL),
      uiOutput("rename_ui"),
      hr(),
      selectInput("metric", "Metric", choices = metric_choices),
      textInput("plot_title", "Plot title", value = "",
                placeholder = "leave blank to use the metric name"),
      sliderInput("bracket_size", "Bracket/asterisk size", min = 2, max = 8,
                  value = 4, step = 0.5),
      sliderInput("bracket_spacing", "Bracket spacing", min = 0.08, max = 0.4,
                  value = 0.15, step = 0.01),
      uiOutput("bracket_pairs_ui"),
      sliderInput("title_size", "Title size", min = 8, max = 28,
                  value = 16, step = 1),
      sliderInput("label_size", "Axis label size", min = 8, max = 24,
                  value = 15, step = 1),
      hr(),
      downloadButton("download_plot", "Download plot (PNG)", width = "100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Plot",
                 br(),
                 DTOutput("counts_table"),
                 br(),
                 plotOutput("morph_plot")
        ),
        tabPanel("Statistics",
                 br(),
                 p("Every pair of selected samples compared directly, BH-corrected across pairs; ",
                   "the same test is used for the significance brackets on the Plot tab (shown ",
                   "there only when significant, and only for ", strong("≤ 8"),
                   " samples selected to stay readable)."),
                 p(em("Why this test: ", strong("Welch's t-test"), " is used for a metric only if ",
                      "every selected sample's values pass a Shapiro-Wilk normality test (p > 0.05) ",
                      "with nonzero variance; otherwise the non-parametric ", strong("Mann-Whitney U"),
                      " test is used instead, since it doesn't assume normally-distributed data.")),
                 fluidRow(
                   column(6, selectizeInput("stats_metric_filter", "Filter by metric",
                                             choices = metric_choices, selected = character(0),
                                             multiple = TRUE,
                                             options = list(plugins = list("remove_button")))),
                   column(3, style = "margin-top: 25px;",
                          actionButton("stats_metric_all", "All metrics")),
                   column(3, style = "margin-top: 25px;",
                          downloadButton("download_pairwise", "Download table (CSV)"))
                 ),
                 DTOutput("pairwise_table")
        )
      )
    )
  )
)


server <- function(input, output, session) {

  df_store <- reactiveVal(df_raw)

  observeEvent(input$reload_data, {
    showNotification("Rescanning CSVs and rebuilding the merged data cache...",
                      id = "reload", duration = NULL, type = "message")
    df_store(load_morphology_data(".", force_rebuild = TRUE))
    removeNotification("reload")
  })

  output$data_summary <- renderUI({
    df <- df_store()
    helpText(strong(n_distinct(df$image_id)), "images,",
             strong(format(nrow(df), big.mark = ",")), "cells loaded from",
             strong(n_distinct(df$source_file)), "CSV file(s) in this folder.")
  })

  levels_all <- reactive({
    req(input$comparison_variable)
    sort(unique(na.omit(df_store()[[input$comparison_variable]])))
  })

  observeEvent(input$comparison_variable, {
    lv <- levels_all()
    updateCheckboxGroupInput(session, "selected_levels", choices = lv, selected = lv)
  })

  observeEvent(input$select_all,  updateCheckboxGroupInput(session, "selected_levels",
                                                             choices = levels_all(), selected = levels_all()))
  observeEvent(input$select_none, updateCheckboxGroupInput(session, "selected_levels",
                                                             choices = levels_all(), selected = character(0)))

 
  selected_sorted <- reactive(sort(input$selected_levels))

  output$rename_ui <- renderUI({
    lv <- selected_sorted()
    if (length(lv) == 0) return(NULL)
    tagList(
      tags$strong("Rename samples (optional)"),
      lapply(seq_along(lv), function(i) {
        textInput(paste0("rename_", i), label = lv[i], value = lv[i])
      })
    )
  })

  rename_map <- reactive({
    lv <- selected_sorted()
    req(length(lv) > 0)
    vals <- vapply(seq_along(lv), function(i) {
      v <- input[[paste0("rename_", i)]]
      if (is.null(v) || v == "") lv[i] else v
    }, character(1))
    setNames(vals, lv)
  })

  filtered_df <- reactive({
    req(length(input$selected_levels) > 0)
    cv <- input$comparison_variable
    rm <- rename_map()
    df_store() %>%
      filter(.data[[cv]] %in% names(rm)) %>%
      mutate(group = factor(unname(rm[as.character(.data[[cv]])]), levels = unique(rm)))
  })

  plot_title <- reactive({
    if (nzchar(trimws(input$plot_title))) input$plot_title else metric_label(input$metric)
  })

  # Which pairwise comparisons are significant for the current metric/groups -
  # the candidate set the "brackets to show" checkboxes are built from.
  sig_pairs <- reactive({
    fd <- filtered_df()
    req(input$metric)
    if (n_distinct(fd$group) < 2 || n_distinct(fd$group) > 8) return(NULL)
    pw <- compare_pairs_metric(fd, input$metric, "group")
    if (is.null(pw)) return(NULL)
    sig <- pw %>% filter(p.adj < 0.05)
    if (nrow(sig) == 0) return(NULL)
    sig
  })

  output$bracket_pairs_ui <- renderUI({
    sig <- sig_pairs()
    if (is.null(sig)) return(NULL)
    labels <- pair_labels(sig)
    tagList(
      tags$strong("Significance brackets to show"),
      lapply(seq_along(labels), function(i) {
        checkboxInput(paste0("bracket_show_", i), label = labels[i], value = TRUE)
      })
    )
  })

  # Individually-named checkboxes (rather than one checkboxGroupInput) so an
  # unchecked box is unambiguous - checkboxGroupInput reports NULL both when
  # the client hasn't echoed its initial value yet and when everything's been
  # unchecked, which briefly hid every bracket on first load.
  bracket_pairs_selected <- reactive({
    sig <- sig_pairs()
    if (is.null(sig)) return(NULL)
    labels <- pair_labels(sig)
    keep <- vapply(seq_along(labels), function(i) {
      v <- input[[paste0("bracket_show_", i)]]
      if (is.null(v)) TRUE else v
    }, logical(1))
    labels[keep]
  })

  # A fixed device height can't fit an arbitrary number of stacked
  # significance brackets without them overlapping - the panel has to grow
  # with the stack. bracket_plot_height() (data_helpers.R) returns however
  # many inches are needed; 96px/in converts that for renderPlot's height.
  plot_height_px <- reactive({
    fd <- filtered_df()
    req(nrow(fd) > 0)
    bracket_plot_height(fd, input$metric, "group", bracket_size = input$bracket_size,
                         bracket_spacing = input$bracket_spacing,
                         selected_pairs = bracket_pairs_selected()) * 96
  })

  output$morph_plot <- renderPlot({
    fd <- filtered_df()
    validate(need(nrow(fd) > 0, "No data for the current sample selection."))
    plot_single_metric_pub(fd, input$metric, group_col = "group",
                            group_label = pretty_var(input$comparison_variable),
                            title = plot_title(), bracket_size = input$bracket_size,
                            bracket_spacing = input$bracket_spacing,
                            title_size = input$title_size, label_size = input$label_size,
                            selected_pairs = bracket_pairs_selected())
  }, height = function() plot_height_px())

  output$counts_table <- renderDT({
    fd <- filtered_df()
    fd %>% count(group, name = "n_cells") %>%
      left_join(fd %>% distinct(group, image_id) %>% count(group, name = "n_images"), by = "group") %>%
      rename(!!pretty_var(input$comparison_variable) := group) %>%
      datatable(rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })

  observeEvent(input$stats_metric_all, {
    updateSelectizeInput(session, "stats_metric_filter", selected = present_metrics)
  })

  pairwise_results <- reactive({
    fd <- filtered_df()
    validate(need(n_distinct(fd$group) >= 2, "Select at least 2 samples to run statistics."))
    metrics <- intersect(input$stats_metric_filter, present_metrics)
    validate(need(length(metrics) > 0, "Select at least 1 metric to filter the table."))
    res <- compare_pairs_all_metrics(fd, metrics, "group")
    validate(need(nrow(res) > 0, "Not enough data to compare the selected samples."))
    res
  })

  output$pairwise_table <- renderDT({
    pairwise_results() %>%
      mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
      datatable(rownames = FALSE)
  })

  output$download_pairwise <- downloadHandler(
    filename = function() paste0("pairwise_comparisons_", Sys.Date(), ".csv"),
    content = function(file) {
      readr::write_csv(pairwise_results(), file)
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() paste0("morphology_", input$metric, "_", Sys.Date(), ".png"),
    content = function(file) {
      fd <- filtered_df()
      p <- plot_single_metric_pub(fd, input$metric, group_col = "group",
                                   group_label = pretty_var(input$comparison_variable),
                                   title = plot_title(), bracket_size = input$bracket_size,
                                   bracket_spacing = input$bracket_spacing,
                                   title_size = input$title_size, label_size = input$label_size,
                                   selected_pairs = bracket_pairs_selected())
      h <- bracket_plot_height(fd, input$metric, "group", bracket_size = input$bracket_size,
                                bracket_spacing = input$bracket_spacing,
                                selected_pairs = bracket_pairs_selected())
      ggsave(file, p, width = max(7, 1.2 + 1.1 * n_distinct(fd$group)), height = h, dpi = 300)
    }
  )
}

shinyApp(ui, server)

