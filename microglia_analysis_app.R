# Interactive Microglia Analysis - Shiny module
#
# Defines microgliaUI(id) / microgliaServer(id) so this app can be embedded
# as a tab in another Shiny app (see PolyLab_Tools). To run it on its own,
# use standalone_app.R in this same folder.
#
# The caller must source data_helpers.R (from wherever this file's copy of
# it lives) before sourcing this file - not done internally here, since this
# file may be sourced from a different working directory than the one it
# lives in (e.g. when embedded from another folder), while data_helpers.R
# always sits right next to this file.
#
# Each visitor uploads their own morphology CSV(s) (and optionally a channel
# intensity CSV) through the browser - nothing is read from disk at startup,
# so a deployed copy of this app doesn't ship with (or share) any one
# person's data. Uploaded files are copied into a per-session temp folder,
# preserving their original filenames, so the existing filename-driven
# parsing in load_morphology_data()/load_intensity_data() (macro-merged CSV
# vs. legacy per-image CSVs, filename metadata extraction, etc.) keeps
# working unchanged.


required_pkgs <- c("shiny", "DT", "tidyverse", "rstatix", "ggpubr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       ".\nInstall with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


intensity_stat_choices <- c(
  "Mean intensity"          = "mean_intensity",
  "Median intensity"        = "median_intensity",
  "Min intensity"           = "min_intensity",
  "Max intensity"           = "max_intensity",
  "Std. dev. intensity"     = "stddev_intensity",
  "Integrated density"      = "integrated_density",
  "Raw integrated density"  = "raw_integrated_density",
  "Area (px)"               = "area_px"
)


microgliaUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Everything Microglia"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput(ns("morph_files"), "Morphology CSV(s)", multiple = TRUE, accept = ".csv",
                  placeholder = "merged_morphology_data_*.csv, or morphology_results_*.csv"),
        fileInput(ns("intensity_file"), "Channel intensity CSV (optional)", accept = ".csv",
                  placeholder = "channel_intensity_results_*.csv"),
        uiOutput(ns("data_summary")),
        uiOutput(ns("comparison_variable_ui")),
        fluidRow(
          column(6, actionButton(ns("select_all"), "Select all", width = "100%")),
          column(6, actionButton(ns("select_none"), "Select none", width = "100%"))
        ),
        br(),
        checkboxGroupInput(ns("selected_levels"), "Samples to plot", choices = NULL, inline = TRUE),
        uiOutput(ns("rename_ui")),
        hr(),
        uiOutput(ns("metric_ui")),
        textInput(ns("plot_title"), "Plot title", value = "",
                  placeholder = "leave blank to use the metric name"),
        sliderInput(ns("bracket_size"), "Bracket/asterisk size", min = 2, max = 8,
                    value = 6, step = 0.5),
        sliderInput(ns("bracket_spacing"), "Bracket spacing", min = 0.08, max = 0.4,
                    value = 0.08, step = 0.01),
        uiOutput(ns("bracket_pairs_ui")),
        sliderInput(ns("title_size"), "Title size", min = 8, max = 28,
                    value = 16, step = 1),
        radioButtons(ns("title_justify"), "Title alignment",
                     choices = c("Left" = "left", "Center" = "center"),
                     selected = "left", inline = TRUE),
        sliderInput(ns("label_size"), "Axis label size", min = 8, max = 24,
                    value = 15, step = 1),
        checkboxInput(ns("log_scale"), "Log-scale y-axis (for extreme outliers)", value = FALSE),
        hr(),
        downloadButton(ns("download_plot"), "Download plot (PNG)", width = "100%")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Morphology",
                   br(),
                   fluidRow(
                     column(4, selectInput(ns("morph_palette_choice"), "Color palette",
                                            choices = palette_choices, selected = "okabe")),
                     column(4, conditionalPanel(
                       condition = sprintf("input['%s'] == 'custom'", ns("morph_palette_choice")),
                       textInput(ns("morph_custom_palette"), "Custom colors",
                                 placeholder = "Comma-separated hex codes, e.g. #FF0000, #00A651")
                     ))
                   ),
                   uiOutput(ns("morph_plot_ui")),
                   br(),
                   fluidRow(
                     column(6, downloadButton(ns("download_pairwise"), "Download this metric's table (CSV)")),
                     column(6, downloadButton(ns("download_pairwise_all"), "Download all metrics (CSV)"))
                   ),
                   br(),
                   DTOutput(ns("pairwise_table")),
                   br(),
                   DTOutput(ns("counts_table"))
          ),
          tabPanel("Channel Intensity",
                   br(),
                   uiOutput(ns("intensity_tab_ui"))
          )
        )
      )
    )
  )
}


microgliaServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- per-session upload staging -----------------------------------
    # Uploaded files land at random tmp paths (input$x$datapath); copying
    # them into a per-session folder under their real name lets
    # load_morphology_data()/load_intensity_data() keep matching filenames
    # by pattern exactly as they do when reading a macro output folder
    # directly. Kept in separate subfolders so re-uploading one doesn't
    # require touching (or accidentally stacking up stale copies of) the
    # other.
    session_dir  <- file.path(tempdir(), paste0("microglia_", session$token))
    morph_dir    <- file.path(session_dir, "morphology")
    intensity_dir <- file.path(session_dir, "intensity")
    dir.create(morph_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(intensity_dir, recursive = TRUE, showWarnings = FALSE)
    session$onSessionEnded(function() unlink(session_dir, recursive = TRUE))

    df_store <- reactiveVal(NULL)
    intensity_store <- reactiveVal(tibble())

    stage_upload <- function(fd, dir) {
      old <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
      if (length(old) > 0) unlink(old)
      file.copy(fd$datapath, file.path(dir, fd$name), overwrite = TRUE)
    }

    observeEvent(input$morph_files, {
      stage_upload(input$morph_files, morph_dir)
      new_morph <- tryCatch(load_morphology_data(morph_dir, force_rebuild = TRUE),
                             error = function(e) {
                               showNotification(paste("Couldn't read morphology CSV(s):", conditionMessage(e)),
                                                 type = "error", duration = NULL)
                               NULL
                             })
      if (!is.null(new_morph)) {
        df_store(new_morph)
        idf <- isolate(intensity_store())
        if (nrow(idf) > 0) intensity_store(attach_cell_metadata(idf, new_morph))
      }
    })

    observeEvent(input$intensity_file, {
      stage_upload(input$intensity_file, intensity_dir)
      new_intensity <- tryCatch(load_intensity_data(intensity_dir),
                                 error = function(e) {
                                   showNotification(paste("Couldn't read intensity CSV:", conditionMessage(e)),
                                                     type = "error", duration = NULL)
                                   NULL
                                 })
      df <- isolate(df_store())
      if (!is.null(new_intensity) && !is.null(df)) {
        intensity_store(attach_cell_metadata(new_intensity, df))
      } else if (!is.null(new_intensity)) {
        intensity_store(new_intensity)
      }
    })

    valid_vars_r <- reactive({
      df <- df_store(); req(df)
      intersect(c("project", "well_position", "site", "wavelength", "image"), names(df))
    })

    present_metrics_r <- reactive({
      df <- df_store(); req(df)
      intersect(all_cell_metrics, names(df))
    })

    metric_choices_r <- reactive({
      pm <- present_metrics_r()
      choices <- setNames(
        lapply(names(metric_categories), function(cat) {
          ms <- intersect(metric_categories[[cat]], pm)
          if (length(ms) == 0) return(NULL)
          setNames(ms, sapply(ms, metric_label))
        }),
        sapply(names(metric_categories), function(cat) str_to_title(gsub("_", " ", cat)))
      )
      choices[!sapply(choices, is.null)]
    })

    output$data_summary <- renderUI({
      df <- df_store()
      if (is.null(df)) {
        return(helpText("Upload your morphology CSV(s) above to begin",
                         "(optionally add the channel intensity CSV too)."))
      }
      idf <- intensity_store()
      helpText(
        strong(n_distinct(df$image_id)), "images,",
        strong(format(nrow(df), big.mark = ",")), "cells loaded from",
        strong(n_distinct(df$source_file)), "CSV file(s).",
        if (nrow(idf) > 0) {
          tagList(br(), strong(format(nrow(idf), big.mark = ",")), "intensity measurements across",
                  strong(n_distinct(idf$channel)), "channel(s).",
                  if ("seeded_status" %in% names(idf)) {
                    counts <- idf %>% distinct(cell_name, seeded_status) %>% count(seeded_status)
                    tagList(br(), paste(paste0(str_to_title(counts$seeded_status), ": ",
                                                format(counts$n, big.mark = ",")), collapse = ", "),
                            "cells.")
                  })
        }
      )
    })

    output$comparison_variable_ui <- renderUI({
      vv <- valid_vars_r()
      req(length(vv) > 0)
      selectInput(ns("comparison_variable"), "Compare by",
                  choices = setNames(vv, sapply(vv, pretty_var)),
                  selected = if ("well_position" %in% vv) "well_position" else vv[1])
    })

    output$metric_ui <- renderUI({
      selectInput(ns("metric"), "Metric", choices = metric_choices_r())
    })

    levels_all <- reactive({
      req(input$comparison_variable)
      sort(unique(na.omit(df_store()[[input$comparison_variable]])))
    })

    observeEvent(input$comparison_variable, {
      lv <- levels_all()
      updateCheckboxGroupInput(session, "selected_levels", choices = lv, selected = lv, inline = TRUE)
    })

    observeEvent(input$select_all,  updateCheckboxGroupInput(session, "selected_levels",
                                                               choices = levels_all(), selected = levels_all(),
                                                               inline = TRUE))
    observeEvent(input$select_none, updateCheckboxGroupInput(session, "selected_levels",
                                                               choices = levels_all(), selected = character(0),
                                                               inline = TRUE))


    selected_sorted <- reactive(sort(input$selected_levels))

    output$rename_ui <- renderUI({
      lv <- selected_sorted()
      if (length(lv) == 0) return(NULL)
      tagList(
        tags$strong("Rename samples (optional)"),
        lapply(seq_along(lv), function(i) {
          textInput(ns(paste0("rename_", i)), label = lv[i], value = lv[i])
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

    apply_group_selection <- function(df, cv, rm) {
      df %>%
        filter(.data[[cv]] %in% names(rm)) %>%
        mutate(group = factor(unname(rm[as.character(.data[[cv]])]), levels = unique(rm)))
    }

    filtered_df <- reactive({
      req(df_store(), input$comparison_variable, length(input$selected_levels) > 0)
      apply_group_selection(df_store(), input$comparison_variable, rename_map())
    })

    intensity_filtered_df <- reactive({
      req(length(input$selected_levels) > 0, input$intensity_channel)
      idf <- intensity_store() %>% filter(channel == input$intensity_channel)
      seeded_filter <- input$intensity_seeded_filter %||% "all"
      if ("seeded_status" %in% names(idf) && seeded_filter != "all") {
        idf <- idf %>% filter(seeded_status == seeded_filter)
      }
      apply_group_selection(idf, input$comparison_variable, rename_map())
    })

    plot_title <- reactive({
      if (nzchar(trimws(input$plot_title))) input$plot_title else metric_label(input$metric)
    })

    # Changing each tab to have its own palette dropdown
    morph_current_palette <- reactive({
      custom_cols <- if (identical(input$morph_palette_choice, "custom")) {
        strsplit(input$morph_custom_palette %||% "", ",")[[1]]
      } else NULL
      resolve_palette(input$morph_palette_choice, length(input$selected_levels), custom_cols)
    })

    intensity_current_palette <- reactive({
      custom_cols <- if (identical(input$intensity_palette_choice, "custom")) {
        strsplit(input$intensity_custom_palette %||% "", ",")[[1]]
      } else NULL
      resolve_palette(input$intensity_palette_choice, length(input$selected_levels), custom_cols)
    })

    # the candidate set the "brackets to show" checkboxes are built from
    # the pairwise comparisons from the current metric groups
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
          checkboxInput(ns(paste0("bracket_show_", i)), label = labels[i], value = TRUE)
        })
      )
    })


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

    # A fixed graph height can't fit an arbitrary number of stacked
    # significance brackets without them overlapping this code allows the panel
    # to grow with the stack.
    plot_height_px <- reactive({
      fd <- filtered_df()
      if (nrow(fd) == 0) return(400)
      bracket_plot_height(fd, input$metric, "group", bracket_size = input$bracket_size,
                           bracket_spacing = input$bracket_spacing, title_size = input$title_size,
                           log_scale = input$log_scale,
                           selected_pairs = bracket_pairs_selected()) * 96
    })

    output$morph_plot_ui <- renderUI({
      plotOutput(ns("morph_plot"), height = paste0(plot_height_px(), "px"))
    })

    output$morph_plot <- renderPlot({
      fd <- filtered_df()
      validate(need(nrow(fd) > 0, "No data for the current sample selection."))
      plot_single_metric_pub(fd, input$metric, group_col = "group",
                              group_label = pretty_var(input$comparison_variable),
                              title = plot_title(), bracket_size = input$bracket_size,
                              bracket_spacing = input$bracket_spacing,
                              title_size = input$title_size, label_size = input$label_size,
                              title_justify = input$title_justify, log_scale = input$log_scale,
                              selected_pairs = bracket_pairs_selected(), palette = morph_current_palette())
    }, res = 96)

    output$counts_table <- renderDT({
      fd <- filtered_df()
      fd %>% count(group, name = "n_cells") %>%
        left_join(fd %>% distinct(group, image_id) %>% count(group, name = "n_images"), by = "group") %>%
        rename(!!pretty_var(input$comparison_variable) := group) %>%
        datatable(rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    morph_pairwise_current <- reactive({
      fd <- filtered_df()
      validate(need(n_distinct(fd$group) >= 2, "Select at least 2 samples to run statistics."))
      res <- compare_pairs_metric(fd, input$metric, "group")
      validate(need(!is.null(res) && nrow(res) > 0, "Not enough data to compare the selected samples."))
      res
    })

    output$pairwise_table <- renderDT({
      morph_pairwise_current() %>%
        mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
        datatable(rownames = FALSE, options = list(pageLength = 10))
    })

    output$download_pairwise <- downloadHandler(
      filename = function() paste0("pairwise_", input$metric, "_", Sys.Date(), ".csv"),
      content = function(file) {
        readr::write_csv(morph_pairwise_current(), file)
      }
    )

    output$download_pairwise_all <- downloadHandler(
      filename = function() paste0("pairwise_all_metrics_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- compare_pairs_all_metrics(filtered_df(), present_metrics_r(), "group")
        readr::write_csv(res, file)
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
                                     title_justify = input$title_justify, log_scale = input$log_scale,
                                     selected_pairs = bracket_pairs_selected(), palette = morph_current_palette())
        h <- bracket_plot_height(fd, input$metric, "group", bracket_size = input$bracket_size,
                                  bracket_spacing = input$bracket_spacing, title_size = input$title_size,
                                  log_scale = input$log_scale,
                                  selected_pairs = bracket_pairs_selected())
        ggsave(file, p, width = max(7, 1.2 + 1.1 * n_distinct(fd$group)), height = h, dpi = 300)
      }
    )

    # Built reactively (not static UI) so the tab's content follows whatever
    # has (or hasn't) been uploaded so far.
    output$intensity_tab_ui <- renderUI({
      idf <- intensity_store()
      if (nrow(idf) == 0) {
        return(div(style = "color: #a00; font-weight: bold;",
                    "Upload a channel_intensity_results_*.csv file in the sidebar to use this tab."))
      }
      tagList(
        fluidRow(
          column(4, selectInput(ns("intensity_channel"), "Channel",
                                 choices = sort(unique(idf$channel)))),
          column(4, selectInput(ns("intensity_statistic"), "Statistic",
                                 choices = intensity_stat_choices)),
          column(4, textInput(ns("intensity_plot_title"), "Plot title", value = "",
                               placeholder = "leave blank to use the metric name"))
        ),
        if ("seeded_status" %in% names(idf)) {
          fluidRow(
            column(4, selectInput(ns("intensity_seeded_filter"), "Seeded status",
                                   choices = c("All cells" = "all",
                                               "Seeded only" = "seeded",
                                               "Unseeded only" = "unseeded"),
                                   selected = "all"))
          )
        },
        fluidRow(
          column(4, selectInput(ns("intensity_palette_choice"), "Color palette",
                                 choices = palette_choices, selected = "okabe")),
          column(4, conditionalPanel(
            condition = sprintf("input['%s'] == 'custom'", ns("intensity_palette_choice")),
            textInput(ns("intensity_custom_palette"), "Custom colors",
                      placeholder = "Comma separated hex codes, e.g. #FF0000, #00A651")
          ))
        ),
        uiOutput(ns("intensity_bracket_pairs_ui")),
        br(),
        uiOutput(ns("intensity_plot_ui")),
        br(),
        downloadButton(ns("download_intensity_plot"), "Download plot (PNG)"),
        br(), br(),
        fluidRow(
          column(6, downloadButton(ns("download_intensity_pairwise"), "Download this statistic's table (CSV)")),
          column(6, downloadButton(ns("download_intensity_pairwise_all"), "Download all statistics (CSV)"))
        ),
        br(),
        DTOutput(ns("intensity_pairwise_table")),
        br(),
        DTOutput(ns("intensity_counts_table"))
      )
    })

    intensity_pairwise_current <- reactive({
      fd <- intensity_filtered_df()
      validate(need(n_distinct(fd$group) >= 2, "Select at least 2 samples to run statistics."))
      res <- compare_pairs_metric(fd, input$intensity_statistic, "group")
      validate(need(!is.null(res) && nrow(res) > 0, "Not enough data to compare the selected samples."))
      res
    })

    output$intensity_pairwise_table <- renderDT({
      intensity_pairwise_current() %>%
        mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
        datatable(rownames = FALSE, options = list(pageLength = 10))
    })

    filename_seeded_suffix <- reactive(switch(input$intensity_seeded_filter %||% "all",
                                               seeded = "_seeded", unseeded = "_unseeded", ""))

    output$download_intensity_pairwise <- downloadHandler(
      filename = function() paste0("intensity_pairwise_", input$intensity_statistic,
                                    filename_seeded_suffix(), "_", Sys.Date(), ".csv"),
      content = function(file) {
        readr::write_csv(intensity_pairwise_current(), file)
      }
    )

    output$download_intensity_pairwise_all <- downloadHandler(
      filename = function() paste0("intensity_pairwise_all_statistics",
                                    filename_seeded_suffix(), "_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- compare_pairs_all_metrics(intensity_filtered_df(), unname(intensity_stat_choices), "group")
        readr::write_csv(res, file)
      }
    )

    seeded_filter_suffix <- reactive(switch(input$intensity_seeded_filter %||% "all",
                                             seeded = " (Seeded)", unseeded = " (Unseeded)", ""))

    intensity_plot_title <- reactive({
      base_title <- if (nzchar(trimws(input$intensity_plot_title))) input$intensity_plot_title
                    else metric_label(input$intensity_statistic)
      paste0(base_title, seeded_filter_suffix())
    })

    intensity_sig_pairs <- reactive({
      fd <- intensity_filtered_df()
      req(input$intensity_statistic)
      if (n_distinct(fd$group) < 2 || n_distinct(fd$group) > 8) return(NULL)
      pw <- compare_pairs_metric(fd, input$intensity_statistic, "group")
      if (is.null(pw)) return(NULL)
      sig <- pw %>% filter(p.adj < 0.05)
      if (nrow(sig) == 0) return(NULL)
      sig
    })

    output$intensity_bracket_pairs_ui <- renderUI({
      sig <- intensity_sig_pairs()
      if (is.null(sig)) return(NULL)
      labels <- pair_labels(sig)
      tagList(
        tags$strong("Significance brackets to show"),
        lapply(seq_along(labels), function(i) {
          tags$label(class = "checkbox-inline",
                     tags$input(id = ns(paste0("intensity_bracket_show_", i)), type = "checkbox",
                                class = "shiny-input-checkbox", checked = "checked"),
                     labels[i])
        })
      )
    })

    intensity_bracket_pairs_selected <- reactive({
      sig <- intensity_sig_pairs()
      if (is.null(sig)) return(NULL)
      labels <- pair_labels(sig)
      keep <- vapply(seq_along(labels), function(i) {
        v <- input[[paste0("intensity_bracket_show_", i)]]
        if (is.null(v)) TRUE else v
      }, logical(1))
      labels[keep]
    })

    intensity_plot_height_px <- reactive({
      fd <- intensity_filtered_df()
      # No req() here: an empty fd must still let renderPlot's own validate()
      # message through, rather than silently killing the height callback and
      # leaving blank space with no explanation.
      if (nrow(fd) == 0) return(400)
      bracket_plot_height(fd, input$intensity_statistic, "group", bracket_size = input$bracket_size,
                           bracket_spacing = input$bracket_spacing, title_size = input$title_size,
                           log_scale = input$log_scale,
                           selected_pairs = intensity_bracket_pairs_selected()) * 96
    })

    # See morph_plot_ui above for why the container is sized explicitly here
    # rather than relying solely on renderPlot's height callback.
    output$intensity_plot_ui <- renderUI({
      plotOutput(ns("intensity_plot"), height = paste0(intensity_plot_height_px(), "px"))
    })

    output$intensity_plot <- renderPlot({
      fd <- intensity_filtered_df()
      validate(need(nrow(fd) > 0, "No data for the current sample selection."))
      plot_single_metric_pub(fd, input$intensity_statistic, group_col = "group",
                              group_label = pretty_var(input$comparison_variable),
                              title = intensity_plot_title(), bracket_size = input$bracket_size,
                              bracket_spacing = input$bracket_spacing,
                              title_size = input$title_size, label_size = input$label_size,
                              title_justify = input$title_justify, log_scale = input$log_scale,
                              selected_pairs = intensity_bracket_pairs_selected(), palette = intensity_current_palette())
    }, res = 96)

    output$intensity_counts_table <- renderDT({
      fd <- intensity_filtered_df()
      fd %>% count(group, name = "n_cells") %>%
        left_join(fd %>% distinct(group, image_id) %>% count(group, name = "n_images"), by = "group") %>%
        rename(!!pretty_var(input$comparison_variable) := group) %>%
        datatable(rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    output$download_intensity_plot <- downloadHandler(
      filename = function() paste0("intensity_", input$intensity_channel, "_",
                                    input$intensity_statistic, filename_seeded_suffix(),
                                    "_", Sys.Date(), ".png"),
      content = function(file) {
        fd <- intensity_filtered_df()
        p <- plot_single_metric_pub(fd, input$intensity_statistic, group_col = "group",
                                     group_label = pretty_var(input$comparison_variable),
                                     title = intensity_plot_title(), bracket_size = input$bracket_size,
                                     bracket_spacing = input$bracket_spacing,
                                     title_size = input$title_size, label_size = input$label_size,
                                     title_justify = input$title_justify, log_scale = input$log_scale,
                                     selected_pairs = intensity_bracket_pairs_selected(), palette = intensity_current_palette())
        h <- bracket_plot_height(fd, input$intensity_statistic, "group", bracket_size = input$bracket_size,
                                  bracket_spacing = input$bracket_spacing, title_size = input$title_size,
                                  log_scale = input$log_scale,
                                  selected_pairs = intensity_bracket_pairs_selected())
        ggsave(file, p, width = max(7, 1.2 + 1.1 * n_distinct(fd$group)), height = h, dpi = 300)
      }
    )
  })
}
