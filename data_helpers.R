# Shared data loading / parsing / plotting / stats helpers for
# analysis_v2.1.Rmd and app.R. Sourced by both so behavior never drifts
# between the static report and the interactive app.
#
# Assumes the caller has already loaded tidyverse (dplyr/stringr/purrr/readr)
# and rstatix.

okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00",
           "#CC79A7", "#56B4E9", "#F0E442", "#000000")

okabe_pal <- grDevices::colorRampPalette(okabe)

# Palette choices for the app's "Color palette" dropdown. All but "custom"
# are generated via base R's hcl.colors() (no extra package dependency).
palette_choices <- c(
  "Okabe-Ito (colorblind-safe)" = "okabe",
  "Viridis"                     = "viridis",
  "Set 2"                       = "set2",
  "Dark 2"                      = "dark2",
  "Pastel 1"                    = "pastel1",
  "Custom"                      = "custom"
)

# Resolves a palette dropdown choice (plus optional user-entered custom hex
# codes) into a vector of n fill colors for plot_single_metric_pub(). Falls
# back to the Okabe-Ito default if "custom" is selected but no valid hex
# codes were entered.
resolve_palette <- function(name, n, custom_colors = NULL) {
  n <- max(as.integer(n), 1L)
  if (identical(name, "custom")) {
    cols <- trimws(custom_colors)
    cols <- cols[nzchar(cols)]
    cols <- ifelse(grepl("^#", cols), cols, paste0("#", cols))
    cols <- cols[grepl("^#[0-9A-Fa-f]{6}$", cols)]
    if (length(cols) >= 1) {
      if (length(cols) == 1) cols <- rep(cols, 2)
      return(grDevices::colorRampPalette(cols)(n))
    }
    name <- "okabe"
  }
  switch(name,
    okabe   = okabe_pal(n),
    viridis = grDevices::hcl.colors(n, "Viridis"),
    set2    = grDevices::hcl.colors(n, "Set 2"),
    dark2   = grDevices::hcl.colors(n, "Dark 2"),
    pastel1 = grDevices::hcl.colors(n, "Pastel 1"),
    okabe_pal(n)
  )
}

theme_pub <- theme_minimal(base_size = 15) +
  theme(
    axis.line          = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks         = element_line(linewidth = 0.6, colour = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.3, colour = "grey88"),
    legend.position    = "right",
    legend.title       = element_text(face = "bold"),
    strip.text         = element_text(face = "bold"),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    plot.title         = element_text(face = "bold", size = rel(1.05)),
    plot.subtitle      = element_text(colour = "grey30", size = rel(0.85))
  )
theme_set(theme_pub)

scale_colour_discrete <- function(...) discrete_scale("colour", palette = okabe_pal, ...)
scale_fill_discrete   <- function(...) discrete_scale("fill", palette = okabe_pal, ...)

# ---- Metric metadata ---------------------------------------------------
# Columns the Fiji macro produces, organized by what they describe
metric_categories <- list(
  size           = c("area_um2", "convexhull_area_um2",
                      "perimeter_um", "convexhull_perimeter_um"),
  ramification   = c("ramification_index", "ramification_index_2d",
                      "circularity", "solidity",
                      "n_branches", "n_junctions", "tree_length_um",
                      "avg_branch_length_um", "max_branch_length_um"),
  arms           = c("n_tips", "n_triple_points", "n_quadruple_points"),
  shape_polarity = c("aspect_ratio", "polarity_offset_um")
)
all_cell_metrics <- unlist(metric_categories, use.names = FALSE)

pretty_var <- function(v) {
  labels <- c(project = "Project", well_position = "Well Position",
              site = "Site/FOV", wavelength = "Wavelength", image = "Image")
  unname(labels[v])
}


metric_label <- function(m) {
  um  <- paste0(intToUtf8(0xB5), "m")
  um2 <- paste0(intToUtf8(0xB5), "m", intToUtf8(0xB2))
  labels <- c(
    area_um2                = paste0("Cell Area (", um2, ")"),
    convexhull_area_um2     = paste0("Convex Hull Area (", um2, ")"),
    perimeter_um            = paste0("Perimeter (", um, ")"),
    convexhull_perimeter_um = paste0("Convex Hull Perimeter (", um, ")"),
    ramification_index      = "Ramification Index",
    ramification_index_2d   = "Ramification Index (2D)",
    circularity             = "Circularity",
    solidity                = "Solidity",
    n_branches              = "Number of Branches",
    n_junctions             = "Number of Junctions",
    tree_length_um          = paste0("Tree Length (", um, ")"),
    avg_branch_length_um    = paste0("Avg. Branch Length (", um, ")"),
    max_branch_length_um    = paste0("Max. Branch Length (", um, ")"),
    n_tips                  = "Number of Tips (Arms)",
    n_triple_points         = "Triple Points",
    n_quadruple_points      = "Quadruple Points",
    aspect_ratio            = "Aspect Ratio",
    polarity_offset_um      = paste0("Polarity Offset (", um, ")"),
    mean_intensity          = "Mean Intensity (a.u.)",
    median_intensity        = "Median Intensity (a.u.)",
    min_intensity           = "Min Intensity (a.u.)",
    max_intensity           = "Max Intensity (a.u.)",
    stddev_intensity        = "Std. Dev. Intensity (a.u.)",
    integrated_density      = "Integrated Density (a.u.)",
    raw_integrated_density  = "Raw Integrated Density (a.u.)",
    area_px                 = "Cell Area (px)"
  )
  if (m %in% names(labels)) unname(labels[m]) else str_to_title(gsub("_", " ", m))
}

parse_cell_name <- function(cell_name) {
  sample_id   <- str_extract(cell_name, "^[^_]+")
  cell_number <- as.integer(str_extract(cell_name, "(?<=_cell)\\d+(?=\\.tif$)"))

  image_id <- str_remove(cell_name, "_cell\\d+\\.tif$")
  image_id <- if_else(str_detect(image_id, "_X\\d+_Y\\d+_Z\\d+\\.tif$"),
                       str_remove(image_id, "_X\\d+_Y\\d+_Z\\d+\\.tif$"),
                       image_id)
  tibble(sample_id = sample_id, image_id = image_id, cell_number = cell_number)
}

# Automatic filename metadata parsing
# For "morphology_results_<project>_<well_position>_<site>_<wavelength>_<image>.csv"
# pulls out project/well_position/site/wavelength/image regardless of how
# many underscore-separated tokens the project name itself contains, by
# recognizing the fixed-format tokens (B07/C07/... = well position - any
# plate row letter(s) followed by digits, s7 = site/FOV, w2 = wavelength,
# trailing digits = image) and treating everything left over as the
# project name. There is no separate "time" field in this naming
# convention - any other token (e.g. "12h") stays folded into project.

parse_filename_metadata <- function(filename) {
  base <- basename(filename)
  base <- str_remove(base, "\\.csv$")
  base <- str_remove(base, "^morphology_results_")

  image <- str_extract(base, "(?<=_)\\d+$")
  rest  <- str_remove(base, "_\\d+$")
  toks  <- str_split(rest, "_")[[1]]

  take <- function(toks, pattern) {
    hit <- str_which(toks, pattern)
    if (length(hit) == 0) return(list(value = NA_character_, toks = toks))
    list(value = toks[hit[1]], toks = toks[-hit[1]])
  }


  r <- take(toks, regex("^[A-Z]{1,2}\\d+$")); well_position <- r$value; toks <- r$toks
  r <- take(toks, regex("^w\\d+$", ignore_case = TRUE)); wavelength    <- r$value; toks <- r$toks
  r <- take(toks, regex("^s\\d+$", ignore_case = TRUE)); site          <- r$value; toks <- r$toks

  project <- if (length(toks) > 0) paste(toks, collapse = "_") else NA_character_

  tibble(project = project, well_position = well_position,
         site = site, wavelength = wavelength, image = image)
}

# Data loading, with a merged-CSV cache

load_morphology_data <- function(dir = ".",
                                  metadata_file = file.path(dir, "sample_metadata.csv"),
                                  cache_file = file.path(dir, "merged_morphology_data.csv"),
                                  force_rebuild = FALSE) {
  csv_files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!basename(csv_files) %in%
                           basename(c(metadata_file, cache_file, "merged_intensity_data.csv"))]
  csv_files <- csv_files[!grepl("^channel_intensity_results_.*\\.csv$", basename(csv_files))]

  cache_exists <- file.exists(cache_file)
  cache_fresh <- cache_exists && (
    length(csv_files) == 0 ||
      file.info(cache_file)$mtime >= max(file.info(csv_files)$mtime)
  )

  if (!force_rebuild && cache_fresh) {
    message("Loading merged data cache '", basename(cache_file), "' (",
            "delete it, or set force_rebuild = TRUE, to re-parse the source CSVs).")
    df <- readr::read_csv(cache_file, show_col_types = FALSE)
  } else {
    if (length(csv_files) == 0)
      stop("No CSV files found in '", normalizePath(dir), "'. ",
           "Point `dir` at the macro's analysis_csv/ output folder.")

    message("Parsing ", length(csv_files), " source CSV(s) - this is slow the ",
            "first time (and whenever new CSVs are added); after this it will ",
            "read the cached '", basename(cache_file), "' instead.")

    read_one <- function(f) {
      d <- readr::read_csv(f, show_col_types = FALSE)
      if (!all(c("cell_name", "area_um2") %in% names(d))) {
        message("Skipping '", basename(f),
                "' - not a morphology CSV (missing cell_name/area_um2).")
        return(NULL)
      }
      dplyr::mutate(d, source_file = basename(f))
    }
    raw <- purrr::map_dfr(csv_files, read_one)
    if (nrow(raw) == 0)
      stop("Found CSVs in '", normalizePath(dir), "' but none had the expected ",
           "columns (need at least cell_name and area_um2). Make sure you are ",
           "pointing at the macro's analysis_csv/ output.")

    file_meta <- purrr::map_dfr(unique(raw$source_file), function(f) {
      bind_cols(tibble(source_file = f), parse_filename_metadata(f))
    })

    df <- bind_cols(parse_cell_name(raw$cell_name), raw) %>%
      left_join(file_meta, by = "source_file") %>%
      mutate(ramification_index_2d = perimeter_um / (2 * sqrt(pi * area_um2)))

    
    if (!all(c("centroid_x_um", "centroid_y_um") %in% names(df))) {
      df <- df %>% mutate(
        centroid_x_um = as.numeric(str_extract(cell_name, "(?<=_X)\\d+")),
        centroid_y_um = as.numeric(str_extract(cell_name, "(?<=_Y)\\d+"))
      )
    }

    tryCatch({
      readr::write_csv(df, cache_file)
      message("Wrote merged data cache: ", normalizePath(cache_file))
    }, error = function(e) {
      warning("Could not write merged data cache '", cache_file, "': ", conditionMessage(e))
    })
  }

  if (file.exists(metadata_file)) {
    meta <- readr::read_csv(metadata_file, show_col_types = FALSE)
    if ("sample_id" %in% names(meta)) {
      meta <- meta %>% select(-any_of(c("project", "well_position", "site", "wavelength", "image")))
      df <- left_join(df, meta, by = "sample_id")
    } else {
      message("sample_metadata.csv found but has no 'sample_id' column - ignoring it.")
    }
  }

  attr(df, "n_csv_loaded") <- n_distinct(df$source_file)
  df
}

# ---- Channel intensity data (from the Channel Intensity Quantification --
# macro's channel_intensity_results_*.csv) ---------------------------------
# Unlike morphology CSVs (one per source image), the quantification macro
# always produces a single combined file per run, long-format: one row per
# cell per channel. There's nothing to merge, so - unlike morphology data -
# this reads the native file directly with no on-disk cache. Returns an
# empty tibble (with a message, not an error) if the macro hasn't been run
# yet, so the app still starts fine; errors if more than one matching file
# is present, since that means a stale run wasn't cleaned up.
load_intensity_data <- function(dir = ".") {
  csv_files <- list.files(dir, pattern = "^channel_intensity_results_.*\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    message("No 'channel_intensity_results_*.csv' files found in '", normalizePath(dir),
            "' - the Channel Intensity tab will be empty until you run the ",
            "Channel Intensity Quantification macro.")
    return(tibble())
  }

  if (length(csv_files) > 1) {
    stop("Found multiple 'channel_intensity_results_*.csv' files in '", normalizePath(dir), "':\n  ",
         paste(basename(csv_files), collapse = "\n  "),
         "\nThe quantification macro only produces one at a time - delete the stale run(s) ",
         "and keep just the current file.")
  }

  message("Reading channel intensity data: ", basename(csv_files))
  d <- readr::read_csv(csv_files, show_col_types = FALSE)
  needed <- c("cell_name", "source_image", "channel", "mean_intensity")
  if (!all(needed %in% names(d))) {
    stop("'", basename(csv_files), "' is missing expected columns (need at least ",
         "cell_name/source_image/channel/mean_intensity). Make sure you're pointing at ",
         "the Channel Intensity Quantification macro's output.")
  }
  d
}

# Attaches the same project/well_position/site/wavelength/image (+
# sample_id/image_id) columns morphology_df already carries, matched by
# cell_name - the two pipelines' cell_name values are identical for the
# same cell, so this is an exact 1:1 join, no re-parsing of filenames needed.
attach_cell_metadata <- function(intensity_df, morphology_df) {
  if (nrow(intensity_df) == 0) return(intensity_df)
  meta <- morphology_df %>%
    distinct(cell_name, sample_id, image_id, cell_number,
             project, well_position, site, wavelength, image)
  intensity_df %>% left_join(meta, by = "cell_name")
}

# Shared analysis / plotting helpers
aggregate_to_image <- function(data, metrics = all_cell_metrics, group_col = "group") {
  metrics <- intersect(metrics, names(data))
  data %>%
    group_by(.data[[group_col]], image_id) %>%
    summarise(
      across(all_of(metrics), \(x) median(x, na.rm = TRUE)),
      n_cells = dplyr::n(),
      project = dplyr::first(project), well_position = dplyr::first(well_position),
      site = dplyr::first(site), wavelength = dplyr::first(wavelength),
      .groups = "drop"
    )
}

compare_metric <- function(data, metric, group_col = "group") {
  data <- data %>% filter(!is.na(.data[[metric]]))
  n_groups <- n_distinct(data[[group_col]])
  if (n_groups < 2) return(NULL)

  group_sizes <- data %>% count(.data[[group_col]])
  group_variances <- data %>% group_by(.data[[group_col]]) %>%
    summarise(v = var(.data[[metric]], na.rm = TRUE), .groups = "drop")
  normal <- FALSE
  if (all(group_sizes$n >= 3) && all(group_variances$v > 0, na.rm = TRUE)) {
    sh <- data %>% group_by(.data[[group_col]]) %>% rstatix::shapiro_test(!!sym(metric))
    normal <- all(sh$p > 0.05)
  }
  f <- as.formula(paste(metric, "~", group_col))

  if (n_groups == 2) {
    res <- if (normal) rstatix::t_test(data, f) %>% mutate(test = "Welch t-test")
           else rstatix::wilcox_test(data, f) %>% mutate(test = "Mann-Whitney U")
  } else {
    res <- if (normal) rstatix::anova_test(data, f) %>% as_tibble() %>% mutate(test = "One-way ANOVA")
           else rstatix::kruskal_test(data, f) %>% mutate(test = "Kruskal-Wallis")
  }
  p_col <- intersect(c("p", "p.adj"), names(res))[1]
  tibble(metric = metric, test = res$test[1], n_groups = n_groups, normal = normal, p = res[[p_col]][1])
}

compare_all_metrics <- function(data, metrics = all_cell_metrics, group_col = "group") {
  metrics <- intersect(metrics, names(data))
  if (length(metrics) == 0) return(tibble())
  purrr::map_dfr(metrics, ~ compare_metric(data, .x, group_col)) %>%
    mutate(p_fdr = if (dplyr::n() > 0) p.adjust(p, method = "BH") else numeric(0)) %>%
    arrange(p_fdr)
}


compare_pairs_metric <- function(data, metric, group_col = "group") {
  data <- data %>% filter(!is.na(.data[[metric]]))
  n_groups <- n_distinct(data[[group_col]])
  if (n_groups < 2) return(NULL)

  group_sizes <- data %>% count(.data[[group_col]])
  group_variances <- data %>% group_by(.data[[group_col]]) %>%
    summarise(v = var(.data[[metric]], na.rm = TRUE), .groups = "drop")
  normal <- FALSE
  if (all(group_sizes$n >= 3) && all(group_variances$v > 0, na.rm = TRUE)) {
    sh <- data %>% group_by(.data[[group_col]]) %>% rstatix::shapiro_test(!!sym(metric))
    normal <- all(sh$p > 0.05)
  }
  f <- as.formula(paste(metric, "~", group_col))

  res <- if (normal) rstatix::t_test(data, f, p.adjust.method = "BH") %>% mutate(test = "Welch t-test")
         else rstatix::wilcox_test(data, f, p.adjust.method = "BH") %>% mutate(test = "Mann-Whitney U")


  if (!"p.adj" %in% names(res)) res <- res %>% mutate(p.adj = p)
  if (!"p.adj.signif" %in% names(res))
    res <- rstatix::add_significance(res, p.col = "p.adj", output.col = "p.adj.signif")

  res %>% mutate(metric = metric, .before = 1) %>% select(-any_of(".y."))
}

compare_pairs_all_metrics <- function(data, metrics = all_cell_metrics, group_col = "group") {
  metrics <- intersect(metrics, names(data))
  if (length(metrics) == 0) return(tibble())
  purrr::map_dfr(metrics, ~ compare_pairs_metric(data, .x, group_col))
}


pair_labels <- function(pw) paste(pw$group1, "vs", pw$group2)

# The significant pairs to bracket, with rstatix's computed y.position for
# each (in the metric's own data units) - the single source of truth both
# plot_single_metric_pub() and bracket_plot_height() build on, so the
# reserved margin always matches where rstatix actually places the brackets
# (which is based on rstatix's own internal range, not necessarily the full
# min/max of the plotted data - see bracket_margin_inches()).
bracket_layout <- function(data, metric, group_col = "group", bracket_spacing = 0.15,
                            max_groups_for_brackets = 8, selected_pairs = NULL) {
  n_grp <- n_distinct(data[[group_col]])
  if (n_grp < 2 || n_grp > max_groups_for_brackets) return(NULL)
  pw <- compare_pairs_metric(data, metric, group_col)
  if (is.null(pw)) return(NULL)
  sig <- pw %>% filter(p.adj < 0.05)
  if (!is.null(selected_pairs)) sig <- sig %>% filter(pair_labels(sig) %in% selected_pairs)
  if (nrow(sig) == 0) return(NULL)
  rstatix::add_xy_position(sig, x = group_col, step.increase = bracket_spacing)
}

# Inches of top margin needed so the highest bracket in sig_pos isn't
# clipped, expressed as "how many panel-heights above the visible axis range
# does it sit" - measured directly from rstatix's own y.position values
# rather than estimated, since rstatix scales step.increase using its own
# internal range (max per group vs min-of-per-group-max) which doesn't match
# the plotted data's full min/max closely enough to predict analytically.
bracket_margin_inches <- function(sig_pos, values, bracket_size, base_height = 5.5) {
  values <- values[!is.na(values)]
  y_range <- range(values)
  y_pad <- diff(y_range) * 0.05
  if (y_pad == 0) y_pad <- 1
  panel_span <- diff(y_range) + 2 * y_pad
  overshoot <- max(0, max(sig_pos$y.position, na.rm = TRUE) - (y_range[2] + y_pad))
  (overshoot / panel_span) * base_height + nrow(sig_pos) * (0.12 + bracket_size * 0.03)
}

# Inches reserved for the plot title, now drawn as its own annotation above
# the top bracket (see plot_single_metric_pub()) rather than ggplot's usual
# title row, which sits immediately above the panel and so got run over by
# any bracket tall enough to overflow into it.
title_margin_inches <- function(title_size) 0.15 + title_size * 0.015

# A few outlier cells (common in intensity channels - saturated pixels etc.)
# can be orders of magnitude above the bulk of the data; since the axis is
# locked to the real range (see plot_single_metric_pub()), that squashes
# every violin into a flat sliver. log10-transforming the plotted values (and
# the bracket y.position, which rstatix computed on the untransformed scale)
# keeps them in sync without changing what the underlying stats test runs on.
log10_sig_pos <- function(sig_pos) {
  sig_pos$y.position <- log10(pmax(sig_pos$y.position, 1e-9))
  sig_pos
}

plot_single_metric_pub <- function(data, metric, group_col = "group", group_label = NULL,
                                    title = NULL, max_groups_for_brackets = 8,
                                    bracket_size = 4, bracket_spacing = 0.15,
                                    title_size = 16, label_size = 15,
                                    title_justify = c("left", "center"),
                                    log_scale = FALSE,
                                    selected_pairs = NULL,
                                    palette = NULL) {
  title_justify <- match.arg(title_justify)
  n_grp <- n_distinct(data[[group_col]])

  plt_df <- data %>%
    select(all_of(group_col), all_of(metric)) %>%
    rename(value = all_of(metric)) %>%
    filter(!is.na(value))
  if (log_scale) plt_df <- plt_df %>% filter(value > 0) %>% mutate(value = log10(value))

  # Fix the y-axis to the real data range (with a small fixed pad) instead of
  # letting it auto-expand to fit the brackets - that expansion is what used
  # to squash the violins down into a fraction of the panel. Brackets are
  # drawn above this range in reserved margin space instead (clip = "off").
  y_range <- range(plt_df$value, na.rm = TRUE)
  y_pad <- diff(y_range) * 0.05
  if (y_pad == 0) y_pad <- 1
  panel_span <- diff(y_range) + 2 * y_pad

  p <- ggplot(plt_df, aes(x = .data[[group_col]], y = value, fill = .data[[group_col]])) +
    geom_violin(alpha = 0.25, colour = "black", linewidth = 0.4, trim = TRUE, na.rm = TRUE) +
    geom_boxplot(width = 0.15, alpha = 0.9, colour = "black", outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.4, shape = 16, colour = "grey20") +
    labs(x = group_label, y = metric_label(metric), fill = group_label) +
    theme(legend.position = if (n_grp > 8) "none" else "right",
          axis.title  = element_text(size = label_size),
          axis.text.x = element_text(angle = 30, hjust = 1, size = label_size * 0.8),
          axis.text.y = element_text(size = label_size * 0.8))

  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }

  if (log_scale) {
    p <- p + scale_y_continuous(labels = function(x) scales::comma(round(10^x)))
  }

  # The title used to be ggplot's own title row, which sits immediately
  # above the panel - so any bracket tall enough to overflow past the panel
  # (via clip = "off" below) ran straight into it. Drawing it as a manual
  # annotation, positioned above the top bracket's own y.position (or just
  # above the axis range if there are no brackets), keeps it clear no matter
  # how many brackets are stacked.
  bracket_margin_in <- 0
  top_y <- y_range[2] + y_pad
  sig_pos <- bracket_layout(data, metric, group_col, bracket_spacing, max_groups_for_brackets, selected_pairs)
  if (!is.null(sig_pos)) {
    if (log_scale) sig_pos <- log10_sig_pos(sig_pos)
    bracket_margin_in <- bracket_margin_inches(sig_pos, plt_df$value, bracket_size)
    top_y <- max(sig_pos$y.position, na.rm = TRUE)
    p <- p + ggpubr::stat_pvalue_manual(sig_pos, label = "p.adj.signif",
                                         tip.length = 0.01, size = bracket_size,
                                         bracket.size = bracket_size / 10)
  }

  title_extra_in <- title_margin_inches(title_size)
  total_margin_in <- bracket_margin_in + title_extra_in
  data_per_inch <- panel_span / 5.5
  # top_y already sits bracket_margin_in's "overshoot" component above the
  # axis top (y_range[2] + y_pad) - anchoring from top_y and then adding the
  # full bracket_margin_in again double-counted that overshoot, pushing the
  # title (and the blank space below it) further up the more brackets were
  # stacked. Anchoring from the axis top instead reserves exactly
  # total_margin_in above the axis, matching plot.margin below.
  title_y <- (y_range[2] + y_pad) + (bracket_margin_in + title_extra_in * 0.5) * data_per_inch

  title_x   <- if (title_justify == "left") 0.5 else (n_grp + 1) / 2
  title_hjust <- if (title_justify == "left") 0 else 0.5

  p <- p +
    annotate("text", x = title_x, y = title_y,
             label = if (!is.null(title)) title else metric,
             hjust = title_hjust, vjust = 0.5, fontface = "bold",
             size = title_size / .pt)

  p + coord_cartesian(ylim = c(y_range[1] - y_pad, y_range[2] + y_pad), clip = "off") +
    theme(plot.margin = margin(t = total_margin_in, r = 0.08, b = 0.08, l = 0.08,
                                unit = "in"))
}


n_significant_pairs <- function(data, metric, group_col = "group", max_groups_for_brackets = 8,
                                 selected_pairs = NULL) {
  sig_pos <- bracket_layout(data, metric, group_col, max_groups_for_brackets = max_groups_for_brackets,
                             selected_pairs = selected_pairs)
  if (is.null(sig_pos)) 0L else nrow(sig_pos)
}


bracket_plot_height <- function(data, metric, group_col = "group", bracket_size = 4,
                                 bracket_spacing = 0.15, title_size = 16,
                                 max_groups_for_brackets = 8, log_scale = FALSE,
                                 base_height = 5.5, selected_pairs = NULL) {
  sig_pos <- bracket_layout(data, metric, group_col, bracket_spacing, max_groups_for_brackets, selected_pairs)
  values <- data[[metric]]
  if (log_scale) values <- log10(values[!is.na(values) & values > 0])
  bracket_in <- if (is.null(sig_pos)) 0 else {
    if (log_scale) sig_pos <- log10_sig_pos(sig_pos)
    bracket_margin_inches(sig_pos, values, bracket_size, base_height)
  }
  base_height + bracket_in + title_margin_inches(title_size)
}


save_morphology_plots_to_disk <- function(data, out_dir = "morphology_plots", group_col = "group",
                                           group_label = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  plot_width <- max(7, 1.2 + 1.1 * n_distinct(data[[group_col]]))

  saved <- character(0)
  for (cat_name in names(metric_categories)) {
    metrics <- intersect(metric_categories[[cat_name]], names(data))
    if (length(metrics) == 0) next
    cat_dir <- file.path(out_dir, cat_name)
    dir.create(cat_dir, showWarnings = FALSE, recursive = TRUE)

    for (m in metrics) {
      p <- plot_single_metric_pub(data, m, group_col, group_label)
      fname <- file.path(cat_dir, paste0(m, ".png"))
      ggsave(fname, p, width = plot_width, height = 5.5, dpi = 300)
      saved <- c(saved, fname)
    }
  }
  invisible(saved)
}
