#' Get assay labels for volcano plots
#'
#' @description
#' Returns the list of assay labels to annotate in a volcano plot.
#' If the current contrast has a list defined in `config$label_on_plots`,
#' that list is used. Otherwise, the function selects the top `n` upregulated
#' and downregulated proteins (by adjusted p-value).
#'
#' @param res A data frame of results containing columns `estimate`,
#'   `Adjusted_pval`, `Assay_Panel`, and optionally `Threshold`.
#' @param contrast Character scalar. Name of the current contrast.
#' @param config Parsed YAML config list that contains `label_on_plots`.
#' @param top_n Integer. Number of top up- and down-regulated proteins
#'   to include when no config entry exists. Default is 5.
#'
#' @return A character vector of assay labels.
#'
#' @export
get_assay_labels <- function(res, contrast, config, top_n = 5) {
  if (contrast %in% names(config$label_on_plots)) {
    # use list from config
    assay_list <- config$label_on_plots[[contrast]]
  } else {
    # default: top 5 up and down
    sig_up <- res %>%
      dplyr::filter(Threshold == "Significant", estimate > 0) %>%
      dplyr::arrange(Adjusted_pval) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(Assay_Panel)

    sig_down <- res %>%
      dplyr::filter(Threshold == "Significant", estimate <= 0) %>%
      dplyr::arrange(Adjusted_pval) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(Assay_Panel)

    assay_list <- unique(c(sig_up, sig_down))
  }

  return(assay_list)
}

#' Build and save a volcano plot with metadata
#'
#' @description
#' Generates a volcano plot for one contrast and test type. Annotates with
#' sample sizes, saves full-range and optional zoomed PDFs, and returns an
#' interactive Plotly object with click-through to the Human Protein Atlas.
#'
#' @param res Data frame with at least `Assay`, `Panel`, `estimate` or `log2FC`,
#'   and `Adjusted_pval`.
#' @param se_list Named list of SummarizedExperiment objects.
#' @param contrast Character, name of the contrast.
#' @param test_type Character test type (e.g. `"anova"`, `"limma"`).
#' @param config Config list with`alpha`.
#' @param out_dir PDF output directory from config.
#' @param parameters Contrast parameter list from \code{get_params()}.
#' @param export_2nd Logical, save zoomed PDF (default TRUE).
#' @param export_2nd_y_lim Numeric, y-axis cap for zoomed plot (default 5).
#' @param seed Random seed for label placement (default 2).
#' @param vert_line Logical, draw vertical line at x=0 (default TRUE).
#'
#' @return Invisibly, list with Plotly object and PDF paths.
#' @export
make_volcano <- function(res,
                         se_list,
                         contrast,
                         test_type,
                         config,
                         out_dir,
                         parameters,
                         export_2nd       = TRUE,
                         export_2nd_y_lim = 5,
                         seed             = 2,
                         vert_line        = TRUE) {
  # Uniquely ID each protein
  res <- dplyr::mutate(res, Assay_Panel = paste(Assay, Panel, sep = "_"))
  alpha <- parameters$alpha

  # decide which column goes on the x-axis and how to label it
  if ("estimate" %in% names(res)) {
    x_var <- "estimate"
    x_lab <- "estimate (adj. log2FC)"
  } else if ("log2FC" %in% names(res)) {
    x_var <- "log2FC"
    x_lab <- "log2FC"
  } else {
    stop("Need either ‘estimate’ or ‘log2FC’ in the results table.")
  }

  # Get the SummarizedExperiment by name
  se_obj <- se_list[[ parameters$se_object_name ]]
  group_col <- parameters$column
  levels    <- c(parameters$treatment, parameters$control)

  # Count samples in each arm
  cd <- SummarizedExperiment::colData(se_obj)
  n1 <- sum(cd[[group_col]] == levels[1], na.rm = TRUE)
  n2 <- sum(cd[[group_col]] == levels[2], na.rm = TRUE)

  # Build a subtitle like "disease_vs_control (disease: 42 vs control: 20)"
  subtitle_txt <- sprintf(
    "%s (%s: n = %d vs %s: n = %d)",
    parameters$name, levels[1], n1, levels[2], n2
  )

  # Which Assay_Panels to label
  assay_list <- get_assay_labels(res, contrast, config)

  brms_flag <- grepl("brms", parameters$name)

  # helper to build ggplot
  build_volcano <- function(df, y_var, x_var, x_lab,
                          shape_var = NULL,
                          title_txt, subtitle_txt) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
            x     = !!rlang::sym(x_var),
            y     = !!rlang::sym(y_var),
            color = Threshold
          )) +
      ggplot2::geom_point() +
      ggplot2::labs(
        x        = x_lab,
        y        = ifelse(brms_flag, "-log10(adj. q-value", "-log10(adj. p-value)"),
        color    = paste("FDR <", alpha),
        title    = title_txt,
        subtitle = subtitle_txt
      ) +
      ggplot2::scale_color_manual(
        values = c(Significant = "red",
                   `Non-significant` = "black")
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(alpha),
        linetype   = "dashed",
        color     = "black"
      ) +
      ggplot2::annotate(
        "text",
        x     = Inf, y = -log10(alpha),
        label = paste0("-log10(", alpha, ")"),
        hjust = 1.1, vjust = -0.5, size = 4
      )

    if (isTRUE(vert_line)) {
      p <- p + ggplot2::geom_vline(
        xintercept = 0,
        linetype   = "dashed",
        color     = "black"
      )
    }

    if (length(assay_list) > 0) {
      set.seed(seed)
      p <- p + ggrepel::geom_label_repel(
        data = subset(df, Assay_Panel %in% assay_list),
        ggplot2::aes(
          x     = estimate,
          y     = !!rlang::sym(y_var),
          label = Assay
        ),
        inherit.aes   = FALSE,
        color         = "black",
        fill          = "white",
        segment.color = "black",
        min.segment.length = 0,
        label.size    = 0.2,
        box.padding   = 0.7,
        point.padding = 0.5,
        max.overlaps  = nrow(df),
        show.legend   = FALSE
      )
    }

    if (!is.null(shape_var)) {
      p <- p +
        ggplot2::aes(shape = !!rlang::sym(shape_var)) +
        ggplot2::scale_shape_manual(
          values = c(circle = 16, triangle = 17),
          guide  = "none"
        )
    }

    p
  }

  # build & save the full-range volcano
  res_full <- dplyr::mutate(res, y_full = -log10(Adjusted_pval))
  p_full <- build_volcano(
    res_full, "y_full",
    x_var, x_lab,
    NULL,
    title_txt    = "Volcano Plot",
    subtitle_txt = subtitle_txt
  )

  # interactive Plotly
  # With results tooltip and onRender callback
  p_plotly <- plotly::ggplotly(
    p_full +
      ggplot2::aes(
        text = paste0(
          "Protein: ", Assay,
          "<br>Panel: ", Panel,
          paste0("<br>", x_lab, ": "), round(.data[[x_var]], 3),
          "<br>Adjusted p-value: ", sprintf("%.2e", Adjusted_pval)
        ),
        customdata = paste(GeneID, Assay, sep = "|")
      ),
    tooltip = "text"
  ) %>%
    plotly::layout(
      shapes = list(
        list(
          type  = "line",
          x0    = 0,
          x1    = 1,
          xref  = "paper",
          y0    = -log10(alpha),
          y1    = -log10(alpha),
          line  = list(dash = "dash", color = "black")
        )
      ),
      annotations = list(
        list(
          x         = 1,
          y         = -log10(alpha),
          text      = paste0("-log10(", alpha, ")"),
          xref      = "paper",
          yref      = "y",
          xanchor   = "right",
          yanchor   = "bottom",
          showarrow = FALSE,
          font      = list(size = 12)
        )
      )
    ) %>%
    htmlwidgets::onRender(
      "function(el, x) {
        el.on('plotly_click', function(data) {
          var pt    = data.points[0];
          var parts = pt.customdata.split('|');
          var gene  = parts[0];
          var assay = parts[1];
          var url   = 'https://www.proteinatlas.org/' + gene + '-' + assay + '/tissue';
          window.open(url, '_blank');
        });
      }"
    )

  pdf_primary <- file.path(out_dir,
                           paste0(test_type, "_", contrast, "_volcano.pdf"))
  ggplot2::ggsave(pdf_primary, p_full, width = 7, height = 5)

  # optionally save zoomed version
  pdf_secondary <- NULL
  if (export_2nd) {
    res_zoom <- res %>%
      dplyr::mutate(
        y_zoom   = pmin(-log10(Adjusted_pval), export_2nd_y_lim),
        shape_id = dplyr::if_else(
          -log10(Adjusted_pval) > export_2nd_y_lim,
          "triangle", "circle"
        )
      )

    p_zoom <- build_volcano(
      res_zoom, "y_zoom",
        x_var, x_lab,
        "shape_id",
        title_txt    = "Volcano Plot (zoomed)",
        subtitle_txt = subtitle_txt
    ) +
      ggplot2::coord_cartesian(ylim = c(0, export_2nd_y_lim + 0.05))

    pdf_secondary <- file.path(
      out_dir,
      paste0(test_type, "_", contrast, "_volcano_zoom.pdf")
    )
    ggplot2::ggsave(pdf_secondary, p_zoom, width = 7, height = 5)
  }

  invisible(list(
    p_plotly      = p_plotly,
    pdf_primary   = pdf_primary,
    pdf_secondary = pdf_secondary
  ))
}


#' Volcano plots for all contrasts and test types
#'
#' @description
#' Iterates over all contrasts and test types, producing volcano plots
#' as PDF and interactive Plotly objects, embedded in Markdown output.
#'
#' @param res_list Nested list of test results (from \code{run_tests()}).
#' @param se_list Named list of SummarizedExperiment objects.
#' @param config Config list with report tests, alpha, and output paths.
#' @param out_dir Directory path for elisa or pea outputs.
#'
#' @return None. Side effect: writes PDFs and prints Markdown.
#' @export
volcano <- function(res_list, se_list, config, out_dir) {
  for (contrast in names(res_list)) {
    mdcat(paste("###", contrast, "{.tabset}"))

    for (test_type in base::setdiff(names(res_list[[contrast]]), c("emm_int", "params"))) {
      test_header <- get_test_name(test_type)
      mdcat(paste("####", test_header, "Volcano Plot"))
      # Get stats table and parameters used for volcano plot annotations
      res <- res_list[[contrast]][[test_type]]
      parameters <- res_list[[contrast]]$params

      # Critical data checks
      stopifnot(all(c("Assay", "estimate", "Adjusted_pval") %in% names(res)))
      required_params <- c("name",
                           "se_object_name",
                           "column",
                           "treatment",
                           "control",
                           "alpha")
      if (!all(required_params %in% names(parameters))) {
        stop("contrast_list[['", contrast, "']] must contain: ",
             paste(required, collapse = ", "))
      }

      # Make the volcano plot and return it along with its PDF path
      plot_n_path <- make_volcano(res, se_list, contrast, test_type, config, out_dir, parameters)

      mdcat(paste0("- [", test_header, " ", contrast, " Volcano Plot PDF](", plot_n_path$pdf_primary, ")"))

      subchunkify(plot_n_path$p_plotly)
    }
  }
}


#' Boxplots for top significant proteins
#'
#' @description
#' For each contrast and test type, selects the top significant proteins
#' (or those configured in \code{config$label_on_plots}) and produces NPX
#' boxplots with significance thresholds annotated.
#'
#' @param res_list Nested list of test results.
#' @param config Config list with alpha thresholds and output paths.
#' @param out_dir Directory path for elisa or pea outputs.
#' @param dependent_var Column name containing the plotting data
#'
#' @return None. Side effect: writes PDFs and prints Markdown.
#' @export
boxplot <- function(res_list, config, out_dir, dependent_var = "NPX") {
  for (contrast in names(res_list)) {
    mdcat(paste("###", contrast, "{.tabset}"))
    produced <- FALSE
    for (test_type in base::setdiff(names(res_list[[contrast]]), c("emm_int", "params"))) {
      plot_n_path <- make_boxplot(res_list, contrast, test_type, config, out_dir, dependent_var)
      if (!is.null(plot_n_path)) {
        produced <- TRUE
        mdcat(paste("####", get_test_name(test_type), "Boxplot {.tabset}"))
        for (protein in names(plot_n_path$p_list)) {
          mdcat(paste("#####", protein))
          idx <- which(names(plot_n_path$p_list) == protein)
          pdf_path <- plot_n_path$path[idx]
          mdcat(paste0(
            "- [", contrast, " ", get_test_name(test_type),
            " Boxplot for ", protein, " PDF](", pdf_path, ")"
          ))
          subchunkify(plot_n_path$p_list[[protein]])
        }
      }
    }
    if (!produced) {
      mdcat(paste("No boxplot to show for contrast:", contrast))
    }
  }
}


#' @title Internal helper: build boxplot
#' @description
#' Constructs boxplots for top proteins in one contrast/test type.
#' Called by \code{boxplot()}.
#' @noRd
make_boxplot <- function(res_list, contrast, test_type, config, out_dir, dependent_var) {
  long_data <- res_list[[contrast]]$params$df_work
  if (is.null(long_data)) {
    mdcat(paste("No data found for contrast:", contrast, test_type))
    return(invisible(NULL))
  }
  column <- res_list[[contrast]]$params$column
  if (is.null(column)) {
    mdcat(paste("No column found in res_list params for contrast:", contrast))
    return(invisible(NULL))
  }

  one_star   <- config$alpha
  two_stars  <- if (isTRUE(all.equal(config$alpha * 0.2, 0.02))) 0.05 else config$alpha * 0.2
  three_stars <- config$alpha * 0.1

  results_df <- res_list[[contrast]][[test_type]]
  if (is.null(results_df) || nrow(results_df) == 0) {
    mdcat(paste("No", get_test_name(test_type), "results for contrast:", contrast))
    return(invisible(NULL))
  }
  # Uniquely ID proteins
  results_df <- dplyr::mutate(
    results_df,
    Assay_Panel = paste(Assay, Panel, sep = "_")
  )

  if (contrast == "disease_vs_control") {
    sel <- results_df %>%
      dplyr::filter(Assay_Panel %in% config$label_on_plots$disease_vs_control)
  } else {
    sel_up <- results_df %>%
      dplyr::filter(Threshold == "Significant", estimate >  0) %>%
      dplyr::arrange(Adjusted_pval) %>%
      dplyr::slice_head(n = 5)

    sel_dn <- results_df %>%
      dplyr::filter(Threshold == "Significant", estimate <= 0) %>%
      dplyr::arrange(Adjusted_pval) %>%
      dplyr::slice_head(n = 5)

    sel <- dplyr::bind_rows(sel_up, sel_dn)
  }

  if (nrow(sel) == 0) {
    warning(paste("No proteins selected for boxplot in contrast:", contrast))
    return(invisible(NULL))
  }
  # Vector for plotting
  olinkid_list <- sel$OlinkID

  overlay_args <- if (grepl("posthoc", test_type)) {
    list(posthoc_results = results_df)
  } else {
    list(ttest_results   = results_df)
  }

  p_gg_list <- do.call(
    olink_boxplot_2,
    c(
      list(
        df         = long_data,
        variable   = column,
        olinkid_list = olinkid_list,
        number_of_proteins_per_plot = 1,
        one_star   = one_star,
        two_stars  = two_stars,
        three_stars = three_stars
      ),
      overlay_args
    )
  )
  names(p_gg_list) <- sel$Assay_Panel

  if (length(p_gg_list) == 0) {
    warning(paste("Boxplot returned an empty list for contrast:", contrast))
    return(invisible(NULL))
  }

  p_plotly_list <- list()
  pdf_paths <- character()
  for (assay_name in names(p_gg_list)) {
    p_ggplot <- p_gg_list[[assay_name]] +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(color = "black"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::labs(title = "Boxplot", x = "", y = dependent_var) +
      ggplot2::guides(color = "none", fill = "none")

    values <- c(
      paste("*     <", one_star),
      paste("**   <", two_stars),
      paste("*** <", three_stars)
    )
    dummy_sig <- data.frame(
      significance_threshold = factor(values, levels = values),
      x = Inf, y = Inf
    )
    p_ggplot <- p_ggplot +
      ggplot2::geom_point(
        data = dummy_sig,
        ggplot2::aes(x = x, y = y, shape = significance_threshold),
        size = 1, color = "black", show.legend = TRUE, inherit.aes = FALSE
      ) +
      ggplot2::scale_shape_manual(
        name = "Significance Threshold",
        values = setNames(rep(19, 3), values)
      )

    pdf_path <- file.path(
      out_dir,
      paste0(test_type, "_", contrast, "_", assay_name, "_boxplot.pdf")
    )
    ggplot2::ggsave(pdf_path, plot = p_ggplot, device = "pdf", width = 7, height = 5)

    # 1) Add tooltip text to the ggplot object
    p_ggplot <- p_ggplot +
      ggplot2::aes(
        text = paste0(
          "SampleID: ", SampleID,
          "<br>Group: ", column,
          "<br>Value: ", dependent_var
        )
      )

    # 2) Convert with tooltip = "text"
    p_plotly_list[[assay_name]] <- plotly::ggplotly(p_ggplot, tooltip = "text")
        pdf_paths <- c(pdf_paths, pdf_path)
      }

  list(p_list = p_plotly_list, path = pdf_paths)
}


#' Extended Olink NPX/ELISA boxplots
#'
#' @description
#' Variant of \code{OlinkAnalyze::olink_boxplot()} supporting NPX or ELISA,
#' significance thresholds, and ggplotly conversion.
#'
#' @inheritParams make_boxplot
#' @param df Long-format NPX or ELISA data.
#' @param variable Grouping column name.
#' @param olinkid_list OlinkIDs to plot.
#' @param config Config list (theme/scales).
#' @param verbose Logical, print plots if TRUE.
#' @param number_of_proteins_per_plot Integer proteins per facet (default 6).
#' @param posthoc_results Optional posthoc result table.
#' @param ttest_results Optional t-test result table.
#' @param one_star,two_stars,three_stars Numeric thresholds for stars.
#' @param ... Extra args to \code{olink_fill_discrete()}.
#'
#' @return List of ggplot objects, one per batch of proteins.
#' @export
olink_boxplot_2 <- function(
  df,
  variable,
  olinkid_list,
  config,
  verbose = FALSE,
  number_of_proteins_per_plot = 6,
  posthoc_results = NULL,
  ttest_results   = NULL,
  one_star        = 0.05,
  two_stars       = 0.01,
  three_stars     = 0.005,
  show_points        = TRUE,
  point_method       = "jitter",
  point_size         = 1.6,
  point_alpha        = 0.55,
  point_jitter_width = 0.15,
  ...
) {
  point_method <- match.arg(point_method)
  value <- "NPX"
  if (!"NPX" %in% names(df)) {
    value <- "ELISA"
    df$NPX <- df$ELISA
  }

  y_lab <- if (identical(value, "ELISA")) "ELISA log10 pg/ml" else "NPX"

  npx_flag <- OlinkAnalyze:::npxCheck(df)
  df_clean <- df %>%
    dplyr::filter(stringr::str_detect(OlinkID, "OID[0-9]{5}"), !OlinkID %in% npx_flag$all_nas)

  if (!variable %in% names(df_clean)) {
    stop("Column '", variable, "' not found in data")
  }
  df_clean[[variable]] <- factor(df_clean[[variable]], levels = unique(df_clean[[variable]]))
  x_levels <- levels(df_clean[[variable]])

  total  <- length(olinkid_list)
  starts <- seq(1, total, by = number_of_proteins_per_plot)
  plots  <- vector("list", length(starts))

  for (i in seq_along(starts)) {
    assays  <- olinkid_list[starts[i]:min(starts[i] + number_of_proteins_per_plot - 1, total)]
    df_plot <- df_clean %>%
      dplyr::filter(OlinkID %in% assays) %>%
      dplyr::mutate(Name_OID = factor(paste(Assay, OlinkID), levels = unique(paste(Assay, OlinkID))))

    # Build plot with points
    p <- ggplot2::ggplot(df_plot, ggplot2::aes_string(x = variable, y = value))

    if (isTRUE(show_points)) {
      if (identical(point_method, "beeswarm") && requireNamespace("ggbeeswarm", quietly = TRUE)) {
        p <- p +
          ggbeeswarm::geom_quasirandom(
            ggplot2::aes_string(color = variable),
            width   = 0.15,  # gentle spread
            alpha   = point_alpha,
            size    = point_size,
            na.rm   = TRUE
          )
      } else {
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes_string(color = variable),
            position = ggplot2::position_jitter(width = point_jitter_width, height = 0),
            alpha    = point_alpha,
            size     = point_size,
            na.rm    = TRUE
          )
      }
    }

    p <- p +
      ggplot2::geom_boxplot(ggplot2::aes_string(fill = variable), outlier.shape = NA) +
      ggplot2::facet_wrap(~Name_OID, scales = "free") +
      OlinkAnalyze::set_plot_theme() +
      OlinkAnalyze::olink_fill_discrete(...) +
      ggplot2::guides(color = "none") +  # hide duplicate color legend for points
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
        legend.text  = ggplot2::element_text(size = 13)
      )

    # (Optional) significance stars
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    res <- ttest_results %||% posthoc_results
    if (!is.null(res)) {
      thr <- res %>%
        dplyr::filter(OlinkID %in% assays) %>%
        dplyr::left_join(df_plot %>% dplyr::distinct(OlinkID, Name_OID), by = "OlinkID") %>%
        dplyr::mutate(
          Threshold = dplyr::case_when(
            Adjusted_pval <= three_stars ~ "***",
            Adjusted_pval <= two_stars   ~ "**",
            Adjusted_pval <= one_star    ~ "*",
            TRUE                         ~ NA_character_
          )
        ) %>%
        dplyr::filter(!is.na(Threshold))

      if (nrow(thr) > 0) {
        stats_df <- df_plot %>%
          dplyr::group_by(Name_OID) %>%
          dplyr::summarise(maxNPX = max(NPX), rangeNPX = diff(range(NPX)), .groups = "drop")

        line_data <- thr %>%
          dplyr::left_join(stats_df, by = "Name_OID") %>%
          tidyr::crossing(x_vals = x_levels) %>%
          dplyr::mutate(
            x_num    = match(x_vals, x_levels),
            y_anchor = maxNPX + 0.1 * rangeNPX
          ) %>%
          dplyr::group_by(Name_OID, Threshold) %>%
          dplyr::mutate(x_mid = mean(x_num)) %>%
          dplyr::ungroup()

        p <- p +
          ggplot2::geom_line(
            data = line_data,
            ggplot2::aes(x = x_num, y = y_anchor, group = Name_OID),
            inherit.aes = FALSE
          ) +
          ggplot2::geom_text(
            data = dplyr::distinct(line_data, Name_OID, Threshold, x_mid, y_anchor, rangeNPX),
            ggplot2::aes(x = x_mid, y = y_anchor + 0.05 * rangeNPX, label = Threshold),
            inherit.aes = FALSE
          ) +
          ggplot2::labs(y = y_lab)
      }
    }

    plots[[i]] <- p
    if (verbose) print(p)
  }

  invisible(plots)
}


#' Dotplots of enrichment results
#'
#' @description
#' For each contrast, test type, method, and ontology/direction, creates
#' dotplots of enriched terms (ORA or GSEA) if enough terms are significant.
#' Outputs PDFs and interactive Plotly objects.
#'
#' @param enrich_list Nested list of enrichment results from \code{enrich()}.
#' @param config Config list with alpha, min_terms, output paths.
#'
#' @return None. Side effect: writes PDFs and prints Markdown.
#' @export
dotplot <- function(enrich_list, config) {
  # minimum-term cutoffs
  default_min <- list(ORA = 10, GSEA = 5)
  min_terms   <- config$enrichment$min_terms %||% default_min

  if (is.null(names(enrich_list))) {
    mdcat("To plot enrichment results, add 'GSEA' and/or 'ORA' to config$enrichment$methods.")
    return(invisible(NULL))
  }

  for (contrast in names(enrich_list)) {
    mdcat(paste0("### ", contrast, " {.tabset}"))

    for (method in names(enrich_list[[contrast]])) {
      mdcat(paste0("#### ", method, " {.tabset}"))
      meth_tbls <- enrich_list[[contrast]][[method]]
      if (is.null(meth_tbls)) {
        mdcat("*No results for method ", method, ".*")
        next
      }

      for (test_type in names(meth_tbls)) {
        header <- get_test_name(test_type)
        mdcat(paste0("##### ", header, " {.tabset}"))
        tt_tbls <- meth_tbls[[test_type]]
        if (is.null(tt_tbls)) {
          mdcat("*No results for ", header, ".*")
          next
        }

        for (direction in names(tt_tbls)) {
          mdcat(paste0("###### ", direction, " {.tabset}"))
          dir_tbls <- tt_tbls[[direction]]
          if (is.null(dir_tbls)) {
            mdcat("*No ", direction, " results.*")
            next
          }

          for (ontology in names(dir_tbls)) {
            mdcat(paste0("####### ", ontology))
            df <- dir_tbls[[ontology]]
            if (is.null(df) || nrow(df) == 0) {
              mdcat(sprintf("*No %s/%s/%s/%s results.*",
                            method, header, direction, ontology))
              next
            }
            # count significant terms
            cutoff  <- min_terms[[method]] %||% default_min[[method]]
            num_sig <- sum(df$p.adjust < config$alpha, na.rm = TRUE)
            if (num_sig < cutoff) {
              mdcat(sprintf(
                "*Insufficient (< %d) terms (%d) for %s/%s/%s/%s; skipping.*",
                cutoff, num_sig, method, header, direction, ontology
              ))
              next
            }
            # make and embed plot
            if (method == "GSEA") {
              pp <- make_dotplot_gsea(df, contrast, test_type, direction, ontology, config)
            } else {
              pp <- make_dotplot_ora(df, contrast, test_type, direction, ontology, config)
            }
            mdcat(sprintf(
              "- [%s %s %s %s %s Dotplot PDF](%s)",
              method, header, direction, ontology, contrast, pp$path
            ))
            subchunkify(pp$p_plotly)
          }
        }
      }
    }
  }
  invisible(NULL)
}


#' @title Internal helper: ORA dotplot
#' @description
#' Builds a dotplot for top ORA terms (up to 20).
#' @noRd
make_dotplot_ora <- function(df, contrast, test_type, direction, ontology, config) {
  df_top <- df %>% head(20) %>%
    mutate(
      GeneRatioNum = as.numeric(sub("/.*","",GeneRatio)) /
                     as.numeric(sub(".*/","",GeneRatio))
    )
  p_gg <- ggplot2::ggplot(df_top, ggplot2::aes(
    x = GeneRatioNum,
    y = stats::reorder(ID, GeneRatioNum),
    size = Count,
    color = -log10(p.adjust)
  )) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste(contrast, get_test_name(test_type),
                    toupper(direction), ontology, "ORA Dotplot"),
      x = "GeneRatio",
      y = "Gene-set ID",
      color = "-log10(p.adjust)",
      size = "Count"
    )

  fname <- paste(test_type, contrast, direction, ontology, "ORA_dotplot.pdf", sep = "_")
  path  <- file.path(config$output_paths$functional_enrichment$plots, fname)
  ggplot2::ggsave(path, p_gg, device = "pdf", width = 7, height = 5)

  p_plotly <- plotly::ggplotly(
    p_gg + ggplot2::aes(
      text = paste0(
        "Term: ", Description,
        "<br>GeneRatio: ", round(GeneRatioNum, 3),
        "<br>p.adjust: ", sprintf("%.2e", p.adjust),
        "<br>Count: ", Count
      )
    ),
    tooltip = "text"
  )
  list(p_plotly = p_plotly, path = path)
}


#' @title Internal helper: GSEA dotplot
#' @description
#' Builds a dotplot for top GSEA terms (up to 20).
#' @noRd
make_dotplot_gsea <- function(df, contrast, test_type, direction, ontology, config) {
  df_top <- df %>% head(20)
  p_gg <- ggplot2::ggplot(df_top, ggplot2::aes(
    x = NES,
    y = stats::reorder(ID, NES),
    size = setSize,
    color = -log10(p.adjust)
  )) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste(contrast, get_test_name(test_type),
                    toupper(direction), ontology, "GSEA Dotplot"),
      x = "Normalized Enrichment Score (NES)",
      y = "Gene-set ID",
      color = "-log10(p.adjust)",
      size = "Gene set size"
    )

  fname <- paste(test_type, contrast, direction, ontology, "GSEA_dotplot.pdf", sep = "_")
  path  <- file.path(config$output_paths$functional_enrichment$plots, fname)
  ggplot2::ggsave(path, p_gg, device = "pdf", width = 7, height = 5)

  p_plotly <- plotly::ggplotly(
    p_gg + ggplot2::aes(
      text = paste0(
        "Term: ", Description,
        "<br>NES: ", round(NES, 3),
        "<br>p.adjust: ", sprintf("%.2e", p.adjust),
        "<br>Size: ", setSize
      )
    ),
    tooltip = "text"
  )
  list(p_plotly = p_plotly, path = path)
}


#' Pathway-level heatmap
#'
#' @description
#' Wrapper around \code{OlinkAnalyze::olink_pathway_heatmap}, saving a PDF
#' and returning interactive Plotly.
#'
#' @param df Data frame or SummarizedExperiment.
#' @param ... Passed to \code{olink_pathway_heatmap()}.
#'
#' @return Plotly object.
#' @export
pathway_heatmap <- function(df, ...) {
  mdcat("### Olink Pathway Heatmap")

  # Call the OlinkAnalyze function
  p_ggplot <- OlinkAnalyze::olink_pathway_heatmap(data, ...)

  # Save to PDF
  path <- file.path(config$output_paths$functional_enrichment$plots, "olink_pathway_heatmap.pdf")
  ggplot2::ggsave(path, p_ggplot, device = "pdf", width = 7, height = 5)
  mdcat(paste0("- [Pathway Heatmap PDF](", path, ")"))

  # Convert to plotly for interactivity
  p_plotly <- plotly::ggplotly(p_ggplot)

  return(p_plotly)
}


#' Pathway-level visualization
#'
#' @description
#' Wrapper around \code{OlinkAnalyze::olink_pathway_visualization}, saving a
#' PDF and returning Plotly.
#'
#' @param df Data frame or SummarizedExperiment.
#' @param ... Passed to \code{olink_pathway_visualization()}.
#'
#' @return Plotly object.
#' @export
pathway_visualization <- function(df, ...) {
  mdcat("### Olink Pathway Visualization")

  p_ggplot <- OlinkAnalyze::olink_pathway_visualization(data, ...)

  path <- file.path(config$output_paths$functional_enrichment$plots, "olink_pathway_visualization.pdf")
  ggplot2::ggsave(path, p_ggplot, device = "pdf", width = 7, height = 5)
  mdcat(paste0("- [Pathway Visualization PDF](", path, ")"))

  p_plotly <- plotly::ggplotly(p_ggplot)
  return(p_plotly)
}


#' NPX quality-control plots
#'
#' @description
#' Wrapper for \code{OlinkAnalyze::olink_qc_plot}. Produces IQR vs median
#' plots for each config color column, with numeric columns shown as
#' continuous gradients. Saves PDFs and embeds ggplot/plotly in Markdown.
#'
#' @param long_data Long-format NPX data.
#' @param config Config list with qc color columns
#' @param out_dir Output path from config.
#' @param path_prefix String prefix for output file names.
#' @param ... Extra args passed to \code{olink_qc_plot()}.
#'
#' @return Invisibly, the last plot produced.
#' @export
qc <- function(long_data, config, out_dir, path_prefix, ...) {
  mdcat("### QC Plots {.tabset}")

  color_cols <- unique(c(
    "QC_Warning",
    config$qc_plot_color_column$plotly,
    config$qc_plot_color_column$ggplot
  ))

  last_plot <- NULL

  for (color_column in color_cols) {
    mdcat(paste0("#### ", color_column, " {.tabset}"))
    message("Making QC plot for ", color_column)

    # If the column doesn't exist in long_data, skip with a warning
    if (!color_column %in% colnames(long_data)) {
      warning("Skipping '", color_column, "' because it is not in long_data.")
      next
    }

    # Check numeric vs discrete
    is_num <- is.numeric(long_data[[color_column]])

    if (!is_num) {
      # Discrete: coerce to factor (now that we know the column exists)
      long_data[[color_column]] <- as.factor(long_data[[color_column]])
    }

    # Build base QC ggplot
    p_gg <- OlinkAnalyze::olink_qc_plot(
      long_data,
      color_g = !!rlang::sym(color_column),
      ...
    ) +
      ggplot2::labs(title = paste(color_column, "IQR vs. Median Sample QC Plot"))

    # Apply continuous or discrete color scale
    if (is_num) {
      p_gg <- p_gg +
        ggplot2::scale_color_gradient(
          low  = "lightblue",
          high = "darkblue",
          name = color_column
        )
    } else {
      p_gg <- p_gg +
        ggplot2::scale_color_discrete(name = color_column)
    }

    # Ensure QC_Warning == "Warning" points are drawn on top of Pass points
    # Use p_gg$data so we have access to the derived columns like sample_median and IQR
    if (color_column == "QC_Warning" && "QC_Warning" %in% names(p_gg$data)) {
      warning_data <- p_gg$data[p_gg$data$QC_Warning == "Warning", , drop = FALSE]

      p_gg <- p_gg +
        ggplot2::geom_point(
          data = warning_data,
          inherit.aes = TRUE,   # reuse the x/y/color aesthetics from the base plot
          size = 2.5,           # optionally make warnings a bit bigger; remove if undesired
          color = "#00bec3"
        )
    }

    # Save static PDF
    pdf_path <- file.path(
      out_dir,
      paste0(path_prefix, "_", color_column, "_qc_plot.pdf")
    )
    ggplot2::ggsave(
      filename = pdf_path,
      plot     = p_gg,
      device   = "pdf",
      width    = 12,
      height   = 12
    )
    mdcat(paste0("- [", color_column, " QC Plot PDF](", pdf_path, ")"))

    interactive <- color_column %in% config$qc_plot_color_column$plotly

    if (!interactive) {
      subchunkify(p_gg, fig.width = 10, fig.height = 10)
      last_plot <- p_gg
    } else {
      # Add hover text aesthetic
      p_gg <- p_gg +
        ggplot2::aes(
          text = paste0(
            "Sample: ", SampleID,
            "<br>Panel: ", Panel,
            "<br>IQR: ", IQR,
            "<br>Median: ", sample_median,
            "<br>", color_column, ": ", !!rlang::sym(color_column),
            "<br>Median Low: ", median_low,
            "<br>Median High: ", median_high,
            "<br>IQR Low: ", iqr_low,
            "<br>IQR High: ", iqr_high,
            "<br>Outlier: ", Outlier
          )
        )

      for (panel in unique(p_gg$data$Panel)) {
        mdcat(paste0("##### Panel: ", panel))
        panel_data <- p_gg$data %>% dplyr::filter(Panel == panel)
        panel_gg <- ggplot2::`%+%`(p_gg, panel_data)
        panel_pl <- plotly::ggplotly(panel_gg, tooltip = "text")
        subchunkify(panel_pl)
        last_plot <- panel_pl
      }
    }
  }

  invisible(last_plot)
}


#' PCA plots
#'
#' @description
#' Wrapper for \code{OlinkAnalyze::olink_pca_plot}, looping over color
#' columns. Saves PDFs, embeds ggplot or plotly plots. Optionally labels
#' selected samples.
#'
#' @inheritParams qc
#' @param label Optional vector of SampleIDs to label.
#'
#' @return None (plots are saved and displayed).
#' @export
pca <- function(long_data,
                config,
                out_dir,
                path_prefix,
                label = NULL) {

  mdcat("### PCA Plots {.tabset}")

  color_cols <- unique(c(
    "QC_Warning",
    config$qc_plot_color_column$plotly,
    config$qc_plot_color_column$ggplot
  ))

  for (color_column in color_cols) {
    mdcat(paste0("#### ", color_column, " {.tabset}"))

    # Check type of color column
    is_num <- is.numeric(long_data[[color_column]])

    if (!is_num) {
      # If not numeric, convert to factor for discrete coloring
      long_data[[color_column]] <- as.factor(long_data[[color_column]])
    }

    # Build base PCA ggplot
    p_gg <- OlinkAnalyze::olink_pca_plot(
      long_data,
      color_g        = color_column,
      outlierDefX    = 3,
      outlierDefY    = 3,
      outlierLines   = TRUE,
      label_outliers = FALSE,
      quiet          = TRUE
    )[[1]] +
      ggplot2::labs(title = paste("PCA colored by", color_column))

      # Remove ALL internal theme settings so ggplot falls back to the global
      # theme setting defined in helpers.R
      p_gg$theme <- ggplot2::theme()

    # Apply continuous or discrete scale
    if (is_num) {
      p_gg <- p_gg +
        ggplot2::scale_color_gradient(
          low  = "lightblue",
          high = "darkblue",
          name = color_column
        )
    } else {
      # default discrete scale is already applied by olink_pca_plot(aes(color = ...))
      # but ensure legend title is blank or the column name
      p_gg <- p_gg +
        ggplot2::scale_color_discrete(name = color_column)
    }

    # If label vector is provided, add geom_label_repel
    if (!is.null(label) && length(label) > 0) {
      label_data <- p_gg$data %>%
        dplyr::filter(SampleID %in% label)

      if (nrow(label_data) > 0) {
        p_gg <- p_gg +
          ggrepel::geom_label_repel(
            data  = label_data,
            ggplot2::aes(x = PCX, y = PCY, label = SampleID),
            inherit.aes        = FALSE,
            color              = "black",
            fill               = "white",
            segment.color      = "black",
            force              = 1,
            min.segment.length = 0,
            label.size         = 0.25,
            box.padding        = 0.7,
            point.padding      = 0.25,
            max.overlaps       = nrow(label_data),
            show.legend        = FALSE
          )
      }
    }

    # Save static PDF
    pdf_path <- file.path(
      out_dir,
      paste0(path_prefix, "_", color_column, "_pca_plot.pdf")
    )
    ggplot2::ggsave(
      filename = pdf_path,
      plot     = p_gg,
      device   = "pdf",
      width    = 7,
      height   = 5
    )
    mdcat(paste0("- [", color_column, " PCA Plot PDF](", pdf_path, ")"))

    # Render interactive or static
    if (color_column %in% config$qc_plot_color_column$plotly) {
      color_df <- long_data %>%
        dplyr::select(SampleID, !!rlang::sym(color_column)) %>%
        dplyr::distinct()

      p_gg$data <- p_gg$data %>%
        dplyr::inner_join(color_df, by = "SampleID")

      # Re-add labels if needed in the interactive version
      if (!is.null(label) && length(label) > 0) {
        label_data <- p_gg$data %>%
          dplyr::filter(SampleID %in% label)
        if (nrow(label_data) > 0) {
          p_gg <- p_gg +
            ggrepel::geom_label_repel(
              data  = label_data,
              ggplot2::aes(x = PCX, y = PCY, label = SampleID),
              inherit.aes        = FALSE,
              color              = "black",
              fill               = "white",
              segment.color      = "black",
              min.segment.length = 0,
              label.size         = 0.25,
              box.padding        = 0.7,
              point.padding      = 0.25,
              max.overlaps       = nrow(label_data),
              show.legend        = FALSE
            )
        }
      }

      p_gg <- p_gg +
        ggplot2::aes(
          text = paste0(
            "Sample: ", SampleID,
            "<br>", color_column, ": ", !!rlang::sym(color_column),
            "<br>PCX: ", PCX,
            "<br>PCY: ", PCY
          )
        )

      subchunkify(plotly::ggplotly(p_gg, tooltip = "text"))
    } else {
      subchunkify(p_gg)
    }
  }

  invisible(NULL)
}


#' UMAP plots
#'
#' @description
#' Wrapper for \code{OlinkAnalyze::olink_umap_plot}, looping over color
#' columns. Saves PDFs, embeds ggplot or plotly plots.
#'
#' @inheritParams qc
#' @param seed Random seed for UMAP reproducibility.
#'
#' @return None (plots are saved and displayed).
#' @export
umap <- function(long_data, config, out_dir, path_prefix, seed = 1) {

  mdcat("### UMAP Plots {.tabset}")

  color_cols <- unique(c(
    "QC_Warning",
    config$qc_plot_color_column$plotly,
    config$qc_plot_color_column$ggplot
  ))

  for (color_column in color_cols) {

    mdcat(paste0("#### ", color_column, " {.tabset}"))
    long_data[[color_column]] <- as.factor(long_data[[color_column]])

    # Set seed equivalent for umap::umap
    conf <- umap::umap.defaults
    conf$random_state <- seed

    p_gg <- OlinkAnalyze::olink_umap_plot(
      long_data,
      color_g = color_column,
      quiet   = TRUE,
      config = conf
    )[[1]] +
      ggplot2::labs(title = paste("UMAP colored by", color_column))

      # Remove ALL internal theme settings so ggplot falls back to the global
      # theme setting defined in helpers.R
      p_gg$theme <- ggplot2::theme()

    pdf_path <- file.path(
      out_dir,
      paste0(path_prefix, "_", color_column, "_umap_plot.pdf")
    )
    ggplot2::ggsave(pdf_path, p_gg, device = "pdf", width = 7, height = 5)
    mdcat(paste0("- [", color_column, " UMAP Plot PDF](", pdf_path, ")"))

    if (color_column %in% config$qc_plot_color_column$plotly) {

      p_gg$data$SampleID <- rownames(p_gg$data)
      color_df <- long_data %>%
        dplyr::select(SampleID, !!rlang::sym(color_column)) %>%
        dplyr::distinct()
      p_gg$data <- p_gg$data %>%
        dplyr::inner_join(color_df, by = "SampleID")

      p_gg <- p_gg + ggplot2::aes(
        text = paste0(
          "Sample: ", SampleID,
          "<br>", color_column, ": ", !!rlang::sym(color_column),
          "<br>UMAPX: ", umapX,
          "<br>UMAPY: ", umapY
        )
      )

      subchunkify(plotly::ggplotly(p_gg, tooltip = "text"))

    } else {

      subchunkify(p_gg)
    }
  }

  invisible(NULL)
}


#' Protein heatmap
#'
#' @description
#' Wrapper for \code{OlinkAnalyze::olink_heatmap_plot}, saving a PDF and
#' embedding ggplot or plotly output.
#'
#' @param long_data Long-format NPX data.
#' @param config Config list with output paths.
#' @param ... Extra args to \code{olink_heatmap_plot()}.
#'
#' @return None.
#' @export
protein_heatmap <- function(long_data, config, ...) {
  # Top-level tabset
  mdcat("### Protein Heatmap {.tabset}")
  # Define the color columns to loop over
  color_cols <- c("QC_Warning", config$qc_plot_color_column)
  # We'll store plotly outputs in a list if desired
  plotly_list <- list()
  # Pass color_g to the heatmap function maybe as row or column variable?
  p_ggplot <- OlinkAnalyze::olink_heatmap_plot(
    long_data,
    ...
  )

  # Save as PDF
  pdf_filename <- "protein_heatmap_plot.pdf"
  pdf_path     <- file.path(config$output_paths$differential_expression$qc, pdf_filename)
  ggplot2::ggsave(pdf_path, p_ggplot, device = "pdf", width = 7, height = 5)
  # Provide a link to the PDF
  mdcat(paste0("- [Protein Heatmap Plot PDF](", pdf_path, ")"))

  if (isFALSE(plotly)) {
    subchunkify(p_ggplot)
  } else {
    # Convert to Plotly
    p_plotly <- plotly::ggplotly(p_ggplot)
    # Print the plotly plot
    subchunkify(p_plotly)
  }
}


#' Bridgeability plot
#'
#' @description
#' Wrapper for \code{OlinkAnalyze::olink_bridgeability_plot}, saving a PDF
#' and returning interactive Plotly.
#'
#' @param long_data NPX data.
#' @param ... Extra args passed to \code{olink_bridgeability_plot()}.
#'
#' @return Plotly object.
#' @export
bridgeability <- function(long_data, ...) {
  mdcat("### NPX Data Bridgeability Plot")

  p_ggplot <- OlinkAnalyze::olink_bridgeability_plot(long_data, ...)

  path <- file.path(config$output_paths$differential_expression$qc, "olink_bridgeability_plot.pdf")
  ggplot2::ggsave(path, p_ggplot, device = "pdf", width = 7, height = 5)
  mdcat(paste0("- [Bridgeability Plot PDF](", path, ")"))

  p_plotly <- plotly::ggplotly(p_ggplot)
  return(p_plotly)
}


#' UpSet plots across contrasts
#'
#' @description
#' For each shared test type, builds UpSet plots of significant proteins
#' across contrasts in up/down/changed directions. Saves PDFs and Excel
#' workbooks with included sets.
#'
#' @param res_list Nested results list.
#' @param config Config with alpha, report_test, and output paths.
#'
#' @return Invisibly, list with `excel` and `pdf` file paths.
#' @export
upset <- function(res_list, config, out_dir) {

  exported <- list(
    excel = list(),
    pdf   = list()
  )

  # Warn and exit if upset_df has 1 column (at least 2 contrasts are required for the upset plot)
  if (length(names(res_list)) < 2) {
    warning(crayon::red("There needs to be at least 2 contrasts to compare in the upset plot"))
    mdcat("At least 2 contrasts are required to make an upset plot")
    return(invisible(NULL))
  }

  # Get the test names shared between all contrasts
  test_names_per_contrast <- lapply(res_list, names)
  common_tests <- Reduce(intersect, test_names_per_contrast)
  common_tests <- base::setdiff(common_tests, c("emm_int", "params"))

  if (length(common_tests) == 0) {
    mdcat("There are no tests shared between two or more contrasts")
    return(invisible(NULL))
  }

  # Loop over each test_type
  for (test_type in common_tests) {
    test_hdr <- get_test_name(test_type)
    mdcat(paste0("### ", test_hdr, " {.tabset}"))

    # Precompute workbook path for this test_type
    wb_path <- file.path(
      out_dir,
      paste0("upset_", test_type, ".xlsx")
    )

    # Build, for this test_type, the sets of proteins per direction
    sets_by_direction <- lapply(
      c("up", "down", "changed"),
      function(direction) {
        lapply(
          names(res_list),
          function(contrast) get_sig(res_list, contrast, direction, test_type)
        ) %>% setNames(names(res_list))
      }
    ) %>% setNames(c("up", "down", "changed"))

    # Prepare sheets storage for Excel
    wb_sheets <- list()

    # Inner tabs: each direction
    for (direction in names(sets_by_direction)) {
      mdcat(paste0("#### ", direction))

      sets <- sets_by_direction[[direction]]
      if (sum(vapply(sets, nrow, 1L) > 0) < 2) {
        mdcat("Not enough sets with significant proteins to create an UpSet plot.")
        next
      }

      upset_df <- from_list_with_names(
        lapply(sets, function(df) df$Assay)   # convert tibbles -> vectors
      )

      wb_sheets[[direction]] <- tibble::rownames_to_column(
        upset_df, var = "Assay_Panel")

      pdf_file <- file.path(out_dir,
                            paste0("upset_", test_type, "_", direction, ".pdf"))
      grDevices::pdf(pdf_file, width = 7, height = 5)
      UpSetR::upset(upset_df, order.by = "freq")
      grDevices::dev.off()
      exported$pdf_path[[paste(test_type, direction, sep = "_")]] <- pdf_file

      mdcat(paste0("- [Excel Workbook](", wb_path, ")  "))
      mdcat(paste0("- [Plot PDF](", pdf_file, ")\n\n"))

      subchunkify(UpSetR::upset(upset_df, order.by = "freq"))
    }

    # save workbook only if something was written
    if (length(wb_sheets)) {
      wb <- openxlsx::createWorkbook()
      for (sheet_name in names(wb_sheets)) {
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet_name, wb_sheets[[sheet_name]], withFilter = TRUE)
        openxlsx::setColWidths(wb, sheet = sheet_name,
                              cols  = seq_len(ncol(wb_sheets[[sheet_name]])),
                              widths = "auto")
      }
      openxlsx::saveWorkbook(wb, wb_path, overwrite = TRUE)
    }
    exported$excel_path[[test_type]] <- wb_path
  }

  invisible(exported)
}


#' Estimate scatter-plot between contrasts
#'
#' @description
#' Creates scatter-plot comparing estimates between two contrasts, labeling
#' significant overlaps (optionally concordant or discordant). Saves PDF,
#' optionally exports table, and can return plotly.
#'
#' @param res Nested list of results.
#' @param config Config with output paths.
#' @param contrasts Character vector of length 2, the contrasts to compare.
#' @param label_concord Logical, label discordant significant overlaps.
#' @param label_sig_both Logical, label all significant overlaps.
#' @param plotly Logical, return plotly interactive plot.
#' @param ... Additional options (y.lim, label_n, export_table, etc.).
#'
#' @return Invisibly NULL or the plotting table (if return_table = TRUE).
#' @export
estimate_plot <- function(res,
                          config,
                          out_dir,
                          contrasts,
                          label_concord   = TRUE,
                          label_sig_both  = FALSE,
                          plotly          = FALSE,
                          color           = "Threshold",
                          y.lim           = c(-2, 2),
                          label_n         = 10,
                          do_plot         = TRUE,
                          export_table    = FALSE,
                          export_path     = NULL,
                          return_table    = FALSE) {

  # Sanity checks
  stopifnot(
    length(contrasts) == 2,
    all(contrasts %in% names(res)),
    is.list(res[[contrasts[1]]]),
    is.list(res[[contrasts[2]]]),
    length(y.lim) == 2,
    is.numeric(y.lim)
  )
  if (isTRUE(label_concord) && isTRUE(label_sig_both)) {
    stop("Only one of label_concord or label_sig_both may be TRUE.")
  }

  alpha <- config$alpha

  # Helper: extract estimate + significance for one contrast
  pull_estimate <- function(ctr) {
    an_tbl      <- res[[ctr]]$anova
    effect_name <- res[[ctr]]$params$column
    an_sub      <- dplyr::filter(an_tbl, Effect == effect_name)

    if (!"estimate" %in% names(an_sub)) {
      stop("Contrast '", ctr, "' ANOVA does not have an 'estimate' column.")
    }

    cols_present <- intersect(
      c("Assay","GeneID","Panel","OlinkID","UniProt","estimate", color, "p.value","Adjusted_pval"),
      names(an_sub)
    )

    tmp <- an_sub %>%
      dplyr::select(dplyr::all_of(cols_present)) %>%
      dplyr::rename(
        value   = estimate,
        sig_col = dplyr::all_of(color)
      )

    # Coerce alpha and p-values to numeric
    alpha_num <- suppressWarnings(as.numeric(alpha))
    if (is.na(alpha_num)) stop("alpha must be numeric and not NA")

    pval_num <- if ("p.value" %in% names(tmp))        suppressWarnings(as.numeric(tmp[["p.value"]]))        else rep(NA_real_, nrow(tmp))
    adjp_num <- if ("Adjusted_pval" %in% names(tmp))   suppressWarnings(as.numeric(tmp[["Adjusted_pval"]]))   else rep(NA_real_, nrow(tmp))

    # Prefer adjusted p-values for significance (aligns with your display columns)
    # Fallback to the Threshold/`color` column ONLY if adjp is missing
    if (!all(is.na(adjp_num))) {
      sig_logical <- adjp_num < alpha_num
    } else if ("sig_col" %in% names(tmp)) {
      if (is.logical(tmp[["sig_col"]])) {
        sig_logical <- tmp[["sig_col"]]
      } else {
        sc <- as.character(tmp[["sig_col"]])
        sig_logical <- sc %in% c("Significant", "TRUE", "True", "true", "1")
      }
    } else {
      sig_logical <- rep(FALSE, nrow(tmp))
    }

    tmp %>%
      dplyr::mutate(
        sig      = sig_logical,
        pval     = pval_num,
        adjp     = adjp_num,
        Contrast = ctr
      ) %>%
      dplyr::select(-dplyr::any_of("sig_col"))
  }

  df1 <- pull_estimate(contrasts[1])
  df2 <- pull_estimate(contrasts[2])

  # Join the two contrasts
  joined <- dplyr::full_join(
    df1, df2,
    by = c("Assay", "Panel"),
    suffix = c(".A", ".B")
  )

  # Build status factor
  base_levels <- c(
    "Not significant",
    sprintf("Significant in %s only", contrasts[1]),
    sprintf("Significant in %s only", contrasts[2]),
    "Significant in both"
  )

  joined <- joined %>%
    dplyr::mutate(
      sigA = dplyr::coalesce(.data[["sig.A"]], FALSE),
      sigB = dplyr::coalesce(.data[["sig.B"]], FALSE),
      status_char = dplyr::case_when(
        sigA &  sigB ~ "Significant in both",
        sigA & !sigB ~ sprintf("Significant in %s only", contrasts[1]),
       !sigA &  sigB ~ sprintf("Significant in %s only", contrasts[2]),
        TRUE         ~ "Not significant"
      ),
      status = factor(
        status_char,
        levels = c(
          "Not significant",
          sprintf("Significant in %s only", contrasts[1]),
          sprintf("Significant in %s only", contrasts[2]),
          "Significant in both"
        )
      )
    ) %>%
    dplyr::select(-status_char) %>%
    # keep sigA/sigB if you find them handy; otherwise drop:
    # dplyr::select(-sigA, -sigB) %>%
    dplyr::mutate(
      Assay_Panel = paste(Assay, Panel, sep = "_"),
      GeneID  = dplyr::coalesce(.data[["GeneID.A"]],  .data[["GeneID.B"]]),
      OlinkID = dplyr::coalesce(.data[["OlinkID.A"]], .data[["OlinkID.B"]]),
      UniProt = dplyr::coalesce(.data[["UniProt.A"]], .data[["UniProt.B"]])
    )

  # Clip extreme y values & choose shapes
  joined <- joined %>%
    dplyr::mutate(
      y_plot = dplyr::case_when(
        value.B < y.lim[1] ~ y.lim[1],
        value.B > y.lim[2] ~ y.lim[2],
        TRUE               ~ value.B
      ),
      shape = dplyr::case_when(
        value.B < y.lim[1] ~ 6,   # upside-down triangle
        value.B > y.lim[2] ~ 2,   # upright triangle
        TRUE               ~ 16   # circle
      )
    ) %>%
    dplyr::arrange(status)

  # ---------------------------
  # Construct a clean plotting/data-export table
  # ---------------------------
  estA_name <- paste0("estimate_", contrasts[1])
  estB_name <- paste0("estimate_", contrasts[2])
  sigA_name <- paste0("sig_",       contrasts[1])
  sigB_name <- paste0("sig_",       contrasts[2])
  adjA_name <- paste0("Adjusted_pval_", contrasts[1])
  adjB_name <- paste0("Adjusted_pval_", contrasts[2])
  pA_name   <- paste0("p.value_",       contrasts[1])
  pB_name   <- paste0("p.value_",       contrasts[2])

  plot_tbl <- joined %>%
    dplyr::transmute(
      Assay, GeneID, Panel, Assay_Panel, OlinkID, UniProt,
      !!estA_name := .data[["value.A"]],
      !!estB_name := .data[["value.B"]],
      !!sigA_name := .data[["sig.A"]],
      !!sigB_name := .data[["sig.B"]],
      !!adjA_name := .data[["adjp.A"]],
      !!adjB_name := .data[["adjp.B"]],
      !!pA_name   := .data[["pval.A"]],
      !!pB_name   := .data[["pval.B"]],
      status,
      significant_in_both =
        dplyr::coalesce(.data[["sig.A"]], FALSE) & dplyr::coalesce(.data[["sig.B"]], FALSE),
      reversal = (.data[["value.A"]] * .data[["value.B"]]) < 0,
      y_plot,
      shape
    )

  # Optional export
  if (isTRUE(export_table)) {
    if (is.null(export_path)) {
      export_path <- file.path(
        file.path(dirname(out_dir), "results"),
        paste0(paste(contrasts, collapse = "_"), "_estimate_plotting_table.tsv")
      )
    }
    dir.create(dirname(export_path), recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(plot_tbl, export_path)
    } else {
      utils::write.table(plot_tbl, export_path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    if (isTRUE(do_plot)) {
      mdcat(paste0("- [estimate plotting table: ", basename(export_path), "](", export_path, ")"))
    }
  }

  # Plot (optional)
  if (isTRUE(do_plot)) {
    pal <- c("#d73027", "#8844aa", "#4575b4", "grey70")
    names(pal) <- c(
      "Significant in both",
      sprintf("Significant in %s only", contrasts[1]),
      sprintf("Significant in %s only", contrasts[2]),
      "Not significant"
    )

    gg <- ggplot2::ggplot(
      joined,
      ggplot2::aes(
        x     = value.A,
        y     = y_plot,
        color = status,
        shape = shape
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(size = 2, alpha = 0.7) +
      ggplot2::scale_color_manual(values = pal, name = NULL) +
      ggplot2::scale_shape_identity(guide = "none") +
      ggplot2::scale_y_continuous(limits = y.lim) +
      ggplot2::labs(
        x     = paste0("estimate (", contrasts[1], ")"),
        y     = paste0("estimate (", contrasts[2], ")"),
        title = paste0("estimate: ", contrasts[1], " vs ", contrasts[2])
      ) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5),
        # Position the legend at the bottom to make more horizontal room for the
        # plot space
        legend.position = "bottom"
      )

    # Labelling logic
    if (isTRUE(label_concord)) {
      # cross-quadrant + significant in both
      to_label <- joined %>%
        dplyr::filter(
          status == "Significant in both",
          (value.A > 0 & value.B < 0) | (value.A < 0 & value.B > 0)
        ) %>%
        dplyr::mutate(rank_abs = abs(value.A) + abs(value.B)) %>%
        dplyr::slice_max(rank_abs, n = label_n, with_ties = FALSE)

      gg <- gg +
        ggrepel::geom_label_repel(
          data              = to_label,
          mapping           = ggplot2::aes(label = Assay),
          inherit.aes       = TRUE,
          color             = "black",
          fill              = "white",
          segment.color     = "black",
          label.size        = 0.25,
          box.padding       = 0.3,
          point.padding     = 0.25,
          max.overlaps      = nrow(joined),
          show.legend       = FALSE
        )

    } else if (isTRUE(label_sig_both)) {
      # significant in both, any sign
      to_label <- joined %>%
        dplyr::filter(status == "Significant in both") %>%
        dplyr::mutate(rank_abs = abs(value.A) + abs(value.B)) %>%
        dplyr::slice_max(rank_abs, n = label_n, with_ties = FALSE)

      gg <- gg +
        ggrepel::geom_label_repel(
          data              = to_label,
          mapping           = ggplot2::aes(label = Assay),
          inherit.aes       = TRUE,
          color             = "black",
          fill              = "white",
          segment.color     = "black",
          label.size        = 0.25,
          box.padding       = 0.7,
          point.padding     = 0.25,
          max.overlaps      = nrow(joined),
          show.legend       = FALSE
        )
    }

    # Save static ggplot to PDF
    out_file <- file.path(
      out_dir,
      paste0(paste(contrasts, collapse = "_"), "_estimate_plot.pdf")
    )
    ggplot2::ggsave(filename = out_file, plot = gg)

    mdcat(paste0(
      "- [estimate Plot: ", contrasts[1], " vs ", contrasts[2], "](", out_file, ")"
    ))

    # Optional Plotly conversion with onRender callback to open Protein Atlas URL on click
    if (isTRUE(plotly)) {
      ply <- plotly::ggplotly(
        gg +
          ggplot2::aes(
            text = paste0(
              "Protein: ", Assay,
              "<br>Panel: ", Panel,
              "<br>Status: ", status,
              "<br>estimate(", contrasts[1], "): ", round(value.A, 3),
              "<br>estimate(", contrasts[2], "): ", round(value.B, 3),
              "<br>Adj. P(", contrasts[1], "): ", formatC(adjp.A, format = "e", digits = 2),
              "<br>Adj. P(", contrasts[2], "): ", formatC(adjp.B, format = "e", digits = 2)
            ),
            customdata = paste(GeneID, Assay, sep = "|")
          ),
        tooltip = "text"
      ) %>%
        plotly::layout(
          legend = list(
            orientation = "h",
            x           = 0.5,
            xanchor     = "center",
            y           = -0.15,
            yanchor     = "top"
          )
        ) %>%
        htmlwidgets::onRender(
          "function(el, x) {
            el.on('plotly_click', function(data) {
              var pt    = data.points[0];
              var parts = pt.customdata.split('|');
              var gene  = parts[0];
              var assay = parts[1];
              var url   = 'https://www.proteinatlas.org/' + gene + '-' + assay + '/tissue';
              window.open(url, '_blank');
            });
          }"
        )

      subchunkify(ply, fig.height = 10, fig.width = 10)
    } else {
      print(gg)
    }
  } # end do_plot

  # Return table if requested
  if (isTRUE(return_table)) {
    return(plot_tbl)
  } else {
    return(invisible(NULL))
  }
}


#' Interaction plot for paired samples
#'
#' @description
#' Plots paired measurements linked by patient, colored by sex and shaped
#' by pre/post treatment. Labels patient IDs. Saves PDF.
#'
#' @param se_list Named list of SummarizedExperiment objects.
#' @param se Name of element in `se_list` to plot.
#' @param y Column name in colData for y-axis.
#' @param x Column name in colData for x-axis.
#'
#' @return Invisibly NULL. Side effect: plot and PDF.
#' @export
interaction_plot <- function(se_list,
                             se,
                             y,
                             x,
                             out_dir) {

  # Fetch metadata
  colData <- tibble::as_tibble(SummarizedExperiment::colData(se_list[[se]]))

  # Checks
  stopifnot(
    is.list(se_list),
    se %in% names(se_list),
    y  %in% colnames(colData),
    x  %in% colnames(colData),
    "PatientID" %in% colnames(colData),
    "Sex" %in% colnames(colData),
    "PrePost" %in% colnames(colData)
  )

  # Mid-point for each patient where we place the label
  label_colData <- aggregate(
    cbind(x_val = colData[[x]], y_val = colData[[y]]),
    by   = list(PatientID = colData$PatientID),
    FUN  = mean,
    na.rm = TRUE
  )

  # Plot
  gg <- ggplot2::ggplot(
    colData,
    ggplot2::aes(
      x      = .data[[x]],
      y      = .data[[y]],
      group  = PatientID,
      color  = Sex
    )
  ) +
    # Draw lines colored by Sex
    ggplot2::geom_line() +
    # Draw points colored by Sex and shaped by PrePost
    ggplot2::geom_point(
      ggplot2::aes(
        shape = PrePost
      ),
      size = 3
    ) +
    # Add patient labels at midpoint
    ggrepel::geom_label_repel(
      data = label_colData,
      ggplot2::aes(
        x     = x_val,
        y     = y_val,
        label = PatientID
      ),
        inherit.aes        = FALSE,
        color              = "black",
        fill               = "white",
        segment.color      = "black",
        min.segment.length = 0,
        label.size         = 0.25,
        box.padding        = 0.7,
        point.padding      = 0.25,
        max.overlaps       = nrow(colData),
        show.legend        = FALSE
    ) +
    # Manual scales for color and shape
    ggplot2::scale_color_manual(
      values = c(
        "M" = "blue",
        "F" = "red"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "pre"  = 16,  # circle
        "post" = 17   # triangle
      )
    ) +
    ggplot2::labs(
      title = sprintf("%s vs %s for paired samples", y, x),
      x     = x,
      y     = y
    )

  # Save
  path <- file.path(out_dir, sprintf("interaction_%s_vs_%s_%s.pdf", y, x, se))
  ggplot2::ggsave(filename = path, plot = gg)
  mdcat("- [", sprintf("Interaction %s vs %s PDF", y, x), "](", path, ")")
  # Render in report
  print(gg)

  invisible(NULL)
}
