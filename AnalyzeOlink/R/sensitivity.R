#' Sensitivity scatter-plot
#'
#' @description
#' Faceted scatter-plot comparing effect sizes (partial eta² or conditional R²)
#' with vs. without a given sample or model term. Colors show significance
#' categories:
#' \itemize{
#'   \item Significant in both
#'   \item Significant with element only
#'   \item Significant without element only
#'   \item Not significant
#' }
#'
#' @param res_list Named list of results from \code{run_tests()}.
#' @param out_dir from config like \code{output_paths$sensitivity_pea$plots}.
#' @param contrasts Character vector of length 2 (with vs without).
#' @param effects Character vector of model terms to plot.
#' @param element Label of the sample/term removed (used in legend & labels).
#' @param type Either `"anova"` (generalized eta²) or `"lmer"` (conditional R²).
#' @param reference Logical. If TRUE, reduce to binary significant vs not.
#' @param label Optional vector of Assay names to annotate.
#' @param column Override statistic column; defaults to `"pes"`.
#' @param plotly Logical. If TRUE, also render an interactive Plotly widget.
#' @param color Column with significance flag (default `"Threshold"`).
#'
#' @return Invisibly NULL. Side effect: saves PDF and prints plot (ggplot or plotly).
#'
#' @examples
#' \dontrun{
#' sensitivity(
#'   res_list,
#'   config     = list(output_paths = list(sensitivity = list(plots = "plots/"))),
#'   contrasts  = c("disease_vs_control_no_A2", "disease_vs_control"),
#'   effects    = c("Treatment","Age_c","Sex"),
#'   element    = "C108",
#'   type       = "anova",
#'   plotly     = TRUE
#' )
#' }
#' @export
sensitivity <- function(res_list,
                        out_dir,
                        contrasts,
                        effects,
                        element,
                        type      = c("anova", "lmer"),
                        reference = FALSE,
                        label     = NULL,
                        column    = NULL,
                        plotly    = FALSE,
                        color     = "Threshold") {

  type <- match.arg(type)
  # Determine statistic column based on type unless overridden
  if (is.null(column)) {
    column <- switch(type,
                     anova = "pes",
                     lmer  = "pes")
  }

  stopifnot(
    length(contrasts) == 2,
    all(contrasts %in% names(res_list)),
    is.list(res_list[[contrasts[1]]]),
    is.list(res_list[[contrasts[2]]])
  )

  pull_contrast <- function(ctr) {
    an_tbl <- res_list[[ctr]]$anova %>%
      dplyr::filter(Effect %in% effects)

    if (!column %in% names(an_tbl)) {
      stop("Column '", column, "' not found in contrast ", ctr)
    }

    an_tbl %>%
      dplyr::select(Assay, Panel, Effect,
                    !!rlang::sym(column), !!rlang::sym(color)) %>%
      dplyr::mutate(Assay_Panel = paste(Assay, Panel, sep = "_")) %>%
      dplyr::rename(
        value = !!rlang::sym(column),
        sig   = !!rlang::sym(color)
      ) %>%
      dplyr::mutate(
        sig      = (sig == "Significant"),
        Contrast = ctr
      )
  }

  df_with    <- pull_contrast(contrasts[1])
  df_without <- pull_contrast(contrasts[2])

  joined <- dplyr::full_join(
    df_with, df_without,
    by     = c("Assay", "Panel", "Effect"),
    suffix = c(".with", ".without")
  ) %>%
    dplyr::mutate( Assay_Panel = paste(Assay, Panel, sep = "_") )

  axis_label <- switch(type,
                       anova = "partial eta²",
                       lmer  = "partial eta²")

  if (isTRUE(reference)) {
    joined <- joined %>%
      dplyr::mutate(
        status = dplyr::if_else(
          sig.with | sig.without,
          "Significant in both",
          "Not significant"
        ),
        status = factor(
          status,
          levels = c("Not significant", "Significant in both")
        )
      ) %>%
      dplyr::arrange(status)

    pal <- c("grey70", "#8844aa")
    names(pal) <- levels(joined$status)

    x_lab      <- axis_label
    y_lab      <- axis_label
    plot_title <- "Reference comparison: effect sizes"

  } else {
    joined <- joined %>%
      dplyr::mutate(
        status = dplyr::case_when(
          sig.with  & sig.without ~ "Significant in both",
          sig.with  & !sig.without ~ sprintf("Significant with %s only",    element),
          !sig.with  & sig.without ~ sprintf("Significant without %s only", element),
          TRUE                     ~ "Not significant"
        ),
        status = factor(
          status,
          levels = c(
            "Not significant",
            sprintf("Significant with %s only",    element),
            sprintf("Significant without %s only", element),
            "Significant in both"
          )
        )
      ) %>%
      dplyr::arrange(status)

    pal <- c(
      "grey70",
      "#d73027",
      "#4575b4",
      "#8844aa"
    )
    names(pal) <- levels(joined$status)

    x_lab      <- sprintf("%s (with %s)",    axis_label, element)
    y_lab      <- sprintf("%s (without %s)", axis_label, element)
    plot_title <- sprintf("Sensitivity analysis: effect sizes with vs. without %s", element)
  }

  gg <- ggplot2::ggplot(
    joined,
    ggplot2::aes(
      x     = value.with,
      y     = value.without,
      color = status
    )
  ) +
    ggplot2::geom_abline(linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 1, alpha = 0.6) +
    ggplot2::scale_color_manual(values = pal, name = NULL) +
    ggplot2::facet_wrap(~ Effect, scales = "free") +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = plot_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text      = ggplot2::element_text(face = "bold"),
      plot.title      = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    )

  if (!is.null(label) && length(label) > 0) {
    gg <- gg +
      ggrepel::geom_label_repel(
        data = joined %>% dplyr::filter(Assay_Panel %in% label),
        ggplot2::aes(label = Assay),
        inherit.aes        = TRUE,
        color             = "black",
        fill               = "white",
        segment.color     = "black",
        min.segment.length = 0,
        label.size         = 0.05,
        box.padding        = 1,
        point.padding      = 0.25,
        max.overlaps       = nrow(joined),
        show.legend        = FALSE
      )
  }

  pdf_path <- file.path(
    out_dir,
    paste0(paste(contrasts, collapse = "_V_"), "_sensitivity_plot.pdf")
  )
  ggplot2::ggsave(pdf_path, gg)
  mdcat(paste0("- [Sensitivity PDF](", pdf_path, ")"))

  if (isTRUE(plotly)) {
    ply <- plotly::ggplotly(
             gg +
               ggplot2::aes(
                 text = paste(
                   "Protein:", Assay,
                   "<br>Panel:", Panel,
                   "<br>Status:", status
                 )
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
           )
    subchunkify(ply)
  } else {
    print(gg)
  }

  invisible(NULL)
}
