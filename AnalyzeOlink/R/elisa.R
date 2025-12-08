#' Plot ELISA Outliers
#'
#' @description
#' Creates a multi‐page PDF (one page per assay) showing boxplots of ELISA
#' concentrations, highlighting outliers detected using Tukey’s fences.
#' Outliers are labeled with SampleID. Fences can be computed per assay or
#' within additional grouping variables.
#'
#' @param df Data frame with columns: Assay, SampleID, and ELISA (log10).
#' @param outdir From config specifying output directory.
#' @param value_col Name of numeric column to analyze (default `"ELISA"`).
#' @param group_col Column name for x-axis grouping in plots (default `"Treatment"`).
#' @param by Columns to compute fences within (default `c("Assay")`).
#' @param k Tukey fence multiplier (default 1.5).
#' @param outfile Output PDF file name (default `"elisa_outliers.pdf"`).
#' @param width,height PDF size in inches.
#' @param detect_outliers Logical; if TRUE (default), compute and label outliers.
#'
#' @return Invisibly returns the PDF path.
#' @export
plot_elisa_outliers <- function(df,
                                outdir,
                                value_col = "ELISA",
                                group_col = "Group",
                                by = c("Assay"),
                                k = 1.5,
                                outfile = "elisa_outliers.pdf",
                                width = 7, height = 5,
                                detect_outliers = TRUE) {

  # Validate inputs
  stopifnot(is.character(value_col), is.character(group_col), is.character(by))
  needed_cols <- unique(c("Assay", "SampleID", value_col, group_col, by))
  if (!all(needed_cols %in% names(df))) {
    stop("Input df is missing required columns: ",
         paste(setdiff(needed_cols, names(df)), collapse = ", "))
  }
  pdf_path <- file.path(outdir, outfile)

  # Symbols
  val    <- rlang::sym(value_col)
  gr_col <- rlang::sym(group_col)

  # Prepare stats (Tukey fences) within `by` if detecting outliers
  if (isTRUE(detect_outliers)) {
    stats <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
      dplyr::summarise(
        n_non_na = sum(!is.na(!!val)),
        q1  = stats::quantile(!!val, 0.25, na.rm = TRUE, names = FALSE),
        q3  = stats::quantile(!!val, 0.75, na.rm = TRUE, names = FALSE),
        iqr = q3 - q1,
        lower = q1 - k * iqr,
        upper = q3 + k * iqr,
        .groups = "drop"
      )

    # Join stats & flag outliers
    flagged <- df %>%
      dplyr::left_join(stats, by = by) %>%
      dplyr::mutate(
        .is_outlier = dplyr::case_when(
          is.na(!!val) ~ FALSE,
          !!val < .data$lower ~ TRUE,
          !!val > .data$upper ~ TRUE,
          TRUE ~ FALSE
        )
      )
  } else {
    # No detection: everyone "looks like" an outlier; no lines needed
    stats <- NULL
    flagged <- df %>%
      dplyr::mutate(
        .is_outlier = TRUE,
        lower = NA_real_,
        upper = NA_real_
      )
  }

  # Open multi-page PDF and plot one page per Assay level (within the chosen `by`)
  assays <- unique(flagged$Assay)
  grDevices::pdf(pdf_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (aa in assays) {
    page_df <- flagged %>% dplyr::filter(.data$Assay == aa)

    p <- ggplot2::ggplot(page_df, ggplot2::aes(x = !!gr_col, y = !!val)) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
      # Style points: when detecting, outliers get alpha 1, others faded.
      # When not detecting, .is_outlier == TRUE for all -> alpha == 1 for all.
      ggplot2::geom_jitter(
        ggplot2::aes(alpha = .is_outlier),
        width = 0.15, height = 0, size = 2, show.legend = FALSE
      ) +
      ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1))

    # Only draw Tukey fences & labels if we are detecting outliers
    if (isTRUE(detect_outliers)) {
      page_stats <- stats %>% dplyr::filter(.data$Assay == aa)
      p <- p +
        ggplot2::geom_hline(data = page_stats, ggplot2::aes(yintercept = lower),
                            linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
        ggplot2::geom_hline(data = page_stats, ggplot2::aes(yintercept = upper),
                            linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
        ggrepel::geom_text_repel(
          data = subset(page_df, .is_outlier),
          ggplot2::aes(label = .data$SampleID),
          size = 3, max.overlaps = Inf, min.segment.length = 0
        )
    }

    p <- p +
      ggplot2::labs(
        title = paste0(
          "Assay: ", aa,
          "  (Tukey k=", k, "; by=", paste(by, collapse = ","), ")"
        ),
        x = group_col,
        y = "ELISA (log10 concentration)"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
      )

    print(p)
  }

  invisible(pdf_path)
}


#' Identify ELISA Outliers by Assay
#'
#' @description
#' Uses Tukey’s fences (k × IQR) to identify outliers in ELISA data grouped by
#' assay (and optionally additional variables). Returns a list of outlier rows
#' per assay.
#'
#' @param df Data frame with columns: Assay, SampleID, and ELISA (log10).
#' @param value_col Name of numeric column to analyze (default `"ELISA"`).
#' @param by Columns to compute fences within (default `c("Assay")`).
#' @param k Tukey fence multiplier (default 1.5).
#'
#' @return Named list of tibbles, one per assay, containing only outlier rows.
#' @export
tukey_outliers_list <- function(df,
                                value_col = "ELISA",
                                by = c("Assay"),
                                k = 1.5) {
  pkgs <- c("dplyr", "rlang")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Please install required packages: ", paste(missing, collapse = ", "))

  stopifnot(is.character(value_col), is.character(by))
  needed_cols <- unique(c("Assay", "SampleID", value_col, by))
  if (!all(needed_cols %in% names(df))) {
    stop("Input df is missing required columns: ",
         paste(setdiff(needed_cols, names(df)), collapse = ", "))
  }

  val <- rlang::sym(value_col)

  stats <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::summarise(
      n_non_na = sum(!is.na(!!val)),
      q1  = stats::quantile(!!val, 0.25, na.rm = TRUE, names = FALSE),
      q3  = stats::quantile(!!val, 0.75, na.rm = TRUE, names = FALSE),
      iqr = q3 - q1,
      lower = q1 - k * iqr,
      upper = q3 + k * iqr,
      .groups = "drop"
    )

  flagged <- df %>%
    dplyr::left_join(stats, by = by) %>%
    dplyr::mutate(
      .is_outlier = dplyr::case_when(
        is.na(!!val) ~ FALSE,
        !!val < .data$lower ~ TRUE,
        !!val > .data$upper ~ TRUE,
        TRUE ~ FALSE
      )
    )

  # Split into a named list by Assay (keep empty tibbles for completeness)
  assays <- sort(unique(flagged$Assay))
  out_list <- lapply(assays, function(a) {
    flagged %>%
      dplyr::filter(.data$Assay == a, .is_outlier) %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(by)))
  })
  names(out_list) <- assays
  out_list
}


#' Filter ELISA Outliers
#'
#' @description
#' Removes or keeps only outlier rows from ELISA data, based on a list of
#' outliers computed with `tukey_outliers_list()`.
#'
#' @param df Data frame with ELISA measurements.
#' @param outliers_by_assay List of outlier tibbles from `tukey_outliers_list()`.
#' @param keys Columns used to match outliers (default `c("Assay","SampleID")`).
#' @param keep Logical; if FALSE (default), remove outliers. If TRUE, keep only outliers.
#'
#' @return Tibble filtered according to `keep`.
#' @export
filter_elisa_outliers <- function(df,
                                  outliers_by_assay,
                                  keys = c("Assay", "SampleID"),
                                  keep = FALSE) {
  # dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package.")
  }

  # validate inputs
  if (!all(keys %in% names(df))) {
    stop("`df` is missing required key columns: ",
         paste(setdiff(keys, names(df)), collapse = ", "))
  }
  if (!is.list(outliers_by_assay)) {
    stop("`outliers_by_assay` must be a named list of tibbles as returned by tukey_outliers_list().")
  }

  # bind all outlier keys across assays (handles empty elements gracefully)
  # keep only the key columns and drop NAs
  if (length(outliers_by_assay) == 0) {
    out_keys <- dplyr::tibble(!!!setNames(vector("list", length(keys)), keys))
    out_keys <- out_keys[0, ] # empty tibble with proper cols
  } else {
    out_keys <- dplyr::bind_rows(outliers_by_assay) %>%
      dplyr::select(dplyr::all_of(keys)) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(keys), ~ !is.na(.x))) %>%
      dplyr::distinct()
  }

  # if no outlier keys, return according to `keep`
  if (nrow(out_keys) == 0) {
    return(if (keep) df[0, , drop = FALSE] else df)
  }

  # filter by join
  if (keep) {
    dplyr::semi_join(df, out_keys, by = keys)
  } else {
    dplyr::anti_join(df, out_keys, by = keys)
  }
}


#' Build Outlier‐Removed SummarizedExperiment
#'
#' @description
#' Creates a copy of a SummarizedExperiment with outliers removed from its
#' `metadata$long_data`. Optionally restricts to samples with non‐outlier rows
#' and removes QC warning rows.
#'
#' @param se A SummarizedExperiment to copy.
#' @param long_df Long‐format ELISA data frame with Assay, SampleID, and ELISA.
#' @param by Columns to compute Tukey fences within (default `c("Assay")`).
#' @param k Tukey fence multiplier (default 1.5).
#' @param value_col Name of numeric column to analyze (default `"ELISA"`).
#' @param outliers_by_assay Optional precomputed list from `tukey_outliers_list()`.
#' @param restrict_samples Logical; if TRUE, drop SE columns for samples with no
#'   remaining non‐outlier rows.
#' @param rm_qc_warn_md Logical; if TRUE, remove rows with QC_Warning/Assay_Warning == "WARN".
#' @param sample_id_coldata Column in colData(se) containing sample IDs (default `"SampleID"`).
#'
#' @return A SummarizedExperiment with updated `metadata$long_data`.
#' @export
build_outlier_removed_se <- function(se,
                                     long_df,
                                     by = c("Assay"),
                                     k = 1.5,
                                     value_col = "ELISA",
                                     outliers_by_assay = NULL,
                                     restrict_samples = FALSE,
                                     rm_qc_warn_md = FALSE,
                                     sample_id_coldata = "SampleID") {
  # deps
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("Please install 'SummarizedExperiment'.")
  if (!requireNamespace("S4Vectors", quietly = TRUE))
    stop("Please install 'S4Vectors'.")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Please install 'dplyr'.")

  # validate inputs
  need_cols <- c("Assay", "SampleID", value_col)
  miss <- setdiff(need_cols, names(long_df))
  if (length(miss)) stop("`long_df` is missing required columns: ", paste(miss, collapse = ", "))

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (!sample_id_coldata %in% names(cd)) {
    stop("`colData(se)` is missing column '", sample_id_coldata, "'.")
  }

  # compute outliers if not supplied
  if (is.null(outliers_by_assay)) {
    outliers_by_assay <- tukey_outliers_list(
      df = long_df,
      value_col = value_col,
      by = by,
      k = k
    )
  }

  # remove outliers from long_df
  cleaned_long <- filter_elisa_outliers(
    df = long_df,
    outliers_by_assay = outliers_by_assay,
    keys = c("Assay", "SampleID"),
    keep = FALSE
  )

  # optionally drop QC warnings from metadata if columns exist
  if (isTRUE(rm_qc_warn_md)) {
    if ("QC_Warning" %in% names(cleaned_long)) {
      cleaned_long <- dplyr::filter(cleaned_long, .data$QC_Warning != "WARN")
    }
    if ("Assay_Warning" %in% names(cleaned_long)) {
      cleaned_long <- dplyr::filter(cleaned_long, .data$Assay_Warning != "WARN")
    }
  }

  # optionally restrict SE columns to samples present in the cleaned metadata
  se_subset <- se
  if (isTRUE(restrict_samples)) {
    keep_samples <- unique(cleaned_long$SampleID)
    if (length(keep_samples) == 0) stop("No samples remain after outlier filtering.")
    keep_idx <- cd[[sample_id_coldata]] %in% keep_samples
    se_subset <- se[, keep_idx, drop = FALSE]
  }

  # replace metadata$long_data and validate
  S4Vectors::metadata(se_subset)$long_data <- cleaned_long
  check_se_object(se_subset)

  se_subset
}
