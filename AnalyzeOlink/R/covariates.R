#' Prepare covariate colData for each contrast
#'
#' @description
#' Builds a named list of per-contrast metadata frames used in covariate analyses.
#' Each element contains:
#'   * `df`  — sample-level metadata for that contrast
#'   * `col` — grouping column to compare (default "Group")
#'   * `g1`  — first group label
#'   * `g2`  — second group label
#'
#' This helper is intentionally explicit and user-configurable. It does **not**
#' auto-detect covariates or groups. Test data expects the universal columns
#' `Group`, `Cov1`, and `Cov2`.
#'
#' @param se_list Named list of SummarizedExperiment objects. The SE named in
#'   each `contrast_specs[[i]]$se_object_name` must exist in `se_list`.
#' @param contrast_specs Named list describing each contrast. Each element must
#'   include:
#'   * `se_object_name`  — name of SE in `se_list` (e.g., "Disease")
#'   * `sampletable_col` — grouping column in colData (default "Group")
#'   * `g1`              — first group label
#'   * `g2`              — second group label
#'   Optionally:
#'   * `columns`         — character vector of colData columns to keep
#'                          (default keeps all).
#'
#' @return Named list. Each element is a list with components `df`, `col`, `g1`, `g2`.
#'
#' @examples
#' contrast_specs <- list(
#'   disease_vs_control = list(
#'     se_object_name  = "Disease",
#'     sampletable_col = "Group",
#'     g1 = "disease",
#'     g2 = "control",
#'     columns = c("SampleID", "Group", "Cov1", "Cov2")
#'   )
#' )
#' coldata_list <- prepare_covariate_coldata(se_list, contrast_specs)
#'
#' @export
prepare_covariate_coldata <- function(se_list, contrast_specs) {
  stopifnot(is.list(se_list), is.list(contrast_specs), length(contrast_specs) > 0)

  out <- purrr::imap(contrast_specs, function(spec, contrast_name) {
    # required fields
    se_object_name  <- spec$se_object_name
    sampletable_col <- spec$sampletable_col %||% "Group"
    g1              <- spec$g1
    g2              <- spec$g2
    columns         <- spec$columns %||% NULL

    if (is.null(se_object_name) || !se_object_name %in% names(se_list)) {
      stop("Contrast '", contrast_name, "' refers to missing SE '", se_object_name, "'.")
    }
    if (is.null(g1) || is.null(g2)) {
      stop("Contrast '", contrast_name, "' must define g1 and g2.")
    }

    cd <- as.data.frame(SummarizedExperiment::colData(se_list[[se_object_name]]))

    if (!sampletable_col %in% names(cd)) {
      stop("Grouping column '", sampletable_col, "' not found in colData for contrast '", contrast_name, "'.")
    }

    if (!is.null(columns)) {
      missing_cols <- setdiff(columns, names(cd))
      if (length(missing_cols) > 0) {
        stop("Requested columns missing in contrast '", contrast_name, "': ",
             paste(missing_cols, collapse = ", "))
      }
      cd <- dplyr::select(cd, dplyr::all_of(columns))
    }

    # make sure grouping col is a factor w/ expected levels present
    cd[[sampletable_col]] <- factor(cd[[sampletable_col]])
    levs <- levels(cd[[sampletable_col]])
    if (!all(c(g1, g2) %in% levs)) {
      stop(
        "Contrast '", contrast_name, "': grouping levels do not include g1/g2. ",
        "Found: {", paste(levs, collapse = ", "), "}; expected g1='", g1, "', g2='", g2, "'."
      )
    }

    list(
      df  = cd,
      col = sampletable_col,
      g1  = g1,
      g2  = g2
    )
  })

  out
}


#' Format Covariate Test Results
#'
#' @description
#' Formats the results of a covariate test:
#' * Adds a `Threshold` column indicating significance at `alpha` on adjusted p-values.
#' * Optionally formats p-values and statistics in scientific notation.
#' * Rounds other numeric columns to 3 decimal places.
#' * Removes row names and returns a plain data.frame.
#'
#' @param df A data.frame or tibble with at least a `p.adj` column.
#' @param alpha Numeric significance threshold for adjusted p-values.
#' @param scientific_notation Logical; if TRUE, format p-values/statistics in scientific notation.
#'
#' @return A formatted data.frame with a `Threshold` column.
#' @export
format_cov_results <- function(df, alpha, scientific_notation = FALSE) {
  df <- df %>%
    dplyr::mutate(
      Threshold = ifelse(.data$p.adj < alpha, "Significant", "Non-significant")
    )

  if (isTRUE(scientific_notation)) {
    df <- df %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::where(is.numeric),
          ~ if (dplyr::cur_column() %in% c("p.value", "statistic")) {
            sprintf("%.2e", .)
          } else {
            round(., 3)
          }
        )
      )
  }

  df <- as.data.frame(df) %>%
    tibble::remove_rownames()

  return(df)
}


#' Plot Covariate Distribution
#'
#' @description
#' Visualizes the distribution of a covariate across two groups. If the covariate
#' column is numeric, a histogram is drawn. If categorical, a bar chart is drawn.
#'
#' @param tbl Data frame containing grouping and covariate columns.
#' @param group_col Character string, name of the grouping column.
#' @param groups Character vector of length 2, the group labels to compare.
#' @param covariate_col Character string, name of the covariate column.
#'
#' @return A ggplot object showing the distribution of the covariate.
#' @export
plot_covariate_dist <- function(tbl, group_col, groups, covariate_col) {
  tbl <- tbl %>%
    dplyr::filter(.data[[group_col]] %in% groups)

  p <- ggplot2::ggplot(tbl, ggplot2::aes(fill = .data[[group_col]]))

  if (is.numeric(tbl[[covariate_col]])) {
    p <- p +
      ggplot2::aes(x = .data[[covariate_col]]) +
      ggplot2::geom_histogram(color = "black", bins = 30) +
      ggplot2::labs(
        x = covariate_col, y = "Count",
        title = paste("Histogram of", covariate_col, "for", paste(groups, collapse = ", ")),
        fill = group_col
      )
  } else {
    p <- p +
      ggplot2::aes(x = .data[[covariate_col]]) +
      ggplot2::geom_bar(position = ggplot2::position_dodge(), color = "black") +
      ggplot2::labs(
        x = covariate_col, y = "Count",
        title = paste("Histogram of", covariate_col, "for", paste(groups, collapse = ", ")),
        fill = group_col
      )
  }

  return(p)
}


#' Print Covariate Analysis Results
#'
#' @description
#' For a single covariate, prints a tabset of contrasts with:
#' * Download link for Excel results
#' * Download link for per-contrast PDF plots
#' * Embedded ggplot and interactive table of test results
#'
#' @param coldata_list List of contrast metadata from `prepare_covariate_coldata()`.
#' @param cov_list Nested list of covariate results from `test_for_covariates()`.
#' @param config Configuration list containing at least `alpha`.
#' @param xlsx_path Path to Excel workbook for download link.
#' @param plots_dir Directory where PDF plots are saved.
#' @param covariate Character string; covariate name.
#' @param plot_fn Function used to generate plots (signature `(df, group_col, groups, covariate)`).
#'
#' @return Invisibly NULL. Called for side-effects.
#' @export
print_covariate_analysis <- function(
  coldata_list,
  cov_list,
  config,
  xlsx_path,
  plots_dir,
  covariate,
  plot_fn
) {
  cov_suffix <- tolower(covariate)
  mdcat(paste0("- [Download ", covariate, " results as Excel](", xlsx_path, ")"))

  for (nm in names(coldata_list)) {
    cd <- coldata_list[[nm]]
    mdcat(paste0("### ", nm))

    plt <- plot_fn(cd$df, cd$col, c(cd$g1, cd$g2), covariate)
    pdf_file <- file.path(plots_dir, paste0(cov_suffix, "_", nm, ".pdf"))
    ggplot2::ggsave(pdf_file, plot = plt, device = "pdf", width = 7, height = 5)

    mdcat(paste0("- [Download ", covariate, " plot](", pdf_file, ")"))
    subchunkify(plt)

    tbl <- cov_list[[covariate]][[nm]] %>%
      format_cov_results(config$alpha, scientific_notation = TRUE)
    subchunkify(DT::datatable(tbl))
  }

  invisible(NULL)
}


#' Save Covariate Results to Excel
#'
#' @description
#' Writes formatted covariate results to an Excel workbook, with one sheet per
#' covariate × contrast combination. Sheet names are truncated to 31 characters
#' if necessary to meet Excel's limitation.
#'
#' @param cov_list Nested list of covariate results (from `test_for_covariates()`).
#' @param covariate_tests Named character vector or list mapping covariate names
#'   to test types (e.g. `list(Age = "wilcox", Sex = "chisq")`).
#' @param out_dir Directory where workbook will be written.
#' @param file_name Name of the Excel file (default `"covariates.xlsx"`).
#'
#' @return Invisibly returns the workbook path.
#' @export
save_covariates <- function(cov_list, covariate_tests, out_dir, file_name = "covariates.xlsx") {
  missing_covs <- setdiff(names(covariate_tests), names(cov_list))
  if (length(missing_covs) > 0) {
    stop("The following covariates are missing in cov_list: ", paste(missing_covs, collapse = ", "))
  }

  wb <- openxlsx::createWorkbook()

  for (cov in names(covariate_tests)) {
    contrast_list <- cov_list[[cov]]
    for (contrast in names(contrast_list)) {
      sheet_name <- paste(cov, contrast, sep = "_")
      if (nchar(sheet_name) > 31) {
        sheet_name <- substr(sheet_name, 1, 31)
      }
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(
        wb,
        sheet = sheet_name,
        x     = contrast_list[[contrast]],
        withFilter = TRUE
      )
    }
  }

  out_path <- file.path(out_dir, file_name)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)

  invisible(out_path)
}


#' Run a Covariate Test
#'
#' @description
#' Runs either a Wilcoxon rank-sum test (for numeric covariates) or a chi-square
#' test (for categorical covariates) between two groups.
#'
#' @param tbl Data frame of metadata.
#' @param test Character; `"wilcox"` or `"chisq"`.
#' @param col Grouping column name.
#' @param var Covariate column name.
#' @param g1 First group label.
#' @param g2 Second group label.
#'
#' @return A tibble of test results from `broom::tidy()`.
#' @export
run_covariate_test <- function(tbl, test, col, var, g1, g2) {
  if (test == "wilcox") {
    test_res <- stats::wilcox.test(
      tbl[[var]][tbl[[col]] == g1],
      tbl[[var]][tbl[[col]] == g2]
    )
    broom::tidy(test_res)
  } else if (test == "chisq") {
    tab <- table(tbl[[col]], tbl[[var]])
    broom::tidy(stats::chisq.test(tab))
  } else {
    stop("Unsupported test: must be 'wilcox' or 'chisq'")
  }
}


#' Run Covariate Tests Across Contrasts
#'
#' @description
#' Applies covariate tests for each specified covariate across all contrasts
#' defined in `coldata_list`. Uses `run_covariate_test()` internally and then
#' adds adjusted p-values with `cov_adjust()`.
#'
#' @param coldata_list List of contrasts from `prepare_covariate_coldata()`.
#' @param covariate_tests Named list mapping covariates to test types.
#' @param alpha Numeric significance threshold.
#'
#' @return A nested list: `cov_list[[covariate]][[contrast]] → tibble of results`.
#' @export
test_for_covariates <- function(coldata_list, covariate_tests, alpha) {
  cov_list <- lapply(names(covariate_tests), function(covariate) {
    test_type <- covariate_tests[[covariate]]
    res_per_contrast <- lapply(coldata_list, function(cd) {
      run_covariate_test(
        tbl   = cd$df,
        test  = test_type,
        col   = cd$col,
        var   = covariate,
        g1    = cd$g1,
        g2    = cd$g2
      )
    })
    names(res_per_contrast) <- names(coldata_list)
    res_per_contrast
  })

  names(cov_list) <- names(covariate_tests)

  cov_list <- cov_adjust(cov_list)
  cov_list
}


#' Adjust P-values in Covariate Results
#'
#' @description
#' Applies multiple testing correction (default `"BH"`) to the p-values in a
#' nested covariate result list. Adds a `p.adj` column to each tibble.
#'
#' @param cov_list Nested list of covariate results (from `test_for_covariates()`).
#' @param method Character string for adjustment method (default `"BH"`).
#'
#' @return Nested list with the same structure as `cov_list`, but with added `p.adj`.
#' @export
cov_adjust <- function(cov_list, method = "BH") {
  lapply(cov_list, function(tables) {
    raw_p <- vapply(tables, `[[`, numeric(1), "p.value")
    adj_p <- stats::p.adjust(raw_p, method = method)
    mapply(function(tbl, p_adj) {
      tbl$p.adj <- p_adj
      tbl
    }, tables, adj_p, SIMPLIFY = FALSE)
  })
}
