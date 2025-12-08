# Global ggplot theme: Times New Roman & large text everywhere
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 14, base_family = "Times New Roman") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Times New Roman", color = "black"),
      axis.title = ggplot2::element_text(size = 16, color = "black"),
      axis.text  = ggplot2::element_text(size = 14, color = "black"),
      legend.title = ggplot2::element_text(size = 16, color = "black"),
      legend.text  = ggplot2::element_text(size = 14, color = "black"),
      strip.text = ggplot2::element_text(size = 14, color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.background = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
)

# Load the extrafonts package to enable Times New Roman
library(extrafont)
# If this is the first time running on your system, you may need to run font_import()
# Here is a work around to make sure font_import() and loadfonts() are only run once
if (!"Times New Roman" %in% extrafont::fonts()) {
  extrafont::font_import(prompt = FALSE)
  extrafont::loadfonts(device = "pdf", quiet = TRUE)
}

#' Summarize Positive Pairing
#'
#' @description
#' Detects proteins where the group with the larger sample size also has the
#' larger SD. Summarizes counts overall, by panel, and returns a full table.
#'
#' @param df Data frame in long format.
#' @param group_var Column with group labels (factor with two levels).
#' @param id_var Column identifying observational units.
#' @param assay_col Assay identifier (default `"Assay"`).
#' @param panel_col Panel identifier (default `"Panel"`).
#' @param outcome_col Numeric measurement column (default `"NPX"`).
#'
#' @return List with tibbles: `overall`, `by_panel`, `tbl`.
#' @export
positive_pair_summary <- function(df,
                                  group_var,          # e.g. "Treatment"
                                  id_var,             # e.g. "SampleID", "PatientID"
                                  assay_col  = "Assay",
                                  panel_col  = "Panel",
                                  outcome_col = "NPX") {

  # one row per sample
  df_id <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(assay_col, panel_col, group_var, id_var)))) |>
    dplyr::summarise(
      !!outcome_col := mean(.data[[outcome_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # sample size & SD within each arm
  grp_stats <- df_id |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(assay_col, panel_col, group_var)))) |>
    dplyr::summarise(
      n  = dplyr::n(),
      sd = stats::sd(.data[[outcome_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # wide form: n_<group>, sd_<group>
  wide_tbl <- tidyr::pivot_wider(
    grp_stats,
    names_from  = !!group_var,
    values_from = c(n, sd),
    names_glue  = "{.value}_{.name}"
  )

  # flag positive pairing
  n_cols  <- grep("^n_",  names(wide_tbl), value = TRUE)
  sd_cols <- grep("^sd_", names(wide_tbl), value = TRUE)

  tbl <- wide_tbl |>
    dplyr::rowwise() |>
    dplyr::mutate(
      grp_big_n  = sub("^n_",  "", n_cols [ which.max(dplyr::c_across(dplyr::all_of(n_cols ))) ]),
      grp_big_sd = sub("^sd_", "", sd_cols[ which.max(dplyr::c_across(dplyr::all_of(sd_cols))) ]),
      positive_pair = grp_big_n == grp_big_sd
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-grp_big_n, -grp_big_sd)

  # summaries
  overall <- tbl |>
    dplyr::summarise(
      proteins          = dplyr::n(),
      positive_pairs    = sum(positive_pair),
      negative_pairs    = proteins - positive_pairs,
      pct_positive_pair = round(100 * positive_pairs / proteins, 1)
    )

  by_panel <- tbl |>
    dplyr::group_by(.data[[panel_col]]) |>
    dplyr::summarise(
      proteins          = dplyr::n(),
      positive_pairs    = sum(positive_pair),
      pct_positive_pair = round(100 * positive_pairs / proteins, 1),
      .groups = "drop"
    )

  list(overall = overall, by_panel = by_panel, tbl = tbl)
}


#' Estimate Student-t ν from Kurtosis
#'
#' @description
#' Estimates Student-t degrees of freedom by inverting the relationship
#' between kurtosis and ν. Returns a fallback if excess kurtosis ≤ 0.
#'
#' @param df Data frame or tibble.
#' @param npx_col Column with numeric data (default `"NPX"`).
#' @param fallback Value returned if kurtosis ≤ 0 (default 10).
#'
#' @return Numeric ν estimate.
#' @export
estimate_nu_from_kurtosis <- function(df,
                                      npx_col  = "NPX",
                                      fallback = 10) {
  # require the moments package for kurtosis()
  if (!requireNamespace("moments", quietly = TRUE)) {
    stop("Please install the 'moments' package to estimate kurtosis.")
  }
  if (!npx_col %in% colnames(df)) {
    stop("Column '", npx_col, "' not found in the data frame.")
  }

  # compute excess kurtosis: E[(x-μ)^4]/σ^4 - 3
  k_excess <- moments::kurtosis(df[[npx_col]], na.rm = TRUE)

  # invert relationship for ν if positive; else fallback
  if (is.na(k_excess) || k_excess <= 0) {
    nu_est <- fallback
  } else {
    nu_est <- 6 / k_excess + 4
  }

  return(nu_est)
}


#' Compute Delta SD
#'
#' @description
#' Computes the pooled SD of per-protein differences between two groups.
#'
#' @param df Long-format data frame with NPX and grouping columns.
#' @param contrast List with elements: `col`, `lvl1`, `lvl2`.
#' @param by Grouping columns (default assay identifiers).
#'
#' @return Numeric SD of deltas.
#' @export
compute_delta_sd <- function(df,
                             contrast,
                             by = c("Assay","OlinkID","UniProt","Panel")) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::summarise(
      delta = mean(NPX[.data[[contrast$col]] == contrast$lvl1], na.rm = TRUE) -
              mean(NPX[.data[[contrast$col]] == contrast$lvl2], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::pull(delta) %>%
    stats::sd(na.rm = TRUE)
}


#' Make Priors for BRMS
#'
#' @description
#' Constructs weakly informative priors for one contrast (e.g., Disease vs control).
#' Priors are tailored to each contrast
#' using empirical SDs and kurtosis estimates.
#'
#' @param df Long-format data frame with NPX and covariates.
#' @param contrast One of `"disease_vs_control"`.
#' @param sd_treat SD for treatment effect prior (default 1).
#'
#' @return List of priors (class `brmsprior`).
#' @export
make_priors <- function(df,
                        contrast,
                        sd_treat = 1) {
  #                                                      disease_vs_control,
  # mu_npx        <-  mean(df$NPX,   na.rm = TRUE)     = 0.7825388,
  # sd_npx        <-  sd(df$NPX,   na.rm = TRUE)       = 1.643428,
  # sd_age        <-  sd(df$Age_c, na.rm = TRUE)       = 13.93882,

  # We must set priors explicitly using strings only (no variables) since
  # brms::prior does not allow variables (seriously).

  if (contrast == "disease_vs_control") {
    # ----------------------- disease vs control -----------------------
    # Compute empirical SD of the disease–control per‐protein contrasts
    sd_treat_disease <- compute_delta_sd(
      df,
      contrast = list(col = "Treatment", lvl1 = "disease", lvl2 = "control")
    )
    # sd_treat_disease = 0.9182764

    # Compute empirical SD of the Male–Female per‐protein contrasts
    sd_sex_effect <- compute_delta_sd(
      df,
      contrast = list(col = "Sex", lvl1 = "M", lvl2 = "F")
    )
    # sd_sex_effect = 0.4133289

    # Estimate between‐protein intercept SD (for the random intercept prior)
    sd_intercept_prot <- df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::summarise(mu = mean(NPX), .groups = "drop") %>%
      dplyr::pull(mu) %>%
      stats::sd(na.rm = TRUE)
    # sd_intercept_prot = 1.134259

    # Estimate ν (degrees-of-freedom) from excess kurtosis of NPX
    nu_est <- estimate_nu_from_kurtosis(df, npx_col = "NPX", fallback = 10)
    # nu_est = 4.760525

    # Assemble the priors list
    pri <- c(
      # global intercept on NPX scale: center = mean(NPX), width = 2 × SD(NPX)
      # mu_npx = 0.7825388, 2 * 1.643428 = 3.286856
      brms::prior(normal(0.7825388, 3.286856), class = "Intercept"),

      # fixed effect: disease vs control treatment slope
      # sd_treat_disease = 0.9182764
      brms::prior(normal(0, 0.9182764), class = "b", coef = "Treatmentdisease"),

      # fixed effect: age slope, ~0.2 NPX per year, scaled by SD(age)
      # sd_age = 13.93882 → 0.2 / 13.93882 = 0.014348
      brms::prior(normal(0, 0.014348), class = "b", coef = "Age_c"),

      # fixed effect: sex difference slope (M vs F)
      # sd_sex_effect = 0.4133289
      brms::prior(normal(0, 0.4133289), class = "b", coef = "SexM"),

      # residual error scale (Student-t): heavy tails, scale = SD(NPX)
      # sd_npx = 1.643428
      brms::prior(student_t(3, 0, 1.643428), class = "sigma"),

      # ---------------------- RANDOM-EFFECT SD PRIORS ----------------------
      # protein-level random INTERCEPT SD
      # sd_intercept_prot = 1.134259 → rate = 1 / 1.134259 = 0.8816329
      brms::prior(exponential(0.8816329), class = "sd",
                  group = "protein", coef = "Intercept"),

      # protein-level random SLOPE SD for disease effect
      # sd_treat_disease = 0.9182764 → rate = 1 / 0.9182764 = 1.089000
      brms::prior(exponential(1.089000), class = "sd",
                  group = "protein", coef = "Treatmentdisease"),

      # (Optional) sample-level random intercept SD
      # left at default or set explicitly, e.g. exponential(1)
      # brms::prior(exponential(1), class = "sd",
      #             group = "SampleID", coef = "Intercept"),

      # prior on degrees-of-freedom ν for Student-t tails
      # nu_est = 4.760525 → rate = 1 / 4.760525 = 0.2100609
      brms::prior(exponential(0.2100609), class = "nu")
    )
  } else {
    stop("contrast must be 'disease_vs_control'")
  }

  pri
}


#' Build BRMS Contrasts
#'
#' @description
#' Creates three contrasts (Disease vs control), generates priors, and fits hierarchical and per-protein
#' models with brms. Adds results to `res_list`.
#'
#' @param se_no_outlier SummarizedExperiment without outliers.
#' @param se_list List of SummarizedExperiments.
#' @param res_list List of results.
#' @param ... Extra arguments passed to `brms_contrast()`.
#'
#' @return List with updated `se_list` and `res_list`.
#' @export
build_brms_contrasts <- function(se_no_outlier,
                                 se_list,
                                 res_list,
                                 ...) {

  # 0. ensure cmdstanr backend
  # on Biowulf, must run module load gcc/13.2.0
  # if TRUE, CmdStan not yet registered
  if (isFALSE(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    dst <- file.path(getwd(), "cmdstan")
    dir.create(dst, showWarnings = FALSE)

    if (!file.exists(file.path(dst, "make", "program"))) {
      Sys.setenv(TBB_CXX_TYPE = "gcc")
      cmdstanr::install_cmdstan(
        dir        = dst,
        cores      = 4,
        overwrite  = TRUE                # rebuild if partially present
      )
    }
    cmdstanr::cmdstan_version()
  }
  # Set backend to cmdstanr for maximum parallel speedup
  options(brms.backend = "cmdstanr")

  # Add unique protein column to base data
  if (!"protein" %in% names(S4Vectors::metadata(se_no_outlier)$long_data)) {
    S4Vectors::metadata(se_no_outlier)$long_data <- S4Vectors::metadata(se_no_outlier)$long_data %>%
      dplyr::mutate(
        protein = interaction(Assay, Panel, OlinkID, UniProt, sep = "|", drop = TRUE)
      )
  }

  # 1. build the three SE subsets
  se_list$no_paired_no_treatment_1 <- subset_se(
      subset_se("Treatment", "treatment_1|treatment_2", "remove") %>%
      center_numerics("Age")
  )

  # 2. extract long data and re-factor levels
  df_disease_control <- S4Vectors::metadata(se_list$no_paired_no_treatment_1)$long_data %>%
                 more_levels("Treatment", "disease", "control") %>%
                 more_levels("Sex",   "M",   "F")

  # 3. formulae
  f1_g <- NPX ~ Treatment      * Age_c + Sex + (1 | SampleID) + (1 + Treatment | protein)
  f1_s <- NPX ~ Treatment      * Age_c + Sex

  # Make priors
  pri_disease_control <- make_priors(df_disease_control,
                              contrast = "disease_vs_control")

  contrasts <- list(
    disease_vs_control = list(df = df_disease_control, form_g = f1_g, form_s = f1_s,
                              pri = pri_disease_control)
  )

  # 4. fit hierarchical & per-protein models ----------------------------------
  for (nm in names(contrasts)) {
    info <- contrasts[[nm]]
    results = brms_contrast(
          df         = info$df,
          formula    = info$form_g,
          priors     = info$pri,
          by_protein = FALSE,
          iter       = 500,
          control    = list(adapt_delta   = 0.99,
                            max_treedepth = 14),
          init       = 0.2,
          backend    = "cmdstanr",
          threads    = 8,
          chains     = 4,
          cores      = 32
        )

    res_list[[paste0(nm, "_brms_1")]] <- list(
        name           = paste0(nm, "_brms_1"),
        se_object_name = "no_paired_no_treatment_1",
        column         = "Treatment",
        treatment      = "disease",
        control        = "control",
        fit            = results$fit,
        anova          = results$anova
    )
  }

  list(se_list = se_list, res_list = res_list)
}


#' Sync a SummarizedExperiment to the samples present in its metadata$long_data
#'
#' @description
#' After you finish creating/updating \code{metadata(se)$long_data} for a contrast,
#' call this to make the \strong{count data (assays)} and \strong{colData} reflect
#' exactly the samples in the updated metadata, in the same order. Optionally copy
#' newly created columns (e.g., \code{treatment_1}, \code{PatientID}) from \code{long_data}
#' into \code{colData}.
#'
#' @param se A \link[SummarizedExperiment]{SummarizedExperiment}.
#' @param add_cols Character vector of column names to copy from
#'   \code{metadata(se)[[metadata_key]]} into \code{colData(se)}. Columns not found
#'   are skipped with a warning. Defaults to \code{NULL} (copy none).
#' @param metadata_key Name of the element in \code{metadata(se)} that holds the
#'   long-format table. Default: \code{"long_data"}.
#' @param id_col Column in \code{long_data} that identifies samples/SE columns.
#'   Default: \code{"SampleID"}.
#' @param strict When \code{TRUE} (default), error if any \code{SampleID}s found
#'   in \code{long_data} are missing from \code{colnames(se)}. When \code{FALSE},
#'   silently drop missing IDs.
#' @param drop_unused_levels Drop unused factor levels in \code{colData} after sync.
#'   Default: \code{TRUE}.
#'
#' @return The updated \code{SummarizedExperiment} with:
#' \itemize{
#'   \item Columns (samples) subset/reordered to match \code{long_data[[id_col]]}.
#'   \item \code{colData} augmented with \code{add_cols}, preserving factor levels
#'         from \code{long_data}.
#' }
#' @export
sync_se_to_metadata <- function(se,
                                add_cols = NULL,
                                metadata_key = "long_data",
                                id_col = "SampleID",
                                strict = TRUE,
                                drop_unused_levels = TRUE) {
  # --- Fetch and validate metadata table ---
  ld <- S4Vectors::metadata(se)[[metadata_key]]
  if (is.null(ld) || !is.data.frame(ld)) {
    stop("metadata(se)$", metadata_key, " must be a data.frame-like object.")
  }
  if (!id_col %in% colnames(ld)) {
    stop("'", id_col, "' not found in metadata(se)$", metadata_key, ".")
  }

  # Determine sample order from long_data (first appearance)
  # Keep order as it appears in ld
  keep_ids <- as.character(ld[[id_col]])
  keep_ids <- keep_ids[!is.na(keep_ids)]
  keep_ids <- unique(keep_ids)

  # Sanity check against SE columns
  se_cols <- colnames(se)
  missing_in_se <- setdiff(keep_ids, se_cols)
  if (length(missing_in_se) > 0) {
    msg <- paste0("These ", id_col, " values are in metadata but not in SE: ",
                  paste(missing_in_se, collapse = ", "))
    if (isTRUE(strict)) stop(msg) else warning(msg, call. = FALSE)
    # If not strict, drop those from keep_ids
    keep_ids <- intersect(keep_ids, se_cols)
  }
  if (length(keep_ids) == 0L) stop("No overlapping samples between metadata and SE.")

  # Subset + reorder SE columns to match metadata order
  se <- se[, keep_ids, drop = FALSE]

  # Build a 1-row-per-sample slice of long_data to map add_cols
  # (Preserve factor levels from ld; keep first occurrence)
  first_idx <- !duplicated(ld[[id_col]])
  ld_once   <- ld[first_idx, , drop = FALSE]
  rownames(ld_once) <- as.character(ld_once[[id_col]])

  # Copy requested columns into colData (preserve factor levels)
  if (!is.null(add_cols) && length(add_cols) > 0) {
    cd_df <- as.data.frame(SummarizedExperiment::colData(se))
    for (nm in add_cols) {
      if (!nm %in% colnames(ld_once)) {
        warning("Column '", nm, "' not found in metadata; skipping.", call. = FALSE)
        next
      }
      vals <- ld_once[keep_ids, nm, drop = TRUE]
      if (is.factor(ld_once[[nm]])) {
        cd_df[[nm]] <- factor(vals, levels = levels(ld_once[[nm]]))
      } else {
        cd_df[[nm]] <- vals
      }
    }
    if (isTRUE(drop_unused_levels)) {
      for (jj in names(cd_df)) {
        if (is.factor(cd_df[[jj]])) cd_df[[jj]] <- base::droplevels(cd_df[[jj]])
      }
    }
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd_df)
  }

  se
}


#' Dynamically render a plot or table in RMarkdown
#'
#' @description Generate and knit a temporary RMarkdown chunk with specified
#'   figure dimensions and options, enabling plots and tables to render inside loops.
#'
#' @param g Expression or function returning a ggplot, plotly, or DT table.
#' @param fig.height Numeric height for the figure.
#' @param fig.width Numeric width for the figure.
#' @param dpi Numeric resolution (dots per inch).
#' @param echo Logical; display code in output.
#' @param message Logical; show messages during knit.
#' @param warning Logical; show warnings during knit.
#'
#' @author Based on code by Andreas Scharmueller, \email{andschar@@proton.me}
#'
#' @return None. Prints the dynamically created chunk to the console for knitting.
#' @export
subchunkify <- function(g,
                       fig.height = 7,
                       fig.width = 10,
                       dpi = 72,
                       echo = FALSE,
                       message = FALSE,
                       warning = FALSE) {
  g_deparsed = paste0(deparse(
    function() {g}
  ), collapse = '')

  sub_chunk = paste0("```{r sub_chunk_", floor(runif(1) * 1e5),
                     ", fig.height=", fig.height,
                     ", fig.width=", fig.width,
                     ", dpi=", dpi,
                     ", echo=", echo,
                     ", warning=", warning,
                     ", message=", message,
                     " }",
                     "\n(",
                     g_deparsed,
                     ")()",
                     "\n```
                     ")

  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}


#' Print Markdown text
#'
#' @description
#' Simple wrapper around \code{cat()} to print text as Markdown.
#' Use inside an R chunk with \code{results="asis"} enabled.
#'
#' @param ... Objects or strings to print.
#'
#' @return No return value; called for side effects.
#' @export
mdcat <- function(...){
  cat('\n\n', ..., ' \n\n', sep='', fill=1500)
}


#' Extract All Paths from Nested List
#'
#' @description
#' Recursively collects all character values from a nested list, typically used
#' to extract directory paths.
#'
#' @param x A nested list containing character vectors.
#'
#' @return Character vector of all extracted paths.
#' @noRd
get_all_paths <- function(x) {
  if (is.list(x)) {
    unlist(lapply(x, get_all_paths), use.names = FALSE)
  } else if (is.character(x)) {
    x
  } else {
    NULL
  }
}


#' Create Directories from Config
#'
#' @description
#' Recursively creates all directories defined in \code{config$output_paths}.
#'
#' @param config List with element \code{output_paths} (nested list of paths).
#'
#' @return Invisibly returns list of created or existing directories.
#' @export
make_dirs <- function(config) {
  all_paths <- get_all_paths(config$output_paths)
  invisible(lapply(all_paths, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive=TRUE)
    }
  }))
}


#' Import Sample Metadata and NPX Data
#'
#' @description
#' Reads the sample metadata TSV and NPX data file. If \code{data_path} contains
#' "pea_test_data.tsv", reads directly with \pkg{readr}; otherwise uses
#' \code{OlinkAnalyze::read_NPX()}.
#'
#' @param data_path Path to NPX data file (TSV or Olink raw format).
#' @param sampletable_path Path to sample metadata TSV.
#'
#' @return List with \code{sampletable} and \code{data}.
#' @export
import_data <- function(data_path, sampletable_path) {
  sample_metadata <- readr::read_tsv(sampletable_path)

  # If pea_test_data.tsv is used, assume preformatted test data;
  # otherwise read raw data using the OlinkAnalyze package.
  if (stringr::str_detect(data_path, "pea_test_data.tsv") |
      stringr::str_detect(data_path, "elisa_test_data.tsv")) {
    long_data <- readr::read_tsv(data_path)
  } else {
    long_data <- OlinkAnalyze::read_NPX(data_path)
  }
  return(list(sampletable = sample_metadata, data = long_data))
}


#' Subset SummarizedExperiment
#'
#' @description
#' Subsets a \linkS4class{SummarizedExperiment} either by:
#' \itemize{
#'   \item Regex match on a column (\code{column}, \code{pattern}, \code{keep_or_remove})
#'   \item Arbitrary filter expression (\code{filter_phrase})
#' }
#' Exactly one mode must be supplied.
#'
#' @param se A \linkS4class{SummarizedExperiment}.
#' @param column Column in \code{colData(se)} to match.
#' @param pattern Regex pattern to filter on.
#' @param keep_or_remove `"keep"` or `"remove"`.
#' @param filter_phrase Dplyr filter expression string.
#' @param rm_qc_warn If TRUE, drop QC and Assay warnings after filtering.
#'
#' @return Filtered \linkS4class{SummarizedExperiment}.
#' @export
subset_se <- function(se,
                      column           = NULL,
                      pattern          = NULL,
                      keep_or_remove   = NULL,
                      filter_phrase    = NULL,
                      rm_qc_warn       = FALSE) {

  # Validate arguments
  mode_filter    <- !is.null(filter_phrase)
  mode_regex     <- !is.null(column) && !is.null(pattern) && !is.null(keep_or_remove)

  if (xor(mode_filter, mode_regex) == FALSE) {
    stop("You must supply exactly one of:\n",
         "  - filter_phrase\n",
         "  - column + pattern + keep_or_remove")
  }

  # Filter by arbitrary expression
  if (mode_filter) {
    expr <- rlang::parse_expr(filter_phrase)

    # colData
    coldf <- as.data.frame(SummarizedExperiment::colData(se))
    keep_samples <- eval(expr, envir = coldf)
    if (!is.logical(keep_samples) || length(keep_samples) != nrow(coldf)) {
      stop("'filter_phrase' must evaluate to a logical vector of length ncol(se)")
    }
    # Filter data and colData
    se_subset <- se[, keep_samples]

    # Filter metadata
    md       <- S4Vectors::metadata(se)$long_data
    se_md    <- dplyr::filter(md, !!expr)

  } else {
  # Filter by single-column regex
    hits <- grepl(pattern, SummarizedExperiment::colData(se)[[column]])
    if (keep_or_remove == "keep") {
      se_subset <- se[, hits]
      se_md     <- S4Vectors::metadata(se)$long_data %>%
                     dplyr::filter(grepl(pattern, .data[[column]]))

    } else if (keep_or_remove == "remove") {
      se_subset <- se[, !hits]
      se_md     <- S4Vectors::metadata(se)$long_data %>%
                     dplyr::filter(!grepl(pattern, .data[[column]]))

    } else {
      stop("`keep_or_remove` must be 'keep' or 'remove', not '", keep_or_remove, "'")
    }
  }

  # Optionally drop QC warnings
  if (isTRUE(rm_qc_warn)) {
    se_md <- dplyr::filter(se_md, QC_Warning != "WARN")
    se_md <- dplyr::filter(se_md, Assay_Warning != "WARN")
  }

  # Attach metadata and validate
  S4Vectors::metadata(se_subset)$long_data <- se_md
  check_se_object(se_subset)

  se_subset
}


#' Add Patient ID Column
#'
#' @description
#' Derives a patient identifier by stripping the trailing letter from
#' \code{SampleID} and adds it as \code{PatientID} to both \code{colData} and
#' \code{metadata$long_data}.
#'
#' @param se A \linkS4class{SummarizedExperiment}.
#' @param sampleid_col Column with sample IDs (default `"SampleID"`).
#' @param pairid_col Name for new patient column (default `"PatientID"`).
#' @param overwrite If TRUE, overwrite existing column.
#'
#' @return Modified SummarizedExperiment.
#' @export
add_PatientID <- function(se,
                        sampleid_col = "SampleID",
                        pairid_col   = "PatientID",
                        overwrite    = TRUE) {

  # Sanity checks
  if (!sampleid_col %in% colnames(SummarizedExperiment::colData(se))) {
    stop(crayon::red("Column '", sampleid_col, "' not found in colData(se)."))
  }

  if (pairid_col %in% colnames(SummarizedExperiment::colData(se)) && !overwrite) {
    stop(crayon::red("Column '", pairid_col,
             "' already exists. Set overwrite = TRUE to replace it."))
  }

  # Derive PairID: only drop a single trailing letter, leave digits intact
  sample_ids <- as.character(SummarizedExperiment::colData(se)[[sampleid_col]])
  pair_ids   <- sub("[A-Za-z]$", "", sample_ids)

  # Add to colData
  SummarizedExperiment::colData(se)[[pairid_col]] <- pair_ids

  # Add to metadata$long_data
  md <- S4Vectors::metadata(se)$long_data

  if (!sampleid_col %in% colnames(md)) {
    stop(crayon::red("Column '", sampleid_col,
             "' not found in metadata$long_data."))
  }
  if (pairid_col %in% colnames(md) && !overwrite) {
    stop(crayon::red("Column '", pairid_col,
             "' already exists in metadata$long_data. ",
             "Set overwrite = TRUE to replace it."))
  }

  md[[pairid_col]] <- sub("[A-Za-z]$", "", md[[sampleid_col]])
  S4Vectors::metadata(se)$long_data <- md

  return(se)
}


#' Center Numeric Columns
#'
#' @description
#' Adds mean-centered versions of numeric columns in both colData and metadata.
#'
#' @param se A \linkS4class{SummarizedExperiment}.
#' @param numeric_cols Character vector of numeric column names.
#' @param overwrite If TRUE, overwrite existing centered columns.
#'
#' @return Modified SummarizedExperiment with new \code{*_c} columns.
#' @export
center_numerics <- function(se,
                            numeric_cols,
                            overwrite = TRUE) {
  # Sanity checks for colData
  cd <- SummarizedExperiment::colData(se)
  missing_cols <- setdiff(numeric_cols, colnames(cd))
  if (length(missing_cols) > 0) {
    stop(crayon::red("Column(s) not found in colData(se): ", paste(missing_cols, collapse=", ")))
  }

  # Prepare metadata long_data
  md <- S4Vectors::metadata(se)$long_data
  missing_md <- setdiff(numeric_cols, colnames(md))
  if (length(missing_md) > 0) {
    stop(crayon::red("Column(s) not found in metadata$long_data: ", paste(missing_md, collapse=", ")))
  }

  # For each numeric column, center and add new column
  for (col in numeric_cols) {
    new_col <- paste0(col, "_c")
    # Check overwrite
    if (new_col %in% colnames(cd) && !overwrite) {
      stop(crayon::red("Column '", new_col,
               "' already exists in colData(se). Set overwrite = TRUE to replace it."))
    }
    if (new_col %in% colnames(md) && !overwrite) {
      stop(crayon::red("Column '", new_col,
               "' already exists in metadata$long_data. Set overwrite = TRUE to replace it."))
    }

    # Center in colData
    vals <- cd[[col]]
    if (!is.numeric(vals)) {
      stop(crayon::red("Column '", col, "' must be numeric to center."))
    }
    cd[[new_col]] <- as.numeric(scale(vals, scale = FALSE))

    # Center in metadata long_data
    md_vals <- md[[col]]
    if (!is.numeric(md_vals)) {
      stop(crayon::red("Column '", col,
               "' in metadata$long_data must be numeric to center."))
    }
    md[[new_col]] <- as.numeric(scale(md_vals, scale = FALSE))
  }

  # Assign back and return
  SummarizedExperiment::colData(se) <- cd
  S4Vectors::metadata(se)$long_data <- md
  return(se)
}


#' Add Baseline-Centered Column
#'
#' @description
#' For each patient, extracts baseline (pre-treatment) values of a numeric
#' variable, centers them, and adds as \code{*_0_c}.
#'
#' @param se A \linkS4class{SummarizedExperiment}.
#' @param numeric_cols Column(s) to baseline-center.
#' @param prepost_col Column indicating pre/post status (default `"PrePost"`).
#' @param patientid_col Column with patient IDs (default `"PatientID"`).
#' @param overwrite If TRUE, overwrite existing columns.
#'
#' @return Modified SummarizedExperiment.
#' @export
add_centered_baseline <- function(se,
                                  numeric_cols,
                                  prepost_col   = "PrePost",
                                  patientid_col = "PatientID",
                                  overwrite     = TRUE) {

  # support passing a single name, or a vector
  numeric_cols <- as.character(numeric_cols)

  cd <- SummarizedExperiment::colData(se)
  md <- S4Vectors::metadata(se)$long_data

  # basic sanity checks
  reqs_cd <- c(numeric_cols, prepost_col, patientid_col)
  reqs_md <- reqs_cd
  miss_cd <- setdiff(reqs_cd, colnames(cd))
  miss_md <- setdiff(reqs_md, colnames(md))
  if (length(miss_cd)) stop("Missing in colData: ", paste(miss_cd, collapse=", "))
  if (length(miss_md)) stop("Missing in metadata$long_data: ", paste(miss_md, collapse=", "))

  for (col in numeric_cols) {
    new_col <- paste0(col, "_0_c")
    if (new_col %in% colnames(cd) && !overwrite) {
      stop("Column '", new_col, "' already exists in colData(se); use overwrite=TRUE")
    }

    # get each patient’s pre‐treatment mean
    sel_pre   <- cd[[prepost_col]] == "pre"
    vals_pre  <- cd[[col]][sel_pre]
    ids_pre   <- cd[[patientid_col]][sel_pre]
    baselines <- tapply(vals_pre, ids_pre, mean)

    # map it back and mean‐center
    cd[[new_col]] <- as.numeric(
      scale(baselines[ as.character(cd[[patientid_col]]) ], scale = FALSE)
    )
    md[[new_col]] <- as.numeric(
      scale(baselines[ as.character(md[[patientid_col]]) ], scale = FALSE)
    )
  }

  SummarizedExperiment::colData(se)             <- cd
  S4Vectors::metadata(se)$long_data <- md
  se
}


#' Prepare SummarizedExperiment
#'
#' @description
#' Merges NPX data with metadata, removes controls, and builds a
#' SummarizedExperiment.
#'
#' @param olink_data Long NPX data frame.
#' @param sampletable Sample metadata table.
#' @param elisa If TRUE, use `"ELISA"` instead of `"NPX"`.
#'
#' @return SummarizedExperiment with assays, colData, and metadata.
#' @export
prepare_data_for_se_list <- function(olink_data,
                                     sampletable,
                                     elisa = FALSE) {
  sampletable <- sampletable %>%
    filter(!stringr::str_detect(SampleID, "CONTROL_SAMPLE"))
  olink_data <- olink_data %>%
    filter(!stringr::str_detect(SampleID, "CONTROL_SAMPLE"))

  if ("Treatment" %in% names(olink_data)) {
    olink_data <- olink_data %>% dplyr::select(-Treatment)
  }

  olink_data <- olink_data %>%
    dplyr::left_join(sampletable, by = "SampleID")

  if (isTRUE(elisa)) {
    value_col <- "ELISA"
  } else {
    value_col <- "NPX"
  }
  olink_wide <- olink_data %>%
    tidyr::pivot_wider(
      id_cols = c(Assay, Panel),
      names_from = SampleID,
      values_from = dplyr::all_of(value_col)
    )

  npx_matrix <- as.matrix(olink_wide[, sampletable$SampleID])
  rownames(npx_matrix) <- paste(olink_wide$Assay, olink_wide$Panel)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(NPX = npx_matrix),
    colData = sampletable,
    metadata = list(long_data = olink_data)
  )

  check_se_object(se)

  return(se)
}


#' Collapse Groups to Two Levels
#'
#' @description
#' Collapses a grouping variable into treatment vs control.
#'
#' @param df Data frame.
#' @param column Column to collapse.
#' @param treatment Values treated as `"treatment"`.
#' @param control Values treated as `"control"`.
#'
#' @return Data frame with two-level factor.
#' @export
two_levels <- function(df, column, treatment, control) {
  df <- df %>%
    dplyr::filter(.data[[column]] %in% c(treatment, control)) %>%
    dplyr::mutate(!!column := factor(
      if_else(.data[[column]] %in% treatment, "treatment", "control",),
      levels = c("treatment", "control")
    ))
  return(df)
}


#' Keep Specific Factor Levels
#'
#' @description
#' Keeps specified treatment and control levels as factor levels without
#' collapsing them.
#'
#' @param df Data frame.
#' @param column Grouping column.
#' @param treatment Treatment levels.
#' @param control control levels.
#'
#' @return Data frame with filtered factor levels.
#' @export
more_levels <- function(df, column, treatment, control) {
  df <- df %>%
    dplyr::filter(.data[[column]] %in% c(treatment, control)) %>%
    dplyr::mutate(!!column := factor(.data[[column]], levels = c(treatment, control))) # Reference level goes first
  return(df)
}


#' Extract parameters for a contrast
#'
#' @description
#' Collects inputs needed to run statistical tests for a contrast, including
#' data, factor/covariate settings, and analysis options.
#'
#' @param contrast_obj List specifying the contrast (fields: \code{se_object_name},
#'   \code{name}, \code{sampletable_column}, \code{treatment}, \code{control},
#'   \code{posthocs}, \code{collapse}, \code{dependent}, \code{pair_id},
#'   \code{covariates}, \code{subject}, \code{run_stats}, etc.).
#' @param se_list Named list of \linkS4class{SummarizedExperiment} objects.
#'
#' @return A named list of parameters ready for downstream testing.
#' @export
get_params <- function(contrast_obj, se_list) {
  # Get the se object for this contrast
  se_object_name <- contrast_obj$se_object_name
  se <- se_list[[se_object_name]]
  df_work <- S4Vectors::metadata(se)$long_data

  params <- rlang::list2(
    name           = contrast_obj$name,
    se_object_name = se_object_name,
    df_work        = df_work,
    alpha          = contrast_obj$alpha,
    posthocs       = contrast_obj$posthocs %||% FALSE,
    collapse       = contrast_obj$collapse %||% FALSE,
    dependent      = contrast_obj$dependent %||% FALSE,
    alt_lfc_col    = contrast_obj$alt_lfc_col,
    pair_id        = contrast_obj$pair_id,
    covariates     = contrast_obj$covariates,
    observed       = contrast_obj$observed,
    formula        = contrast_obj$formula,
    subject        = contrast_obj$subject,
    column         = contrast_obj$sampletable_column,
    treatment      = contrast_obj$treatment,
    control        = contrast_obj$control,
    more_factors   = contrast_obj$more_factors,
    emm            = contrast_obj$emm,
    report_lfc     = contrast_obj$report_lfc,
    robust         = contrast_obj$robust,
    run_stats      = contrast_obj$run_stats
  )

  # Messages of collapsed factors and dependent data setting outcome
  if (isTRUE(params$collapse) &&
      (length(params$treatment) > 1 || length(params$control) > 1)) {
    message(crayon::red("collapse = TRUE -> multiple conditions will be collapsed"))
  }
  if (isTRUE(params$dependent)) {
    message(crayon::red("dependent = TRUE -> using Friedman or LMM for repeated measures"))
  }

  return(params)
}


#' Run core statistical tests
#'
#' @description
#' Executes the specified set of core tests (ANOVA (may contain
#' configured posthocs), linear mixed models,
#' t-test, Wilcoxon, Kruskal–Wallis, Friedman, ordinal regression).
#'
#' @param params Parameter list from \code{get_params()}.
#'
#' @return Named list of test results (tables).
#' @export
choose_frameworks <- function(params) {
  run_stats    <- params$run_stats
  df_work      <- params$df_work
  df_orig      <- params$df_orig
  column       <- params$column
  treatment    <- params$treatment
  control      <- params$control
  factor       <- params$factor
  more_factors <- params$more_factors
  pair_id      <- params$pair_id
  covariates   <- params$covariates
  dependent    <- params$dependent
  subject      <- params$subject
  alpha        <- params$alpha

  # List to store results from this contrast_obj/set of params
  res <- list()

  if (is.null(run_stats) || "anova" %in% run_stats) {
    message(crayon::green("Running anova ..."))

    anova_args <- list(
      long_df     = params$df_work,
      formula     = params$formula,
      treatment   = params$treatment,
      control     = params$control,
      variable    = params$column,
      observed    = params$observed,
      alpha       = params$alpha,
      report_lfc  = params$report_lfc,
      emm         = params$emm,
      alt_lfc_col = params$alt_lfc_col,
      lmer        = FALSE
    )
    # Return both anova table and estimated marginal means summary for interaction terms.
    # These results will be saved in res_list[[contrast_name]]
    all_res <- do.call(afex_test, anova_args)
    res$anova   <- all_res$anova
    res$emm_int <- all_res$emm_int
  }

  if (is.null(run_stats) || "lmer" %in% run_stats) {
    message(crayon::green("Running lmer ..."))

    anova_args <- list(
      long_df     = params$df_work,
      formula     = params$formula,
      treatment   = params$treatment,
      control     = params$control,
      variable    = params$column,
      alpha       = params$alpha,
      report_lfc  = params$report_lfc,
      emm         = params$emm,
      alt_lfc_col = params$alt_lfc_col,
      lmer        = TRUE
    )
    all_res <- do.call(afex_test, anova_args)
    res$anova   <- all_res$anova
    res$emm_int <- all_res$emm_int
  }

  if (isTRUE(factor) && grepl("ttest|wilcox", run_stats)) {
    # collapse to two levels
    olink_df_collapsed <- two_levels(df_orig, column, treatment, control)
  }

  ## t-test
  if (is.null(run_stats) || "ttest" %in% run_stats) {
    ttest_args <- list(
      long_df  = olink_df_collapsed,
      variable = column,
      alpha    = alpha
    )
    if (!is.null(pair_id)) ttest_args$pair_id <- pair_id
    res$ttest <- do.call(t_test, ttest_args)
  }

  ## Wilcoxon
  if (is.null(run_stats) || "wilcox" %in% run_stats) {
    wilcox_args <- list(
      long_df  = olink_df_collapsed,
      variable = column,
      alpha    = alpha
    )
    if (!is.null(pair_id)) wilcox_args$pair_id <- pair_id
    res$wilcox <- do.call(wilcox_test, wilcox_args)
  }

  ## Non-parametric omnibus (Kruskal–Wallis or Friedman)
  if (is.null(run_stats) || "nonparam" %in% run_stats) {
    np_args <- list(
      long_df    = df_factor,
      variable   = column,
      alpha      = alpha,
      dependence = dependent
    )
    # for Friedman, subject ID is required
    if (!is.null(subject) || dependent) {
      np_args$subject <- subject %||% pair_id
    }
    res$nonparam <- do.call(kruskal_test, np_args)
  }

  ## Ordinal regression (not currently supported)
  if (is.null(run_stats) || "ordinal" %in% run_stats) {
    ordinal_args <- list(
      long_df          = df_factor,
      variable         = column,
      covariates       = covariates,
      alpha            = params$alpha,
      outcome          = params$outcome %||% "NPX",
      return.covariates= params$return.covariates %||% FALSE,
      verbose          = params$verbose %||% FALSE
    )
    res$ordinal <- do.call(ordinal_regression, ordinal_args)
  }

  return(res)
}


#' Run all tests for a contrast
#'
#' @description
#' Wrapper to run tests for one contrast, returning
#' results together with the parameter list.
#'
#' @param contrast_obj List specifying a contrast.
#' @param se_list Named list of \linkS4class{SummarizedExperiment} objects.
#'
#' @return List of test results and parameters for the contrast.
#' @export
tests <- function(contrast_obj, se_list) {
  params <- get_params(contrast_obj, se_list)

  res <- choose_frameworks(params)

  # If collapse is TRUE, normalize treatment/control labels
  if (isTRUE(params$collapse)) {
    params$treatment <- "treatment"
    params$control   <- "control"
  }

  params$data <- params$df_factor
  params$df_factor <- NULL
  result <- rlang::list2(
    params = params,    # Includes $df_factor as $data
    !!!res  # Unpack res like **kwargs in python so we don't have another
            # layer to iterate over when using results for future steps
  )

  return(result)
}


#' Run tests for multiple contrasts
#'
#' @description
#' Iterates over a list of contrast specifications and runs \code{tests()}
#' for each.
#'
#' @param contrast_list Named list of contrast specification lists.
#' @param se_list Named list of \linkS4class{SummarizedExperiment} objects.
#'
#' @return Named list of results for all contrasts.
#' @export
run_tests <- function(contrast_list, se_list) {
  results <- list()
  for (contrast_name in names(contrast_list)) {
    contrast_obj <- contrast_list[[contrast_name]]
    message(crayon::green(paste0("Running tests on ", contrast_name, ": ", contrast_obj$name)))
    results[[contrast_obj$name]] <- tests(contrast_obj, se_list)
  }
  return(results)
}


#' Validate contrast parameter object(s)
#'
#' @description
#' Checks that contrast parameter objects are correctly structured, names are
#' unique, and required fields are present. Also validates factor or numeric
#' column types against the SummarizedExperiment data.
#'
#' @param params_list Single contrast list or list of contrast lists.
#' @param se_list Named list of \linkS4class{SummarizedExperiment} objects
#' @param sensitivity Boolean if TRUE, don't worry about the Excel sheet name
#' length requirement.
#'
#' @return Invisibly returns the original \code{params_list} if validation passes.
#' @export
validate_parameters <- function(params_list, se_list, sensitivity = FALSE) {
  # Detect if a *single* contrast object was supplied (it owns a `name` field)
  if (is.list(params_list) && !is.null(params_list$name)) {
    objs   <- list(params_list)       # wrap for internal checks
  } else {
    objs   <- params_list             # already a list of objects
  }

  validate_params(objs, sensitivity)
  validate_columns(objs, se_list)

  msg <- paste("Validated contrast parameter object(s):",
               paste(vapply(objs, `[[`, character(1), "name"), collapse = ", "))
  message(crayon::green(msg))

  return(params_list)  # return *exactly* what the user passed in
}


#' Validate essential parameters
#'
#' @inheritParams validate_parameters
#' @noRd
validate_params <- function(objs, sensitivity) {
  contrast_names <- vapply(objs, `[[`, character(1), "name")

  if (anyDuplicated(contrast_names)) {
    dupes <- unique(contrast_names[duplicated(contrast_names)])
    stop(crayon::red("Duplicated contrast names detected: ",
             paste(dupes, collapse = ", ")))
  }

  if (isFALSE(sensitivity) && any(nchar(contrast_names) > 31)) {
    long <- contrast_names[nchar(contrast_names) > 31]
    stop(crayon::red("Excel sheet names must be ≤31 characters. Please shorten: ",
             paste(long, collapse = ", ")))
  }

  base_required <- c("name", "se_object_name", "sampletable_column", "factor")
  for (p in objs) {
    missing_base <- setdiff(base_required, names(p))
    if (length(missing_base)) {
      stop(crayon::red("Missing essential parameters in contrast '", p$name,
               "': ", paste(missing_base, collapse = ", ")))
      if (isTRUE(p$factor)) {
        # When factor=TRUE, require treatment and control in params
        req_fc <- c("treatment", "control")
        missing_fc <- setdiff(req_fc, names(p))
        if (length(missing_fc)) {
          stop(crayon::red("Missing factor parameters in contrast '", p$name,
                   "': ", paste(missing_fc, collapse = ", ")))
        }
      }
    }
  }
}


#' Validate columns & optionally check factor levels or numeric type
#'
#' @inheritParams validate_parameters
#' @noRd
validate_columns <- function(objs, se_list) {
  for (p in objs) {
    se_object_name <- p$se_object_name
    se <- se_list[[se_object_name]]
    check_se_object(se)
    long_data <- S4Vectors::metadata(se)$long_data

    # Column present?
    if (!p$sampletable_column %in% colnames(long_data)) {
      stop(crayon::red("Column '", p$sampletable_column,
               "' not found in Olink data of contrast '", p$name, "'."))
    }

    col_vals <- long_data[[p$sampletable_column]]

    if (isTRUE(p$factor)) {
      # Verify column is character
      if (!is.factor(col_vals)) {
        stop(crayon::red(
          "Column '", p$sampletable_column,
          "' in contrast '", p$name,
          "' must be factor when factor=TRUE."
        ))
      }
      # Verify treatment & control values exist (as specified in params)
      missing_vals <- setdiff(c(p$treatment, p$control), unique(col_vals))
      if (length(missing_vals)) {
        stop(crayon::red(
          "Missing values in column '", p$sampletable_column,
          "' for contrast '", p$name, "': ",
          paste(missing_vals, collapse = ", ")
        ))
      }
      # Optional post-hoc: only if provided
      if (!is.null(p$posthocs) && !p$posthocs %in% colnames(long_data)) {
        stop(crayon::red("Post-hoc column '", p$posthocs,
                 "' not found in Olink data of contrast '", p$name, "'."))
      }
    } else {
      # Continuous covariate: ensure numeric
      if (!is.numeric(col_vals)) {
        stop(crayon::red(
          "Column '", p$sampletable_column,
          "' in contrast '", p$name,
          "' must be numeric when factor=FALSE."
        ))
      }
      # Skip treatment/control validation
    }
  }
}

#' Validate a SummarizedExperiment object
#'
#' @description
#' Ensures an object is a \linkS4class{SummarizedExperiment}, with non-empty
#' colData, an NPX assay, and metadata containing \code{long_data}.
#'
#' @param se Object to validate.
#'
#' @return Invisible \code{NULL}; errors if validation fails.
#' @export
check_se_object <- function(se) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("The object `se` is not a SummarizedExperiment.")
  }

  if (ncol(SummarizedExperiment::colData(se)) == 0) {
    stop("The SummarizedExperiment `se` does not contain any `colData`.")
  }

  if (!"NPX" %in% names(SummarizedExperiment::assays(se))) {
    stop("The SummarizedExperiment `se` does not contain an assay named `NPX`.")
  }

  if (is.null(S4Vectors::metadata(se)$long_data)) {
    stop("The metadata of `se` does not contain `long_data`.")
  }

  message(crayon::green("All checks passed, SE object is valid."))
}


#' Map internal test codes to names
#'
#' @description
#' Converts short test identifiers (e.g. "ttest") into human-readable names.
#'
#' @param test_type Short test identifier string.
#'
#' @return Character string with a descriptive test name.
#' @export
get_test_name <- function(test_type) {
  test_names <- list(
    ttest            = "t-test",
    wilcox           = "Wilcoxon",
    anova            = "ANOVA",
    anova_posthoc    = "ANOVA (post hoc)",
    lmer             = "Linear mixed model",
    lmer_posthoc     = "Linear mixed model (post hoc)",
    nonparam         = "Non-parametric",
    nonparam_posthoc = "Non-parametric (post hoc)",
    ordinal          = "Ordinal regression",
    ordinal_posthoc  = "Ordinal regression (post hoc)"
  )
  return(test_names[[test_type]] %||% test_type)
}


#' Extract significant proteins for UpSet plots
#'
#' @description
#' Selects proteins that are both significant and directionally changed
#' in a given contrast/test. Direction is defined using the `estimate` column.
#'
#' @param res_list Nested results list from `run_tests()` and `add_estimate_column()`.
#' @param contrast Contrast name within `res_list`.
#' @param direction One of `"up"`, `"down"`, or `"changed"`.
#' @param test_type Test type name (e.g. `"ttest"`, `"anova"`).
#'
#' @return Tibble with one column `Assay` listing significant proteins.
#' @export
get_sig <- function(res_list, contrast, direction,
                    test_type)
{
  res   <- res_list[[contrast]][[test_type]]
  alpha <- res_list[[contrast]]$params$alpha        # contrast-specific FDR

  res <- switch(direction,
                up      = res[ res$estimate >= 0, ],
                down    = res[ res$estimate <  0, ],
                changed = res,
                stop("direction must be 'up', 'down', or 'changed'"))

  sel <- !is.na(res$Adjusted_pval) & res$Adjusted_pval < alpha
  sig_vec <- paste(res$Assay[sel], res$Panel[sel], sep = "_")

  # ALWAYS return a tibble, even when empty
  tibble::tibble(Assay = unique(sig_vec))
}


#' Convert list of Assay tibbles to incidence matrix
#'
#' @description
#' Converts a named list of tibbles (each with an `Assay` column) into a
#' binary 0/1 incidence matrix suitable for UpSet plotting.
#'
#' @param lst Named list of tibbles with an `Assay` column.
#'
#' @return Data frame of 0/1 indicators (rows = assays, columns = list names).
#' @export
from_list_with_names <- function(lst) {
  elements <- unique(unlist(lst, use.names = FALSE))

  mat <- vapply(lst,
                function(x) as.integer(elements %in% x),
                integer(length(elements)))
  df  <- as.data.frame(mat, row.names = elements, check.names = FALSE)

  df[rowSums(df) > 0, , drop = FALSE]
}


#' Format results for export or display
#'
#' @description
#' Cleans and formats a results data frame:
#' - Drops all-NA columns and unnecessary analysis columns.
#' - Rounds numeric columns or formats p-values in scientific notation.
#' - Supports special handling of therapeutic-reversal tables.
#'
#' @param df Results data frame.
#' @param one_effect Logical, if `TRUE` assumes single effect.
#' @param drop_enrich_cols Logical, drop enrichment-specific columns.
#' @param scientific_notation Logical, format p-values in scientific notation.
#'
#' @return Data frame formatted for printing or Excel export.
#' @export
format_results <- function(df,
                           one_effect = FALSE,
                           drop_enrich_cols = FALSE,
                           scientific_notation = FALSE) {

  # Drop any columns that contain only NAs
  df <- df %>%
    dplyr::select(where(~ !all(is.na(.x))))

  # ---------------------------------------------------------------------------
  # Branch: therapeutic_reversal-style tables (joined two-contrast results)
  # Detect by presence of estimate_/sig_ cols + reversal flag
  # ---------------------------------------------------------------------------
  is_tr_table <- (
    "reversal" %in% names(df) &&
    any(startsWith(names(df), "estimate_")) &&
    any(startsWith(names(df), "sig_"))
  )

  if (is_tr_table) {
    # Remove plotting-only helpers
    df2 <- df %>% dplyr::select(-dplyr::any_of(c("y_plot", "shape")))

    # Capture estimate_*/sig_* columns in the order they appear
    est_cols <- names(df2)[startsWith(names(df2), "estimate_")]
    sig_cols <- names(df2)[startsWith(names(df2), "sig_")]

    # Aim for a column order similar to single-contrast results,
    # but with the paired-contrast fields up front after IDs.
    base_cols <- c(
      "Assay","GeneID","Panel","UniProt","OlinkID","Assay_Panel",
      est_cols, sig_cols,
      "status","significant_in_both","reversal","reversal_sig_both"
    )

    formatted_df <- df2 %>%
      dplyr::select(-Threshold, dplyr::any_of(base_cols), dplyr::everything())

    # Numeric formatting: round estimates to 3 decimals
    if (length(est_cols)) {
      formatted_df <- formatted_df %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(est_cols), ~ round(.x, 3)))
    }

    if (isTRUE(drop_enrich_cols)) {
      formatted_df <- formatted_df %>%
        dplyr::select(-dplyr::any_of(c("ID","leading_edge","core_enrichment","geneID")))
    }

    return(as.data.frame(formatted_df))
  }

  # ---------------------------------------------------------------------------
  # Default branch: standard single-contrast results
  # ---------------------------------------------------------------------------
  formatted_df <- as.data.frame(df) %>%
    # Remove intermediary analysis columns from final Excel data
    dplyr::select(-dplyr::any_of(c(
      "resids","AD_padj","LT_padj","AIC","BIC","AICc","AD_pval","LT_pval"
    ))) %>%
    # Reorder final Excel data columns
    dplyr::select(
      dplyr::any_of(c(
        "Assay","GeneID","Panel","UniProt","OlinkID",
        "Effect","MajorRegions","BrainStructures","Cerebullum","Brainstem","Cerebrum",
        "pes","p.value","Adjusted_pval","Threshold",
        "log2FC","estimate","formula","df","MSE","SE","F","Chisq","t.ratio",
        "SW_pval","SW_padj","FT_pval","FT_padj"
      )),
      dplyr::everything() # Raw NPX / sample-level columns, if present
    )

  if (isTRUE(scientific_notation)) {
    # Include only these base columns for embedding in the report
    base_cols <- c(
      "Assay","Panel","Effect","p.value","Adjusted_pval","Threshold",
      "log2FC","estimate","SW_padj","FT_padj","formula"
    )

    formatted_df <- formatted_df %>%
      # Keep only the “base” columns, dropping every SampleID column
      dplyr::select(dplyr::any_of(base_cols)) %>%
      # Format numeric columns (Warning: becomes character)
      dplyr::mutate(
        dplyr::across(
          where(is.numeric),
          ~ if (cur_column() %in%
                c("SW_padj","FT_padj","p.value","Adjusted_pval","pvalue",
                  "p.adjust","adj.pval","qvalue")) {
              sprintf("%.2e", .)
            } else {
              round(., 3)
            }
        )
      ) %>%
      tibble::remove_rownames()
  }

  if (isTRUE(drop_enrich_cols)) {
    formatted_df <- formatted_df %>%
      dplyr::select(-dplyr::any_of(c("ID", "leading_edge", "core_enrichment", "geneID")))
  }

  return(formatted_df)
}


#' Prepare metadata and counts from SummarizedExperiment
#'
#' @description
#' Extracts `colData`, long-format `olink_data`, and NPX assay counts
#' from a SummarizedExperiment, rounding numeric columns.
#'
#' @param se SummarizedExperiment object.
#'
#' @return List with `metadata`, `olink_data`, and `counts`.
#' @export
prep_data <- function(se) {
  check_se_object(se)

  colData   <- SummarizedExperiment::colData(se)
  counts    <- SummarizedExperiment::assays(se)$NPX
  olink_data <- S4Vectors::metadata(se)$long_data

  # Convert counts: always add Protein_Panel column
  counts_df <- as.data.frame(counts)
  colnames(counts_df) <- make.unique(colnames(counts_df))
  counts_df <- tibble::rownames_to_column(counts_df, var = "Protein_Panel")

  # Build list
  data_list <- list(
    metadata   = colData,
    olink_data = olink_data,
    counts     = counts_df
  )

  # Clean metadata and olink_data formatting
  for (name in c("metadata", "olink_data")) {
    df <- tibble::as_tibble(data_list[[name]]) %>%
      dplyr::slice(gtools::mixedorder(SampleID)) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))) %>%
      tibble::remove_rownames() %>%
      as.data.frame()
    data_list[[name]] <- df
  }

  return(data_list)
}


#' Retrieve top significant proteins
#'
#' @description
#' Gets top 5 up-regulated and top 5 down-regulated proteins by adjusted
#' p-value.
#'
#' @param df Results data frame with `estimate`, `Adjusted_pval`, `Threshold`.
#'
#' @return List with `OlinkID` and `Assay` vectors.
#' @export
top_sig <- function(df) {
  top5_pos <- df %>%
    dplyr::filter(estimate >= 0, Threshold == 'Significant') %>%
    dplyr::arrange(Adjusted_pval) %>%
    dplyr::slice_head(n = 5)

  top5_neg <- df %>%
    dplyr::filter(estimate < 0, Threshold == 'Significant') %>%
    dplyr::arrange(Adjusted_pval) %>%
    dplyr::slice_head(n = 5)

  combined_df <- bind_rows(top5_pos, top5_neg) %>%
    dplyr::arrange(Adjusted_pval)

  return(list(OlinkID = combined_df$OlinkID, Assay = combined_df$Assay))
}


#' Select proteins for plot labeling
#'
#' @description
#' Returns proteins to label on plots. Prefers user-specified labels in
#' `config$label_on_plots$disease_vs_control`, falling back to top significant proteins.
#'
#' @param result Results data frame with `Assay` and significance columns.
#' @param config Configuration list (may contain `label_on_plots`).
#'
#' @return Data frame of rows to label.
#' @export
label_plot <- function(result, config) {
  # Compute the fallback list of top significant assays
  top_assays <- top_sig(result)$Assay

  # If user-specified labels exist, try to use those
  user_labels <- config$label_on_plots$disease_vs_control %||% character()
  if (length(user_labels) > 0) {
    found    <- intersect(user_labels, result$Assay)
    missing  <- base::setdiff(user_labels, result$Assay)

    if (length(found) > 0) {
      if (length(missing) > 0) {
        warning(
          crayon::red(
            "The following assays from config$label_on_plots$disease_vs_control were not found in the results: ",
            paste(missing, collapse = ", "),
            ".\nLabeling available assays: ",
            paste(found, collapse = ", ")
          )
        )
      }
      selected_assays <- found

    } else {
      # none of the user labels were found
      fallback <- if (length(top_assays) > 0) {
        paste(top_assays, collapse = ", ")
      } else {
        "None"
      }
      warning(
        crayon::red(
          "None of the assays in config$label_on_plots$disease_vs_control were found in the results.",
          if (fallback != "None") {
            paste0(" Falling back to top significant proteins: ", fallback, ".")
          } else {
            " No significant proteins available to label."
          }
        )
      )
      selected_assays <- top_assays
    }

  } else {
    # No custom labels specified: use the top significant proteins
    selected_assays <- top_assays
  }

  # Return only the rows for the selected assays
  result %>% dplyr::filter(Assay %in% selected_assays)
}


#' Report number of significant proteins across contrasts (with up/down counts)
#'
#' @description
#' Counts significant proteins (`Threshold == "Significant"`) for every
#' contrast and test type, and adds directional counts using the `estimate`
#' column: `up` = count of significant rows with `estimate > 0`;
#' `down` = count of significant rows with `estimate < 0`.
#'
#' @param res_list Named list of contrast results.
#'
#' @return Tibble with columns `contrast`, `test_type`, `n_significant`, `up`, `down`.
#' @export
report_num_significant <- function(res_list) {
  purrr::map_dfr(
    names(res_list),
    function(contrast) {
      purrr::map_dfr(
        base::setdiff(names(res_list[[contrast]]), "params"),
        function(test_type) {
          df <- res_list[[contrast]][[test_type]]

          # Significant flag
          sig <- df$Threshold == "Significant"

          # Totals
          n_sig <- sum(sig, na.rm = TRUE)
          up    <- sum(sig & df$estimate > 0, na.rm = TRUE)
          down  <- sum(sig & df$estimate < 0, na.rm = TRUE)

          tibble::tibble(
            contrast      = contrast,
            test_type     = test_type,
            n_significant = n_sig,
            up            = up,
            down          = down
          )
        }
      )
    }
  )
}


#' Get number of significant proteins for one contrast/test
#'
#' @description
#' Returns the count of significant proteins for a specific contrast and test type.
#'
#' @inheritParams report_num_significant
#' @param contrast Contrast name.
#' @param test_type Test type name.
#'
#' @return Integer count.
#' @export
get_num_sig <- function(res_list, contrast, test_type) {
  report_num_significant(res_list) %>%
    dplyr::filter(
      .data$contrast  == {{contrast}},
      .data$test_type == {{test_type}}
    ) %>%
    dplyr::pull(n_significant)
}


#' Get number of proteins lost after adjustment
#'
#' @description
#' Counts proteins significant in one contrast but not in another
#' (e.g. ANOVA vs ANCOVA).
#'
#' @param res_list Results list.
#' @param contrast_a Baseline contrast.
#' @param contrast_b Adjusted contrast.
#' @param test_type Test type to check (default `"anova"`).
#'
#' @return Integer count of “lost” proteins.
#' @export
get_num_sig_lost <- function(res_list,
                             contrast_a,
                             contrast_b,
                             test_type = "anova") {

  # Helper to pull the tibble and filter significant rows
  pull_sig <- function(contrast) {
    if (!contrast %in% names(res_list)) return(NULL)
    tbl <- res_list[[contrast]][[test_type]]
    if (is.null(tbl)) return(NULL)
    dplyr::filter(tbl, .data$Threshold == "Significant")
  }

  a_sig <- pull_sig(contrast_a)
  b_sig <- pull_sig(contrast_b)

  if (is.null(a_sig) || nrow(a_sig) == 0) return(integer(0))

  # Identify an ID column present in both tables
  id_col <- intersect(
    c("OlinkID", "UniProt", "Protein", "Gene", "Accession"),
    intersect(colnames(a_sig), if (!is.null(b_sig)) colnames(b_sig) else character())
  )[1]

  if (is.na(id_col)) {
    warning("No common ID column found; returning integer(0).")
    return(integer(0))
  }

  if (is.null(b_sig) || nrow(b_sig) == 0) {
    # Nothing significant after adjustment -> all are “lost”
    return(nrow(a_sig))
  }

  lost <- dplyr::anti_join(a_sig, b_sig, by = id_col)
  nrow(lost)
}


#' Export lost-protein tables to Excel
#'
#' @description
#' For contrast pairs, exports proteins significant in ANOVA but not
#' ANCOVA to Excel.
#'
#' @param res_list Results list.
#' @param contrast_pairs Named list of contrast pairs.
#' @param config Config with `output_paths`.
#'
#' @return Invisibly, named character vector of Excel paths.
#' @export
export_lost_proteins <- function(res_list, contrast_pairs, config) {
  outdir <- config$output_paths$differential_expression$results

  wb_paths <- character()

  for (name in names(contrast_pairs)) {
    anova_con  <- contrast_pairs[[name]]["anova"]
    ancova_con <- contrast_pairs[[name]]["ancova"]

    a_sig <- res_list[[anova_con]][["anova"]]  %>%
      dplyr::filter(Threshold == "Significant")

    c_sig <- res_list[[ancova_con]][["anova"]] %>%
      dplyr::filter(Threshold == "Significant")

    if (is.null(a_sig) || nrow(a_sig) == 0) next

    id_col <- if ("OlinkID" %in% colnames(a_sig)) "OlinkID" else "UniProt"

    lost_df <- dplyr::anti_join(a_sig, c_sig, by = id_col)
    if (nrow(lost_df) == 0) next

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Lost_ANOVA_hits")
    openxlsx::writeData(
      wb, sheet = 1,
      x   = format_results(
        lost_df,
        drop_enrich_cols     = FALSE,
        scientific_notation  = FALSE
      ),
      withFilter = TRUE
    )

    file_name <- paste0(name, "_lost_proteins.xlsx")
    path      <- file.path(outdir, file_name)
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    wb_paths[name] <- path
  }

  invisible(wb_paths)
}


#' Print lost-protein tables with download links
#'
#' @description
#' Displays lost-protein tables as DT datatables with Excel download links.
#'
#' @inheritParams export_lost_proteins
#' @export
print_lost_proteins <- function(res_list, contrast_pairs, config) {
  wb_paths <- export_lost_proteins(res_list, contrast_pairs, config)

  for (name in names(contrast_pairs)) {
    mdcat(paste("###", name))

    anova_con  <- contrast_pairs[[name]]["anova"]
    ancova_con <- contrast_pairs[[name]]["ancova"]

    a_sig <- res_list[[anova_con]][["anova"]]  %>%
      dplyr::filter(Threshold == "Significant")

    c_sig <- res_list[[ancova_con]][["anova"]] %>%
      dplyr::filter(Threshold == "Significant")

    id_col <- if ("OlinkID" %in% colnames(a_sig)) "OlinkID" else "UniProt"

    lost_df <- dplyr::anti_join(a_sig, c_sig, by = id_col)

    mdcat("#### Lost proteins (in ANOVA not in ANCOVA)")

    if (!is.na(wb_paths[name])) {
      mdcat(paste0("- [Download Excel workbook](", wb_paths[name], ")"))
    }

    if (nrow(lost_df) == 0) {
      mdcat("No proteins were lost after Age adjustment.")
    } else {
      lost_fmt <- format_results(
        lost_df,
        drop_enrich_cols    = FALSE,
        scientific_notation = TRUE
      )
      subchunkify(DT::datatable(lost_fmt))
    }
  }
}


#' Run all tests in parallel
#'
#' @description
#' Runs `run_tests()` in parallel if enabled in config, otherwise sequentially.
#'
#' @param contrast_list Contrast parameters list.
#' @param se_list SummarizedExperiment list.
#' @param config Config list with `parallel` and `cores`.
#'
#' @return Results list.
#' @export
parallel_test_runner <- function(contrast_list, se_list, config) {
  if (isTRUE(config$parallel)) {
    ncores <- config$cores

    # Helper that runs *one* contrast at a time. We wrap the single‐element list
    # around contrast_list[[cn]] so that run_tests still returns a “list” of length 1.
    .run_one_contrast <- function(cn) {
      single_param <- contrast_list[cn]         # e.g. list( contrast_list[[cn]] )
      out <- run_tests(single_param, se_list)
      # run_tests() always returns a named list; extract its first (and only) element
      out[[1]]
    }

    # Dispatch each contrast to its own core:
    res_list <- parallel::mclapply(
      X       = names(contrast_list),
      FUN     = .run_one_contrast,
      mc.cores = ncores
    )

    names(res_list) <- names(contrast_list)
    return(res_list)

  } else {

    #Run all tests on all contrasts in contrast_list
    res_list <- run_tests(contrast_list, se_list)
    return(res_list)
  }
}


#' Keep primary effects and best hit per assay
#'
#' @description
#' Filters result tables to the main effect and sorts by ascending Adjusted_pval.
#'
#' @param res_list Results list.
#' @param test_types Character vector of test table names (default `"anova"`).
#'
#' @return Results list with filtered tables.
#' @export
keep_primary_effects <- function(res_list, test_types = c("anova")) {
  # Exclude contrasts lacking Effect column
  res_list <- res_list[setdiff(names(res_list), "therapeutic_reversal")]

  stopifnot(is.list(res_list), is.character(test_types), length(test_types) >= 1)

  lapply(res_list, function(contrast_res) {
    # Identify which column holds the main effect name
    main_effect <- contrast_res$params$column

    for (tt in test_types) {
      if (!tt %in% names(contrast_res)) {
        stop(
          "Contrast '", names(contrast_res), "' does not contain a table named '",
          tt, "'."
        )
      }

      # Filter to the main effect, then sort by ascending Adjusted_pval
      contrast_res[[tt]] <- contrast_res[[tt]] %>%
        dplyr::filter(Effect == main_effect) %>%
        sort_results()
    }

    contrast_res
  })
}


#' Record standard-assumption violations
#'
#' @description
#' For each contrast, identifies assays failing assumption tests based on
#' adjusted-p columns.
#'
#' @param res_list Results list.
#' @param column Adjusted-p column name.
#'
#' @return Nested list of violation tables by contrast/test.
#' @export
record_standard_assumption_violations <- function(res_list, column) {
  out_list <- setNames(vector("list", length(res_list)), names(res_list))
  for (contrast_name in names(res_list)) {
    contrast_res <- res_list[[contrast_name]]
    alpha        <- contrast_res$params$alpha

    # pull the ANOVA table every time
    tbl <- contrast_res[["anova"]]
    if (!all(c("Assay", "Panel", column) %in% names(tbl))) {
      stop("missing columns in ANOVA for ", contrast_name)
    }

    viol_df <- tbl %>%
      dplyr::filter(.data[[column]] < alpha) %>%
      dplyr::select(Assay, Panel) %>%
      dplyr::distinct()

    # override with a single-element list named "anova"
    out_list[[contrast_name]] <- list(anova = viol_df)
  }
  out_list
}


#' Find overlap between paired and unpaired treatment_1 contrasts
#'
#' @description
#' Identifies significant overlapping assays between paired and unpaired treatment_1
#' contrasts and exports them to Excel.
#'
#' @param res_list Results list with treatment_1 contrasts.
#' @param config Config with output path.
#' @param filename Excel filename.
#'
#' @return Invisibly returns overlap table.
#' @export
find_treatment_1_overlap <- function(res_list, config, filename = "overlap_unpaired_vs_paired.xlsx") {

  # Extract significant results for unpaired and paired treatment_1 contrasts
  sig_unpaired <- res_list[[3]]$anova %>%
    dplyr::filter(Effect == res_list[[3]]$params$column, Threshold == "Significant")
  sig_paired <- res_list[[4]]$anova %>%
    dplyr::filter(Effect == res_list[[4]]$params$column, Threshold == "Significant")

  # Find overlapping Assay + Panel
  overlap <- dplyr::inner_join(
    sig_unpaired,
    sig_paired,
    by = c("Assay", "Panel"),
    suffix = c(".unpaired", ".paired")
  )

  # Arrange columns for output
  overlap_table <- overlap %>%
    dplyr::select(
      Assay, Panel,
      Effect.unpaired, pes.unpaired, log2FC.unpaired, estimate.unpaired, p.value.unpaired, Adjusted_pval.unpaired, Threshold.unpaired,
      Effect.paired, pes.paired, log2FC.paired, estimate.paired, p.value.paired, Adjusted_pval.paired, Threshold.paired
    )

  # Construct the output file path
  out_path <- file.path(config$output_paths$differential_expression$results, filename)

  # Write to Excel
  openxlsx::write.xlsx(overlap_table, file = out_path, asTable = TRUE, overwrite = TRUE)

  message("Overlap table written to: ", out_path)
  invisible(overlap_table)
}


#' Assess sample fit with Cook’s distance or residuals
#'
#' @description
#' Generic driver to compute sample diagnostics (Cook’s distance or residuals)
#' and return plots per contrast.
#'
#' @param res_list Results list.
#' @param metric_fun Function computing diagnostics.
#' @param plot_fun Function plotting diagnostics.
#' @param extra_pkgs Required packages.
#'
#' @return Named list of plots.
#' @export
assess_sample_fit <- function(res_list, metric_fun, plot_fun, extra_pkgs = character()) {
  stopifnot(requireNamespace("dplyr"), requireNamespace("afex"))
  if (length(extra_pkgs)) for (pkg in extra_pkgs) stopifnot(requireNamespace(pkg))

  plots <- vector("list", length(res_list))
  names(plots) <- names(res_list)

  for (contrast in names(res_list)) {
    params   <- res_list[[contrast]]$params
    df_long  <- params$df_work
    fml_chr  <- params$formula
    test_ty  <- params$run_stats

    # remove incomplete rows based on predictors in the model
    form_obj <- stats::as.formula(fml_chr)
    vars     <- all.vars(form_obj)
    outcome  <- vars[1]
    preds    <- setdiff(vars, outcome)

    idx <- stats::complete.cases(df_long[preds])
    if (any(!idx)) {
      warning("Removed ", sum(!idx), " incomplete rows.")
      df_long <- df_long[idx, ] %>% droplevels()
    }

    # keep model info on df for plotting functions
    attr(df_long, "formula")   <- fml_chr
    attr(df_long, "test_type") <- test_ty

    # unique assay_panel names
    assay_panel_names <- df_long %>%
      dplyr::group_by(Assay, Panel) %>%
      dplyr::group_keys() %>%
      dplyr::mutate(assay_panel = paste(Assay, Panel, sep = "_")) %>%
      dplyr::pull(assay_panel)

    # compute metric per protein
    metrics <- df_long %>%
      dplyr::group_by(Assay, Panel) %>%
      dplyr::group_map(~ metric_fun(.x, fml_chr, test_ty))
    names(metrics) <- assay_panel_names

    # generate and collect the plot object (contrast is now just a label)
    plots[[contrast]] <- plot_fun(metrics, df_long, contrast)
  }

  plots
}



# helper: choose a grouping variable from df + formula
.get_group_var <- function(df) {
  # An explicit Group column wins
  if ("Group" %in% names(df)) return("Group")

  # Or fall back to first predictor in formula that is present in df
  fml_chr <- attr(df, "formula", exact = TRUE)
  if (!is.null(fml_chr)) {
    vars <- all.vars(stats::as.formula(fml_chr))
    preds <- setdiff(vars, vars[1])
    preds <- preds[preds %in% names(df)]
    if (length(preds) > 0) return(preds[1])
  }

  # If no formula is found
  NULL
}


#' @noRd
get_cooks <- function(df_one, fml, test_ty) {
  fit <- if (test_ty == "anova") {
    afex::aov_car(as.formula(fml), data = df_one, factorize = FALSE, return = "lm")
  } else {
    afex::mixed(as.formula(fml), data = df_one, progress = FALSE, verbose = FALSE)$full_model
  }
  cd <- stats::cooks.distance(fit)
  names(cd) <- df_one$SampleID
  cd
}


#' @noRd
get_residuals <- function(df_one, fml, test_ty) {
  rownames(df_one) <- df_one$SampleID
  fit <- if (test_ty == "anova") {
    afex::aov_car(as.formula(fml), data = df_one, factorize = FALSE, return = "lm")
  } else {
    afex::mixed(as.formula(fml), data = df_one, progress = FALSE, verbose = FALSE)$full_model
  }
  res <- stats::residuals(fit)
  tibble::tibble(SampleID = df_one$SampleID, residual = res)
}


#' @noRd
plot_cooks <- function(cook_list, df, contrast) {
  # assemble Cook's distance matrix
  cook_mat   <- do.call(cbind, cook_list)
  sample_ids <- sort(unique(df$SampleID))
  rownames(cook_mat) <- sample_ids
  cook_mat <- cook_mat[sample_ids, , drop = FALSE]

  # determine grouping variable and values
  group_var <- .get_group_var(df)
  if (!is.null(group_var)) {
    grp <- df[[group_var]][match(sample_ids, df$SampleID)]
  } else {
    grp <- rep("all_samples", length(sample_ids))
  }
  grp <- as.factor(grp)

  # automatic color mapping for groups
  levs <- levels(grp)
  cols <- grDevices::rainbow(length(levs))
  names(cols) <- levs
  col_map <- cols

  row_anno <- ComplexHeatmap::rowAnnotation(
    Group = grp,
    col   = list(Group = col_map)
  )

  ComplexHeatmap::Heatmap(
    cook_mat,
    name             = "Cook",
    left_annotation  = row_anno,
    show_row_names   = TRUE,
    cluster_columns  = FALSE,
    show_column_names = FALSE,
    row_names_gp     = grid::gpar(fontsize = 6)
  )
}


#' @noRd
plot_residuals <- function(res_list, df, contrast) {
  # bind the residuals
  res_tbl <- dplyr::bind_rows(res_list)

  # determine grouping variable
  group_var <- .get_group_var(df)
  if (!is.null(group_var)) {
    grp_vals <- df[[group_var]][match(res_tbl$SampleID, df$SampleID)]
  } else {
    grp_vals <- rep("all_samples", nrow(res_tbl))
  }
  res_tbl$Group <- factor(grp_vals)

  # summarize per sample
  sample_diag <- res_tbl %>%
    dplyr::group_by(SampleID, Group) %>%
    dplyr::summarise(
      mean_res = mean(residual),
      sd_res   = sd(residual),
      prop_hi  = mean(residual > 0),
      prop_low = mean(residual < 0),
      .groups  = "drop"
    )

  # color mapping for whatever groups exist
  levs <- levels(sample_diag$Group)
  cols <- grDevices::rainbow(length(levs))
  names(cols) <- levs
  col_map <- cols

  p <- ggplot2::ggplot(
    sample_diag,
    ggplot2::aes(
      x = prop_hi,
      y = mean_res,
      color = Group,
      text = paste0(
        "SampleID: ", SampleID, "<br>",
        "Mean residual: ", round(mean_res, 3), "<br>",
        "Prop. positive residuals: ", round(prop_hi, 3), "<br>",
        "Group: ", Group
      )
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2) +
    ggplot2::labs(
      x     = "Proportion of positive residuals",
      y     = paste(
        "Mean residual across",
        length(unique(paste(df$Assay, df$Panel))),
        "proteins"
      ),
      color = "Group"
    ) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::theme_minimal()

  plotly::ggplotly(p, tooltip = "text")
}


# wrappers
assess_sample_fit_with_cooks <- function(res_list) {
  assess_sample_fit(res_list, get_cooks, plot_cooks, extra_pkgs = c("ComplexHeatmap"))
}

assess_sample_fit_with_residuals <- function(res_list) {
  assess_sample_fit(res_list, get_residuals, plot_residuals, extra_pkgs = c("ggplot2"))
}


#' Parallel assessment of sample fit with Cook's distance
#'
#' @description
#' Runs \code{assess_sample_fit_with_cooks()} in parallel if
#' \code{config$parallel} is TRUE, otherwise sequentially.
#'
#' @param res_list Results list, with one element per contrast.
#' @param config Config list with \code{parallel} (logical) and
#'   \code{cores} (integer, number of cores).
#'
#' @return Named list of Cook’s distance heatmap objects (one per contrast).
#' @export
parallel_assess_sample_fit_with_cooks <- function(res_list, config) {
  if (isTRUE(config$parallel)) {
    ncores <- config$cores
    .run_one <- function(cn) {
      single <- res_list[cn]; names(single) <- cn
      assess_sample_fit(single, get_cooks, plot_cooks, extra_pkgs = "ComplexHeatmap")[[cn]]
    }
    plots <- parallel::mclapply(names(res_list), .run_one, mc.cores = ncores)
    names(plots) <- names(res_list)
    plots
  } else {
    assess_sample_fit_with_cooks(res_list)
  }
}


#' Parallel assessment of sample fit with residuals
#'
#' @description
#' Runs \code{assess_sample_fit_with_residuals()} in parallel if
#' \code{config$parallel} is TRUE, otherwise sequentially.
#'
#' @inheritParams parallel_assess_sample_fit_with_cooks
#'
#' @return Named list of residual diagnostic plots (one per contrast).
#' @export
parallel_assess_sample_fit_with_residuals <- function(res_list, config) {
  if (isTRUE(config$parallel)) {
    ncores <- config$cores
    .run_one <- function(cn) {
      single <- res_list[cn]; names(single) <- cn
      assess_sample_fit(single, get_residuals, plot_residuals, extra_pkgs = "ggplot2")[[cn]]
    }
    plots <- parallel::mclapply(names(res_list), .run_one, mc.cores = ncores)
    names(plots) <- names(res_list)
    plots
  } else {
    assess_sample_fit_with_residuals(res_list)
  }
}


#' Append GeneID annotations
#'
#' @description
#' Maps assay symbols to ENSEMBL GeneIDs and appends them to result tables.
#'
#' @param res_list Results list.
#'
#' @return Results list with GeneIDs.
#' @export
append_geneid <- function(res_list) {
  # Collect all unique assay symbols actually present in res_list
  # There are 2843 unique protein (Assay) names
  assays <- unique(unlist(lapply(res_list, function(df) df$anova$Assay), use.names = FALSE))
  # We dont lose any after checking for NAs or empty strings
  assays <- assays[!is.na(assays) & nzchar(assays)]

  # Build SYMBOL -> ENSEMBL mapping only for those assays
  mapping_tbl <- AnnotationDbi::select(
    x       = org.Hs.eg.db,
    keys    = assays,
    columns = c("SYMBOL", "ENSEMBL"),
    keytype = "SYMBOL"
  )

  # Clean & reduce to one GeneID per Assay (take the first)
  mapping_tbl <- mapping_tbl |>
    base::subset(!is.na(ENSEMBL)) |>
    base::unique() |>
    tibble::as_tibble() |>
    dplyr::rename(Assay = SYMBOL, GeneID = ENSEMBL) |>
    dplyr::group_by(Assay) |>
    dplyr::summarise(GeneID = dplyr::first(GeneID), .groups = "drop")

  # Append GeneID to each result table without changing row counts
  res_list <- append_cols(res_list, mapping_tbl, by = "Assay")

  # Reorder columns/rows as your helper already does
  res_list <- sort_all_results(res_list)

  res_list
}


#' Sort all contrast result tables
#'
#' @description
#' Iterates through each element of a results list and applies
#' \code{sort_results()} to the \code{$anova} table, ensuring
#' consistent ordering across all contrasts.
#'
#' @param res_list Named list of contrast results. Each element must
#'   contain a tibble \code{$anova}.
#'
#' @return The same list structure, with each \code{$anova} table sorted.
#' @export
sort_all_results <- function(res_list) {
  for (nm in names(res_list)) {
    res_list[[nm]]$anova <-
      res_list[[nm]]$anova %>%
      sort_results()
  }
  res_list
}


#' Append annotation columns to all contrasts
#'
#' @description
#' Adds extra annotation columns (e.g. gene IDs) to each contrast’s
#' \code{$anova} table by performing a left join with a provided
#' annotation tibble.
#'
#' @param res_list Named list of contrast results. Each element must
#'   contain a tibble \code{$anova}.
#' @param annot_tbl Tibble or data frame with annotation columns to join.
#' @param by Character name(s) of key columns to match between
#'   \code{res_list[[nm]]$anova} and \code{annot_tbl}.
#'
#' @return The same results list, with annotation columns appended
#'   to each \code{$anova} table.
#' @export
append_cols <- function(res_list, annot_tbl, by) {
  for (nm in names(res_list)) {
    res_list[[nm]]$anova <-
      res_list[[nm]]$anova %>%
      dplyr::left_join(annot_tbl, by = by)
  }
  res_list
}


#' Sort results consistently
#'
#' @description
#' Sorts result tables by adjusted p-value and/or effect size.
#' Special handling for contrasts with multiple `estimate_*` columns
#' (e.g. therapeutic reversal).
#'
#' @param results Results data frame.
#'
#' @return Tibble sorted by significance and effect size.
#' @export
sort_results <- function(results) {
  results %>%
    # Clean-up: drop temporary columns created by emmeans
    dplyr::select(-matches("(eff|_estimate|est_raw|p_raw|_emms_pval|_emms_padj)$")) %>%

    # Re-order the remaining columns for easier viewing
    dplyr::select(
      Assay, any_of("GeneID"), Panel,
      any_of(c("Effect", "pes", "log2FC", "estimate")),
      starts_with("estimate_"),
      any_of(c("p.value", "Adjusted_pval", "Threshold", "formula")),
      any_of(c("SW_padj", "AD_padj", "LT_padj", "FT_padj")),
      everything()
    ) %>%

    # Sort with branching logic
    {
      estimate_cols <- grep("^estimate_", names(.), value = TRUE)

      if (length(estimate_cols) >= 2) {
        # Therapeutic reversal case: combine multiple estimates
        dplyr::mutate(., rank_abs = rowSums(abs(dplyr::across(all_of(estimate_cols))))) %>%
          dplyr::arrange(dplyr::desc(rank_abs))
      } else if ("estimate" %in% names(.)) {
        dplyr::arrange(., Adjusted_pval, dplyr::desc(abs(estimate)))
      } else {
        dplyr::arrange(., Adjusted_pval)
      }
    }
}


#' Build therapeutic-reversal results
#'
#' @description
#' Joins two contrasts into a therapeutic-reversal result object using
#' `estimate_plot()`. Flags concordance and reversal events.
#'
#' @param res Results list.
#' @param config Config with output paths.
#' @param contrasts Two contrast names.
#' @param color Threshold column (default `"Threshold"`).
#' @param y.lim Numeric y-axis limits.
#' @param export_table Logical, export joined table to file.
#' @param export_path Optional explicit file path.
#'
#' @return Result object shaped like `res_list` entries.
#' @export
make_therapeutic_reversal_results <- function(res,
                                              config,
                                              contrasts,
                                              long_df,
                                              color         = "Threshold",
                                              y.lim         = c(-2, 2),
                                              export_table  = FALSE,
                                              export_path   = NULL) {
  plot_tbl <- estimate_plot(
    res           = res,
    config        = config,
    contrasts     = contrasts,
    label_concord = FALSE,
    label_sig_both= FALSE,
    plotly        = FALSE,
    color         = color,
    y.lim         = y.lim,
    label_n       = 0,
    do_plot       = FALSE,
    export_table  = export_table,
    export_path   = export_path,
    return_table  = TRUE
  )

  estA_name <- paste0("estimate_", contrasts[1])
  estB_name <- paste0("estimate_", contrasts[2])
  sigA_name <- paste0("sig_",       contrasts[1])
  sigB_name <- paste0("sig_",       contrasts[2])
  adjA_name <- paste0("Adjusted_pval_", contrasts[1])
  adjB_name <- paste0("Adjusted_pval_", contrasts[2])
  pA_name   <- paste0("p.value_",       contrasts[1])
  pB_name   <- paste0("p.value_",       contrasts[2])

  required_cols <- c(
    "Assay","OlinkID","UniProt","Panel",
    estA_name, estB_name, sigA_name, sigB_name, adjA_name, adjB_name, pA_name, pB_name
  )
  missing_cols  <- setdiff(required_cols, names(plot_tbl))
  if (length(missing_cols) > 0) stop("Missing expected columns in plotting table: ", paste(missing_cols, collapse = ", "))

  res_tbl <- plot_tbl %>%
    dplyr::mutate(
      rank_abs          = abs(.data[[estA_name]]) + abs(.data[[estB_name]]),
      reversal          = (.data[[estA_name]] * .data[[estB_name]]) < 0,
      reversal_sig_both = significant_in_both & reversal,
      Threshold         = dplyr::if_else(reversal_sig_both, "Significant", "Non-significant")
    )

  # Add raw NPX columns wide-by-SampleID
  res_tbl <- add_raw_npx(results = res_tbl, long_df = long_df) %>%
    dplyr::arrange(dplyr::desc(reversal_sig_both), dplyr::desc(rank_abs))

  res_obj <- list(
    params = list(
      name        = "therapeutic_reversal",
      contrast_A  = contrasts[1],
      contrast_B  = contrasts[2],
      column      = paste0(contrasts[1], " vs ", contrasts[2]),
      color       = color,
      source      = "estimate_plot(joined_table)+add_raw_npx"
    ),
    anova = res_tbl
  )
  return(res_obj)
}


#' Add raw NPX values to results
#'
#' Pivot NPX (or ELISA) values wide by SampleID and join onto results.
#'
#' @param results Results tibble.
#' @param long_df Long-format data with NPX/ELISA.
#' @return Results tibble with added sample columns.
#' @export
add_raw_npx <- function(results, long_df) {
  if ("NPX" %in% names(long_df)) {
    data_col <- "NPX"
  } else {
    data_col <- "ELISA"
  }
  # Ensure NPX and SampleID exist
  if (!all(c("Assay", "OlinkID", "UniProt", "Panel", "SampleID", data_col) %in% names(long_df))) {
    stop("`long_df` must contain columns: Assay, OlinkID, UniProt, Panel, SampleID, and NPX/ELISA.", call. = FALSE)
  }

  # Select only the columns needed to build the wide NPX table
  npx_wide <- long_df %>%
    dplyr::select(Assay, OlinkID, UniProt, Panel, SampleID, dplyr::all_of(data_col)) %>%
    # Pivot so each SampleID becomes its own column
    tidyr::pivot_wider(
      names_from  = SampleID,
      values_from = dplyr::all_of(data_col)
    )

  # Join the wide NPX data onto the results tibble
  results_with_npx <- results %>%
    dplyr::left_join(
      npx_wide,
      by = c("Assay", "OlinkID", "UniProt", "Panel")
    )

  # Return the augmented results
  results_with_npx
}
