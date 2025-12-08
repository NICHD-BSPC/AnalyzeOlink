#' Compute Spearman Correlations
#'
#' @description
#' For each metadata variable in `vars`, computes Spearman rank correlations
#' between NPX values and the variable of interest. Age is treated separately:
#' if the variable name matches `"Age"` or any name in `age_vars`, the test is
#' run using `long_data_age`; otherwise `long_data_others`. P-values are
#' adjusted per-variable using the Benjamini–Hochberg method (FDR).
#'
#' @param long_data_age    Long-format NPX data frame used for age-based tests.
#' @param long_data_others Long-format NPX data frame used for all other tests.
#' @param vars             Character vector of metadata variable names to test.
#' @param alpha            Numeric significance threshold for adjusted p-values.
#' @param age_vars         Character vector of variable names to treat as age
#'                         (default `"Age"`).
#'
#' @return A named list of tibbles (one per variable) with columns:
#'   * `Assay`, `OlinkID`, `UniProt`, `Panel`
#'   * `spearman_coefficient` — correlation coefficient (ρ)
#'   * `p.value` — raw p-value
#'   * `adj.pval` — FDR-adjusted p-value
#'   * `Threshold` — significance label ("Significant"/"Non-significant")
#'
#' @export
compute_spearman <- function(long_data_age,
                             long_data_others,
                             vars,
                             alpha,
                             age_vars = "Age") {

  calc_one <- function(df, var) {
    df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::summarise(
        spearman_coefficient =
          stats::cor(NPX, .data[[var]], method = "spearman", use = "complete.obs"),
        p.value = stats::cor.test(
          NPX, .data[[var]], method = "spearman", exact = FALSE
        )$p.value,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        adj.pval  = stats::p.adjust(p.value, method = "fdr"),
        Threshold = ifelse(adj.pval < alpha, "Significant", "Non-significant")
      ) %>%
      dplyr::arrange(adj.pval)
  }

  res <- lapply(vars, function(var) {
    data_to_use <- if (var %in% age_vars) long_data_age else long_data_others
    calc_one(data_to_use, var)
  })
  names(res) <- vars
  res
}


#' Format Correlation Results
#'
#' @description
#' Cleans and formats a single correlation results tibble:
#' - Drops all-NA columns.
#' - Renames `adj.pval` → `Adjusted_pval`.
#' - Orders columns consistently.
#' - Optionally converts numeric values to strings for compact reporting.
#' - Optionally drops enrichment-related columns.
#'
#' @param df A correlation results tibble.
#' @param alpha Numeric significance threshold (unused; kept for compatibility).
#' @param drop_enrich_cols Logical; drop enrichment columns if present.
#' @param scientific_notation Logical; if TRUE, display p-values in scientific
#'   notation and round coefficients.
#'
#' @return A tibble with standardized columns, ready for export or reporting.
#' @export
format_correlations_results <- function(df,
                                        alpha = 0.1,
                                        drop_enrich_cols = FALSE,
                                        scientific_notation = FALSE
) {
  # Note: `alpha` kept for call-compatibility but not used (Threshold already computed)

  df <- df %>%
    tibble::as_tibble() %>%
    # Drop columns that contain only NA
    dplyr::select(where(~ !all(is.na(.x)))) %>%
    # Rename adjusted p-value column if present (keep p.value as-is)
    dplyr::rename_with(~ "Adjusted_pval", .cols = dplyr::any_of("adj.pval"))

  # Preferred correlation-specific column order (others follow)
  preferred_order <- c(
    "Assay", "Panel", "UniProt", "OlinkID",
    "spearman_coefficient", "p.value", "Adjusted_pval", "Threshold"
  )

  formatted_df <- df %>%
    dplyr::select(dplyr::any_of(preferred_order), dplyr::everything()) %>%
    tibble::remove_rownames()

  # Optional: compact display for reports (formats numerics as strings)
  if (isTRUE(scientific_notation)) {
    base_cols <- c("Assay", "Panel", "spearman_coefficient", "p.value", "Adjusted_pval", "Threshold")
    formatted_df <- formatted_df %>%
      dplyr::select(dplyr::any_of(base_cols)) %>%
      dplyr::mutate(
        dplyr::across(
          where(is.numeric),
          ~ if (dplyr::cur_column() %in% c("p.value", "Adjusted_pval")) {
              sprintf("%.2e", .)
            } else {
              round(., 3)
            }
        )
      )
  }

  # Optional: drop enrichment-related columns if present
  if (isTRUE(drop_enrich_cols)) {
    formatted_df <- formatted_df %>%
      dplyr::select(-dplyr::any_of(c("ID", "leading_edge", "core_enrichment", "geneID")))
  }

  formatted_df
}


#' Export Spearman Correlations to Excel
#'
#' @description
#' Combines all per-variable correlation tables into one long-format tibble,
#' adds a `variable` column, and writes a single Excel file with one sheet
#' (`"correlations"`).
#'
#' @param spearman_res Named list of tibbles from `compute_spearman()`.
#' @param out_dir Directory where the workbook will be saved.
#'
#' @return Invisibly returns the path to the written workbook.
#' @export
export_spearman <- function(spearman_res, out_dir) {
  # Format each table and tag with its variable name
  combined_long <- purrr::imap_dfr(
    spearman_res,
    function(df, var) {
      format_correlations_results(df, scientific_notation = FALSE) |>
        dplyr::mutate(variable = var, .before = 1)
    }
  )

  # Write one-sheet workbook
  file_out <- file.path(out_dir, "spearman_correlations.xlsx")
  openxlsx::write.xlsx(
    x          = list(correlations = combined_long),
    file       = file_out,
    rowNames   = FALSE,
    overwrite  = TRUE
  )

  invisible(file_out)
}


#' Print Spearman Correlation Tables in RMarkdown
#'
#' @description
#' Creates a tabset of DT datatables, one per metadata variable. Only significant
#' results (adj.pval < alpha) are shown. Includes a download link to the full
#' Excel workbook.
#'
#' @param spearman_res Named list of tibbles from `compute_spearman()`.
#' @param out_xlsx Path to the Excel workbook from `export_spearman()`.
#'
#' @return Invisibly returns NULL. Called for side-effects.
#' @export
print_spearman_tables <- function(spearman_res, out_xlsx) {
  mdcat(paste0("- [Download all Excel tables](", out_xlsx, ")"))
  for (var in names(spearman_res)) {
    mdcat(paste0("### ", var))
    df <- spearman_res[[var]]
    subchunkify(
      DT::datatable(
        format_correlations_results(df, scientific_notation = TRUE) %>%
        dplyr::filter(Threshold == "Significant")
      ),
      fig.height = 5, fig.width = 8
    )
  }
  invisible(NULL)
}


#' Print Top N Significant Correlation Plots
#'
#' @description
#' For each metadata variable and each Assay, produces scatterplots of NPX vs
#' the variable for the top N most significant correlations (optionally filtered
#' by sign). Saves each plot to PDF and embeds interactive Plotly versions.
#'
#' @param spearman_res Named list of tibbles from `compute_spearman()`.
#' @param long_data_age Long-format NPX data frame for age variables.
#' @param long_data_other Long-format NPX data frame for other variables.
#' @param vars Character vector of metadata variable names.
#' @param out_dir Directory where plots will be saved.
#' @param n Integer; number of top hits per variable (default 10).
#' @param sign One of `NULL` (both), `"pos"`, or `"neg"`; restricts by ρ sign.
#' @param age_vars Character vector of variable names treated as age.
#'
#' @return Invisibly returns NULL.
#' @export
print_top_corr_plots <- function(spearman_res,
                                 long_data_age,
                                 long_data_other,
                                 color_col,
                                 vars,
                                 out_dir,
                                 n = 10,
                                 sign = NULL,
                                 age_vars = "Cov1") {

  # Validate sign argument
  if (!is.null(sign) && !sign %in% c("pos", "neg"))
    stop("'sign' must be NULL, 'pos', or 'neg'")

  title_prefix <- if (is.null(sign)) {
    ""
  } else if (sign == "pos") {
    "Positive "
  } else {
    "Negative "
  }
  mdcat(paste0("## Top ", n, " ", title_prefix, "Significant Correlations {.tabset}"))

  # Loop over metadata variable
  for (var in vars) {
    mdcat(paste("###", var, "{.tabset}"))
    df <- spearman_res[[var]]

    if (identical(sign, "pos")) df <- dplyr::filter(df, spearman_coefficient > 0)
    if (identical(sign, "neg")) df <- dplyr::filter(df, spearman_coefficient < 0)

    df20 <- df %>% dplyr::arrange(adj.pval) %>% dplyr::slice_head(n = n)

    if (nrow(df20) == 0) {
      mdcat(paste0("*No ", title_prefix, "significant correlations for ", var, ".*"))
      next
    }

    # Choose appropriate long-data frame
    long_df <- if (var %in% age_vars) long_data_age else long_data_other

    for (assay in unique(df20$Assay)) {

      meta_cols <- setdiff(
        colnames(long_df),                          # long_df = long_data_age or _other
        c("Assay","OlinkID","UniProt","Panel",
          "NPX","SampleID", "Panel_Lot_Nr","QC_Warning", "Assay_Warning",
          "Normalization", "Missing_Freq", "LOD", "Index", var)                    # columns we don’t want duplicated
      )


      df_plot <- long_df %>%
        dplyr::filter(Assay == assay) %>%
        # 2) build an 'extras' list-column row-wise
        dplyr::mutate(
          extras = purrr::pmap(
            dplyr::select(., all_of(meta_cols)),
            ~setNames(list(...), meta_cols)
          )
        ) %>%
        # 3) now zip exactly the four needed vectors
        dplyr::mutate(
          tooltip = purrr::pmap_chr(
            list(
              SampleID = SampleID,
              NPX       = NPX,
              var_val   = .data[[var]],
              extras    = extras
            ),
            function(SampleID, NPX, var_val, extras) {
              lines <- c(
                paste0("SampleID: ", SampleID),
                paste0("NPX: ", round(NPX, 3)),
                paste0(var, ": ", var_val),
                # flatten extras as "Name: value"
                purrr::map_chr(names(extras), ~paste0(.x, ": ", extras[[.x]]))
              )
              paste(lines, collapse = "<br>")
            }
          )
        ) %>%
        dplyr::select(-extras)  # clean up

        # Choose color mapping & legend title
        df_plot      <- df_plot %>% dplyr::mutate(color_value = as.character(color_col ))
        legend_title <- color_col

      mdcat(paste("####", assay))
      plt <- df_plot %>%
        dplyr::filter(Assay == assay) %>%
        ggplot2::ggplot(ggplot2::aes(
          x    = .data[[var]],
          y    = NPX,
          text = tooltip,
          color = color_value
        )) +
        ggplot2::geom_point() +
        ggplot2::labs(
          x = var, y = "NPX", color = legend_title,
          title = paste0(
            "rho = ", signif(df20$spearman_coefficient[df20$Assay == assay][1], 3),
            "; Adj. p = ", signif(df20$adj.pval[df20$Assay == assay][1], 3)
          )
        )

      suffix   <- ifelse(is.null(sign), "", paste0("_", sign))
      pdf_file <- file.path(out_dir, paste0(var, "_", assay, "_top20", suffix, ".pdf"))
      ggplot2::ggsave(pdf_file, plot = plt, device = "pdf", width = 7, height = 5)

      mdcat(paste0("- [Download PDF for ", assay, "](", pdf_file, ")"))
      subchunkify(plotly::ggplotly(plt, tooltip = "text"), fig.height = 5, fig.width = 7)
    }
  }
  invisible(NULL)
}


#' Print Top 10 Absolute-Magnitude Correlation Plots
#'
#' @description
#' For each metadata variable and Assay, plots NPX vs variable for the top 10
#' correlations ranked by |ρ| (absolute Spearman coefficient). Saves each plot
#' to PDF and embeds Plotly versions.
#'
#' @param spearman_res Named list of tibbles from `compute_spearman()`.
#' @param long_data_age Long-format NPX data for age variables.
#' @param long_data_other Long-format NPX data for other variables.
#' @param vars Character vector of metadata variable names.
#' @param out_dir Directory path to save plots.
#' @param age_vars Variables treated as age.
#'
#' @return Invisibly returns NULL.
#' @export
print_top10_abs <- function(spearman_res,
                            long_data_age,
                            long_data_other,
                            vars,
                            out_dir,
                            age_vars = "Age") {

  mdcat("## Top 10 Largest-Magnitude Correlations {.tabset}")

  for (var in vars) {
    mdcat(paste("###", var, "{.tabset}"))

    df10 <- spearman_res[[var]] %>%
      dplyr::mutate(abs_rho = abs(spearman_coefficient)) %>%
      dplyr::arrange(desc(abs_rho)) %>%
      dplyr::slice_head(n = 10)

    long_df <- if (var %in% age_vars) long_data_age else long_data_other

    for (assay in unique(df10$Assay)) {
      mdcat(paste("####", assay))
      plt <- long_df %>%
        dplyr::filter(Assay == assay) %>%
        ggplot2::ggplot(ggplot2::aes(
          x = .data[[var]], y = NPX, color = treatment_1lustat, text = SampleID)) +
        ggplot2::geom_point() +
        ggplot2::labs(
          x = var, y = "NPX",
          title = paste0("ρ = ", signif(df10$spearman_coefficient[df10$Assay == assay][1], 3))
        )

      pdf_file <- file.path(out_dir, paste0(var, "_", assay, "_top10abs.pdf"))
      ggplot2::ggsave(pdf_file, plot = plt, device = "pdf", width = 7, height = 5)

      mdcat(paste0("- [Download PDF for ", assay, "](", pdf_file, ")"))
      subchunkify(plotly::ggplotly(plt), fig.height = 5, fig.width = 7)
    }
  }
  invisible(NULL)
}


#' Summarize Shared Significant Correlations
#'
#' @description
#' Identifies proteins that are significant (adj.pval < alpha) across all
#' `include_factors`, optionally excluding those significant in an
#' `exclude_factor`. Computes mean adj.pval for shared proteins and writes a TSV.
#'
#' @inheritParams compute_spearman
#' @param include_factors Character vector of variable names to intersect.
#' @param exclude_factor Optional single variable to exclude.
#' @param output_path Path for the TSV export.
#' @param export Logical; whether to write TSV (default TRUE).
#'
#' @return A tibble summarizing shared-significant proteins. Invisibly returned.
#' @export
summarize_shared_significant <- function(spearman_res,
                                         alpha,
                                         include_factors,
                                         exclude_factor = NULL,
                                         output_path,
                                         export = TRUE) {
  # Validate inputs
  missing_incl <- base::setdiff(include_factors, names(spearman_res))
  if (length(missing_incl) > 0) {
    stop("include_factors not found in spearman_res: ",
         paste(missing_incl, collapse = ", "))
  }
  if (!is.null(exclude_factor) && !(exclude_factor %in% names(spearman_res))) {
    stop("exclude_factor '", exclude_factor, "' is not a name in spearman_res")
  }

  # Helper to get significant OlinkIDs per factor
  get_sig_ids <- function(tbl) {
    tbl %>%
      dplyr::filter(adj.pval < alpha) %>%
      dplyr::pull(OlinkID) %>%
      unique()
  }

  # Find IDs significant in all include_factors
  sig_lists  <- lapply(include_factors, function(f) get_sig_ids(spearman_res[[f]]))
  shared_ids <- Reduce(intersect, sig_lists)

  # Optionally exclude
  if (!is.null(exclude_factor)) {
    excl_ids  <- get_sig_ids(spearman_res[[exclude_factor]])
    shared_ids <- base::setdiff(shared_ids, excl_ids)
  }

  # Build **complete** base table of metadata from first factor
  base_tbl <- spearman_res[[ include_factors[1] ]] %>%
    dplyr::select(Assay, OlinkID, UniProt, Panel) %>%
    dplyr::distinct()

  # Compute mean adj.pval only for the shared_ids
  adj_tbl <- dplyr::bind_rows(
    lapply(include_factors, function(f) {
      spearman_res[[f]] %>%
        dplyr::filter(OlinkID %in% shared_ids) %>%
        dplyr::select(OlinkID, adj.pval)
    })
  ) %>%
    dplyr::group_by(OlinkID) %>%
    dplyr::summarise(adj.pval = mean(adj.pval), .groups = "drop")

  # Join back to the full set, fill non-shared with 1
  summary_tbl <- base_tbl %>%
    dplyr::left_join(adj_tbl, by = "OlinkID") %>%
    dplyr::mutate(
      adj.pval = tidyr::replace_na(adj.pval, 1),
      Threshold = ifelse(adj.pval < alpha, "Significant", "Non-significant")
    ) %>%
    dplyr::arrange(adj.pval)

  # Export if requested
  if (isTRUE(export)) {
    readr::write_tsv(summary_tbl, output_path)
  }
  invisible(summary_tbl)
}

#' Print Shared Significant Correlations Table
#'
#' @description
#' Prints an interactive DT table of significant shared proteins and a link to
#' the TSV file.
#'
#' @param summary_tbl Tibble from `summarize_shared_significant()`.
#' @param tsv_path Path to TSV file.
#' @param export Logical; whether to print the download link.
#'
#' @export
print_summary_table <- function(summary_tbl, tsv_path, export = TRUE) {
  if (isTRUE(export)) {
    mdcat(paste0("- [Download TSV](", tsv_path, ")"))
  }
  summary_tbl <- summary_tbl %>%
    dplyr::filter(Threshold == "Significant") %>%
    dplyr::mutate(adj.pval = sprintf("%.2e", adj.pval))

  subchunkify(
    DT::datatable(summary_tbl)
  )
  invisible(NULL)
}


#' Run Enrichment for Correlation Results
#'
#' @description
#' For each metadata variable, runs enrichment analysis with methods specified
#' in `config$enrichment$methods`. Supports `"ORA"` and `"GSEA"`. Returns a
#' nested list of results.
#'
#' @param res_list Named list of correlation tibbles.
#' @param config List containing enrichment configuration.
#'
#' @return Nested list: enrich_list[[variable]][[method]] → tibble of results.
#' @export
enrich_correlations <- function(res_list, config) {
  methods <- toupper(config$enrichment$methods %||% character())
  if (length(methods) == 0) {
    stop("config$enrichment$methods is empty; specify at least one method, e.g. ORA, GSEA")
  }
  out <- list()
  for (variable in names(res_list)) {
    tbl <- res_list[[variable]]
    out[[variable]] <- list()
    for (m in methods) {
      if (m == "ORA") {
        out[[variable]][[m]] <- enrich_ora_corr(tbl, config)
      } else if (m == "GSEA") {
        out[[variable]][[m]] <- enrich_gsea_corr(tbl, config)
      } else {
        warning("Unknown enrichment method '", m, "' – skipping for variable '", variable, "'.")
        out[[variable]][[m]] <- tibble::tibble()
      }
    }
  }
  out
}


#' ORA Enrichment for Correlation Results
#'
#' @description
#' Performs over-representation analysis (ORA) using
#' `clusterProfiler::enricher()` on significant genes from correlation results.
#' Gene sets are sourced from `msigdbr` according to `config$enrichment$ontology`.
#'
#' @param test_results Tibble of correlation results with at least column `Assay`.
#' @param config Configuration list containing:
#'   * `alpha` — numeric significance cutoff
#'   * `enrichment$organism` — species string for msigdbr
#'   * `enrichment$ontology` — ontology or MSigDb collections (e.g. `"C2"`, `"C5"`)
#'
#' @return Tibble of ORA results (columns include `ID`, `adj.pval`, `Threshold`),
#'   or an empty tibble if no results.
#' @export
enrich_ora_corr <- function(test_results, config) {
  # Install msigdbdf when running ORA, since ORA relies on msigdbdf’s MSigDB gene sets
    if (!requireNamespace("msigdbdf", quietly = TRUE)) {
    install.packages(
      "msigdbdf",
      repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org")
    )
  }



  sig <- test_results %>%
    dplyr::filter(adj.pval < config$alpha) %>%
    dplyr::distinct(Assay) %>%
    dplyr::pull(Assay)

  if (length(sig) == 0) {
    warning("No significant results to enrich – returning empty tibble.")
    return(tibble::tibble())
  }

  universe <- test_results %>%
    dplyr::distinct(Assay) %>%
    dplyr::pull(Assay)

  cols <- toupper(config$enrichment$ontology)
  if (identical(cols, "MSIGDB")) cols <- c("C2", "C5")
  msig <- lapply(cols, function(c) {
    msigdbr::msigdbr(species = config$enrichment$organism, collection = c)
  }) %>% dplyr::bind_rows()
  term2gene <- msig %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(TERM = gs_name, GENE = gene_symbol)

  # identify any genes from results table that are not present in the TERM→GENE mapping
  missing_genes <- base::setdiff(test_results$Assay, term2gene$GENE)

  if (length(missing_genes) > 0) {
    # only show the first 10 to avoid a wall of text
    to_show <- if (length(missing_genes) > 10) {
      c(head(missing_genes, 10), paste0("…(+", length(missing_genes) - 10, " more)"))
    } else {
      missing_genes
    }

    warning(
      "The following ", length(missing_genes),
      " gene", if (length(missing_genes) > 1) "s" else "",
      " were not found in your TERM2GENE mapping and will be skipped: ",
      paste(to_show, collapse = ", ")
    )
  }

  enriched <- tryCatch(
    clusterProfiler::enricher(
      gene          = sig,
      TERM2GENE     = term2gene,
      universe      = universe,
      pvalueCutoff  = config$alpha,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      warning("ORA enrichment error: ", e$message)
      NULL
    }
  )
  if (is.null(enriched)) return(tibble::tibble())
  tibble::as_tibble(enriched@result) %>%
    dplyr::rename(adj.pval = `p.adjust`) %>%
    dplyr::arrange(adj.pval) %>%
    dplyr::mutate(
      Description_short = if_else(
        stringr::str_length(Description) <= 25,
        Description,
        paste0(stringr::str_sub(Description, 1, 25), "…")
      )
    ) %>%
    mutate(
      ID = paste0(row_number(), "_", Description_short),
      Threshold = ifelse(adj.pval < config$alpha, "Significant", "Non-significant")
    ) %>%
    dplyr::select(ID, dplyr::everything()) %>%
    dplyr::select(-Description_short)
}


#' GSEA Enrichment for Correlation Results
#'
#' @description
#' Performs gene set enrichment analysis (GSEA) using `clusterProfiler::GSEA()`.
#' Ranks genes by signed -log10(p-value) from correlation results and tests for
#' enrichment across MSigDb collections.
#'
#' @param tbl Tibble with columns `Assay`, `adj.pval`, and `spearman_coefficient`.
#' @param config Configuration list containing:
#'   * `alpha` — numeric significance cutoff
#'   * `enrichment$organism` — species string for msigdbr
#'   * `enrichment$ontology` — ontology or MSigDb collections
#'
#' @return Tibble of GSEA results (columns include `ID`, `adj.pval`, `Threshold`),
#'   or an empty tibble if no results.
#' @export
enrich_gsea_corr <- function(tbl, config) {
  req <- c("Assay", "adj.pval", "spearman_coefficient")
  if (!all(req %in% names(tbl))) {
    warning("tbl missing required columns – returning empty tibble.")
    return(tibble::tibble())
  }
  pvals <- tbl$adj.pval
  pvals[pvals == 0] <- .Machine$double.xmin
  ranks <- (-log10(pvals)) * sign(tbl$spearman_coefficient)
  names(ranks) <- tbl$Assay
  ranks <- ranks[!is.na(ranks)]
  ranks <- unlist(
    lapply(split(ranks, names(ranks)), mean),
    use.names = TRUE
  )
  if (!is.numeric(ranks)) {
  stop("`ranks` must be a numeric vector, but is class: ", paste(class(ranks), collapse = ", "))
  }
  ranks <- sort(ranks, decreasing = TRUE)
  if (length(ranks) < 20) {
    warning("Fewer than 20 ranked genes—GSEA may be uninformative.")
  }
  stopifnot(!any(is.na(ranks)))
  stopifnot(length(ranks) == length(unique(names(ranks))))

  cols <- toupper(config$enrichment$ontology)
  if (identical(cols, "MSIGDB")) cols <- c("C2", "C5")
  msig <- lapply(cols, function(c) {
    msigdbr::msigdbr(species = config$enrichment$organism, collection = c)
  }) %>% dplyr::bind_rows()
  term2gene <- msig %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(TERM = gs_name, GENE = gene_symbol)

  # identify any ranked genes not present in the TERM2GENE mapping
  missing_genes <- base::setdiff(names(ranks), term2gene$GENE)

  if (length(missing_genes) > 0) {
    # only show the first 10 to avoid a wall of text
    to_show <- if (length(missing_genes) > 10) {
      c(head(missing_genes, 10), paste0("…(+", length(missing_genes) - 10, " more)"))
    } else {
      missing_genes
    }

    warning(
      "The following ", length(missing_genes),
      " gene", if (length(missing_genes) > 1) "s" else "",
      " were not found in your TERM2GENE mapping and will be skipped: ",
      paste(to_show, collapse = ", ")
    )
  }

  enriched <- tryCatch(
    clusterProfiler::GSEA(
      geneList      = ranks,
      TERM2GENE     = term2gene,
      pvalueCutoff  = config$alpha,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      warning("GSEA enrichment error: ", e$message)
      NULL
    }
  )
  if (is.null(enriched)) return(tibble::tibble())
  tibble::as_tibble(enriched@result) %>%
    dplyr::rename(adj.pval = `p.adjust`) %>%
    dplyr::arrange(adj.pval) %>%
    dplyr::mutate(
      Description_short = if_else(
        stringr::str_length(Description) <= 25,
        Description,
        paste0(stringr::str_sub(Description, 1, 25), "…")
      )
    ) %>%
    dplyr::mutate(
      ID = paste0(row_number(), "_", Description_short),
      Threshold = ifelse(adj.pval < config$alpha, "Significant", "Non-significant")
    ) %>%
    dplyr::select(ID, everything()) %>%
    dplyr::select(-Description_short)
}


#' Export Enrichment Results to Excel
#'
#' @description
#' Saves one metadata variable × method enrichment table to a single-sheet Excel
#' workbook. Uses `format_correlations_results()` for consistent formatting.
#'
#' @param enrich_tbl Tibble of enrichment results.
#' @param config Configuration list with `output_paths$correlations$results`.
#' @param variable Metadata variable name (string).
#' @param method Enrichment method name (`"ORA"` or `"GSEA"`).
#'
#' @return File path to the written workbook.
#' @export
export_enrichment_corr <- function(enrich_tbl, config, variable, method) {
  outdir <- config$output_paths$correlations$results
  wb     <- openxlsx::createWorkbook()
  df_fmt <- enrich_tbl %>% format_correlations_results(drop_enrich_cols = FALSE, scientific_notation = FALSE)
  openxlsx::addWorksheet(wb, method)
  openxlsx::writeData(wb, sheet = method, x = df_fmt, withFilter = TRUE)
  fname <- file.path(outdir, paste0(variable, "_", method, "_enrichment.xlsx"))
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
  fname
}


#' Print and Export One Enrichment Table
#'
#' @description
#' Writes one enrichment results workbook to disk, prints a download link, and
#' displays a DT datatable of significant terms (adj.pval < alpha).
#'
#' @param enrich_tbl Tibble of enrichment results.
#' @param config Configuration list with `alpha` and output paths.
#' @param variable Metadata variable name.
#' @param method Enrichment method name.
#'
#' @export
print_enrichment_corr <- function(enrich_tbl, config, variable, method) {
  wb_path <- export_enrichment_corr(enrich_tbl, config, variable, method)
  mdcat(paste0(
    "- [Download ", method, " enrichment results for ", variable, " Excel](",
    wb_path, ")"
  ))
  subchunkify(
    DT::datatable(
      enrich_tbl %>%
        format_correlations_results(drop_enrich_cols = TRUE, scientific_notation = TRUE) %>%
        dplyr::filter(Threshold == "Significant")
    )
  )
  invisible(NULL)
}


#' Print All Enrichment Results
#'
#' @description
#' Iterates over all metadata variables × methods in an enrichment result list,
#' printing either a datatable of results or a “no results” message.
#'
#' @param enrich_list Named list: `enrich_list[[variable]][[method]] → tibble`.
#' @param config Configuration list with `enrichment$methods`, `alpha`, and output paths.
#'
#' @export
print_enrichment_results_corr <- function(enrich_list, config) {
  vars    <- names(enrich_list)
  methods <- toupper(config$enrichment$methods %||% character())
  if (length(vars) == 0 || length(methods) == 0) {
    mdcat("No metadata variables or enrichment methods specified; check your config.")
    return(invisible(NULL))
  }
  for (var in vars) {
    mdcat(paste0("### ", var, "{.tabset}"))
    for (m in methods) {
      mdcat(paste0("#### ", m))
      tbl <- enrich_list[[var]][[m]]
      if (is.null(tbl) || !"adj.pval" %in% names(tbl) || nrow(tbl) == 0) {
        mdcat(paste0("*No ", m, " enrichment results for ", var, ".*"))
      } else {
        print_enrichment_corr(tbl, config, variable = var, method = m)
      }
    }
  }
  invisible(NULL)
}


#' Dotplot for GSEA Correlation Results
#'
#' @description
#' Creates a dotplot of top 20 GSEA results, colored by -log10(adj.pval) and
#' sized by set size. Saves as PDF and returns interactive Plotly version.
#'
#' @param df Tibble of GSEA results.
#' @param config Configuration list with `plots_dir`.
#' @param tag Title/filename prefix.
#'
#' @return List with elements `p_plotly` (Plotly object) and `path` (PDF path).
#' @export
make_dotplot_gsea_corr <- function(df, config, tag) {
  df_top <- df %>%
    dplyr::slice_head(n = 20)

  p_gg <- ggplot2::ggplot(df_top, ggplot2::aes(
    x = NES, y = stats::reorder(ID, NES),
    size = setSize, color = -log10(adj.pval)
  )) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste(tag, "GSEA Dotplot"),
      x = "Normalized Enrichment Score (NES)",
      y = "Gene-set ID",
      color = "-log10(adj.pval)",
      size = "Size"
    )
  path <- file.path(config$plots_dir, paste0(tag, "_GSEA_dotplot.pdf"))
  ggplot2::ggsave(path, p_gg, device = "pdf", width = 7, height = 5)
  p_plotly <- plotly::ggplotly(
    p_gg + ggplot2::aes(
      text = paste0(
        "Term: ", Description,
        "<br>NES: ", round(NES,3),
        "<br>adj.pval: ", sprintf("%.2e", adj.pval),
        "<br>Size: ", setSize
      )
    ), tooltip = "text"
  )
  list(p_plotly = p_plotly, path = path)
}


#' Dotplot for ORA Correlation Results
#'
#' @description
#' Creates a dotplot of top 20 ORA results, colored by -log10(adj.pval) and
#' sized by hit count. Saves as PDF and returns interactive Plotly version.
#'
#' @param df Tibble of ORA results.
#' @param config Configuration list with `plots_dir`.
#' @param tag Title/filename prefix.
#'
#' @return List with elements `p_plotly` and `path`.
#' @export
make_dotplot_ora_corr <- function(df, config, tag) {
  df_top <- df %>%
    dplyr::slice_head(n = 20)

  p_gg <- ggplot2::ggplot(df_top, ggplot2::aes(
    x = RichFactor, y = stats::reorder(ID, RichFactor),
    size = Count, color = -log10(adj.pval)
  )) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste(tag, "ORA Dotplot"),
      x = "RichFactor", y = "Gene-set ID",
      color = "-log10(adj.pval)", size = "Count"
    )
  path <- file.path(config$plots_dir, paste0(tag, "_ORA_dotplot.pdf"))
  ggplot2::ggsave(path, p_gg, device = "pdf", width = 7, height = 5)
  p_plotly <- plotly::ggplotly(
    p_gg + ggplot2::aes(
      text = paste0(
        "Term: ", Description,
        "<br>RichFactor: ", round(RichFactor,3),
        "<br>adj.pval: ", sprintf("%.2e", adj.pval),
        "<br>Count: ", Count
      )
    ), tooltip = "text"
  )
  list(p_plotly = p_plotly, path = path)
}


#' Dotplot Wrapper for Correlation Enrichment
#'
#' @description
#' Wrapper to call the appropriate dotplot function (`make_dotplot_ora_corr` or
#' `make_dotplot_gsea_corr`) based on `method`.
#'
#' @param enrich_tbl Tibble of enrichment results.
#' @param config Configuration list with output paths.
#' @param method Method string (`"ORA"` or `"GSEA"`).
#' @param tag Title/filename prefix.
#'
#' @return List with Plotly object and PDF path.
#' @export
dotplot_corr <- function(enrich_tbl, config, method, tag) {
  method <- match.arg(toupper(method), c("ORA","GSEA"))
  if (method == "ORA") {
    make_dotplot_ora_corr(enrich_tbl, config, tag = tag)
  } else {
    make_dotplot_gsea_corr(enrich_tbl, config, tag = tag)
  }
}


#' Print All Enrichment Dotplots
#'
#' @description
#' Creates a nested tabset of dotplots for each metadata variable and method.
#' Prints download links, embeds interactive Plotly figures, or prints messages
#' when insufficient results are available.
#'
#' @param enrich_list Named list of enrichment results.
#' @param config Configuration list with enrichment methods, alpha, and output paths.
#'
#' @export
dotplot_corr_list <- function(enrich_list, config) {
  vars    <- names(enrich_list)
  methods <- toupper(config$enrichment$methods %||% character())

  # Define default minimum term cutoffs
  default_min <- list(ORA = 10, GSEA = 5)
  min_terms   <- config$enrichment$min_terms %||% default_min

  if (length(vars) == 0 || length(methods) == 0) {
    mdcat("No metadata variables or enrichment methods specified; check your config.")
    return(invisible(NULL))
  }

  for (var in vars) {
    mdcat(paste0("### ", var, " {.tabset}"))
    for (m in methods) {
      mdcat(paste0("#### ", m))
      tbl <- enrich_list[[var]][[m]]

      # No table or no rows
      if (is.null(tbl) || !"adj.pval" %in% names(tbl) || nrow(tbl) == 0) {
        mdcat(sprintf("*No %s enrichment results to plot for %s.*", m, var))
        next
      }

      # Check cutoff for ORA/GSEA
      cutoff <- min_terms[[m]] %||% default_min[[m]]
      num_sig <- sum(tbl$adj.pval < config$alpha, na.rm = TRUE)
      if (m %in% names(default_min) && num_sig < cutoff) {
        mdcat(sprintf("*Insufficient (< %d) %s enrichment results to plot for %s.*", cutoff, m, var))
        next
      }

      # Plot: PDF link + embed
      tag <- paste(var, m, sep = "_")
      plt_info <- dotplot_corr(tbl, config, method = m, tag = tag)
      mdcat(sprintf("- [Download %s %s Dotplot PDF](%s)", m, var, plt_info$path))
      subchunkify(plt_info$p_plotly)
    }
  }

  invisible(NULL)
}


#' Report Significant Correlation Counts
#'
#' @description
#' Counts the number of significant proteins (`Threshold == "Significant"`)
#' for each metadata variable.
#'
#' @param spearman_res Named list of correlation result tibbles.
#'
#' @return Tibble with columns `variable` and `n_significant`.
#' @export
report_num_sig_corr <- function(spearman_res) {
  purrr::map_dfr(
    names(spearman_res),
    function(variable) {
      df    <- spearman_res[[variable]]
      n_sig <- sum(df$Threshold == "Significant", na.rm = TRUE)
      tibble::tibble(
        variable      = variable,
        n_significant = n_sig
      )
    }
  )
}


#' Get Number of Significant Correlations for One Variable
#'
#' @description
#' Retrieves the number of significant proteins for a single metadata variable
#' from `report_num_sig_corr()`.
#'
#' @param spearman_res Named list of correlation result tibbles.
#' @param variable Metadata variable name.
#'
#' @return Integer count of significant proteins.
#' @export
get_num_sig_corr <- function(spearman_res, variable) {
  report_num_sig_corr(spearman_res) %>%
    dplyr::filter(.data$variable == {{variable}}) %>%
    dplyr::pull(n_significant)
}


#' Plot Spearman Coefficient Distributions by Category
#'
#' @description
#' For selected metadata categories, draws violin plots of Spearman coefficients
#' across proteins. Significant results are highlighted within each violin. Saves
#' PDF and prints download link.
#'
#' @param spearman_res Named list of correlation result tibbles.
#' @param outdir Directory path to save plots.
#' @param categories Character vector of variable names to include.
#' @param y_lim Numeric vector of y-axis limits.
#' @param outfile File name for the PDF.
#' @param width,height PDF dimensions in inches.
#' @param order_by How to order categories: `"given"`, `"median_abs"`, `"pct_sig"`.
#' @param show_counts Logical; annotate counts above violins.
#' @param seed RNG seed for reproducibility.
#'
#' @return Invisibly returns the ggplot object with attribute `"pdf_path"`.
#' @export
plot_spearman_violins <- function(spearman_res,
                                  outdir,
                                  categories,
                                  y_lim = c(-1, 1),
                                  outfile = "spearman_violins.pdf",
                                  width = 6, height = 4.5,
                                  order_by = c("given", "median_abs", "pct_sig"),
                                  show_counts = TRUE,
                                  seed = 1L) {
  order_by <- base::match.arg(order_by)

  # validate/save path
  if (base::is.null(outdir) || !base::nzchar(outdir)) base::stop("config$output_path$correlations$plots must be a valid directory path (string).")
  if (!base::dir.exists(outdir)) base::dir.create(outdir, recursive = TRUE)
  pdf_path <- base::file.path(outdir, outfile)

  # collect and tag categories
  cats_present <- base::intersect(categories, base::names(spearman_res))
  if (!base::length(cats_present)) base::stop("None of the requested categories found in 'spearman_res'.")
  df <- dplyr::bind_rows(lapply(cats_present, function(nm) {
    x <- spearman_res[[nm]]
    x$Category <- nm
    x
  }))

  # significance flag from "Threshold"
  thr <- df$Threshold
  thr[base::is.na(thr)] <- ""
  df$significant <- base::startsWith(base::tolower(base::as.character(thr)), "signific")

  # summaries for ordering & labels
  summ <- df |>
    dplyr::group_by(.data$Category) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_sig = base::sum(.data$significant, na.rm = TRUE),
      n_pos_sig = base::sum(.data$significant & .data$spearman_coefficient > 0, na.rm = TRUE),
      n_neg_sig = base::sum(.data$significant & .data$spearman_coefficient < 0, na.rm = TRUE),
      med_abs = stats::median(base::abs(.data$spearman_coefficient), na.rm = TRUE),
      pct_sig = n_sig / n,
      .groups = "drop"
    )

  # x order
  if (order_by == "median_abs") {
    x_levels <- summ$Category[base::order(summ$med_abs, decreasing = TRUE)]
  } else if (order_by == "pct_sig") {
    x_levels <- summ$Category[base::order(summ$pct_sig, decreasing = TRUE)]
  } else {
    x_levels <- categories[categories %in% cats_present]
  }
  df$Category <- base::factor(df$Category, levels = x_levels)
  summ$Category <- base::factor(summ$Category, levels = x_levels)

  # thresholds where "significant proteins start appearing" per Category
  thresh <- df |>
    dplyr::group_by(.data$Category) |>
    dplyr::summarise(
      t_pos = if (base::any(.data$significant & .data$spearman_coefficient > 0, na.rm = TRUE))
                 base::min(.data$spearman_coefficient[.data$significant & .data$spearman_coefficient > 0], na.rm = TRUE)
               else NA_real_,
      t_neg = if (base::any(.data$significant & .data$spearman_coefficient < 0, na.rm = TRUE))
                 base::max(.data$spearman_coefficient[.data$significant & .data$spearman_coefficient < 0], na.rm = TRUE)
               else NA_real_,
      .groups = "drop"
    )

  # y-range per category (respect y_lim if provided)
  if (!base::is.null(y_lim) && base::all(!base::is.na(y_lim))) {
    yrng <- dplyr::tibble(Category = base::factor(x_levels, levels = x_levels),
                          y_min = y_lim[1], y_max = y_lim[2])
  } else {
    yrng <- df |>
      dplyr::group_by(.data$Category) |>
      dplyr::summarise(y_min = base::min(.data$spearman_coefficient, na.rm = TRUE),
                       y_max = base::max(.data$spearman_coefficient, na.rm = TRUE),
                       .groups = "drop")
  }

  # rectangles to color (masked to violin shape later)
  rect_top <- dplyr::left_join(thresh, yrng, by = "Category") |>
    dplyr::filter(!base::is.na(.data$t_pos)) |>
    dplyr::mutate(xmin = base::as.numeric(.data$Category) - 0.45,
                  xmax = base::as.numeric(.data$Category) + 0.45,
                  ymin = .data$t_pos, ymax = .data$y_max) |>
    dplyr::select(.data$xmin, .data$xmax, .data$ymin, .data$ymax, .data$Category)

  rect_bot <- dplyr::left_join(thresh, yrng, by = "Category") |>
    dplyr::filter(!base::is.na(.data$t_neg)) |>
    dplyr::mutate(xmin = base::as.numeric(.data$Category) - 0.45,
                  xmax = base::as.numeric(.data$Category) + 0.45,
                  ymin = .data$y_min, ymax = .data$t_neg) |>
    dplyr::select(.data$xmin, .data$xmax, .data$ymin, .data$ymax, .data$Category)

  rects <- dplyr::bind_rows(rect_top, rect_bot)

  # label y
  y_top <- if (!base::is.null(y_lim) && base::all(!base::is.na(y_lim))) y_lim[2] else base::max(df$spearman_coefficient, na.rm = TRUE)
  y_offset <- 0.01 * base::diff(base::range(c(-1, 1)))
  summ$y_label <- y_top - y_offset
  summ$label <- base::paste0("n=", summ$n_sig, " (", summ$n_pos_sig, "+/", summ$n_neg_sig, "-)")

  base::set.seed(seed)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Category, y = .data$spearman_coefficient)) +
    # use violin as a MASK reference
    ggfx::as_reference(
      ggplot2::geom_violin(fill = "grey85", color = NA, width = 0.9, trim = TRUE),
      id = "violin_mask"
    ) +
    # color only the parts INSIDE the violin using the mask
    ggfx::with_mask(
      ggplot2::geom_rect(data = rects,
                         ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                      ymin = .data$ymin, ymax = .data$ymax),
                         inherit.aes = FALSE, fill = "#B22222", alpha = 0.6),
      mask = "violin_mask"
    ) +
    # draw the violin outline on top
    ggplot2::geom_violin(fill = NA, color = "black", width = 0.9, trim = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    (if (!base::is.null(y_lim) && base::all(!base::is.na(y_lim))) ggplot2::coord_cartesian(ylim = y_lim) else ggplot2::guides()) +
    ggplot2::labs(x = NULL, y = "spearman_coefficient") +
    ggplot2::theme(legend.position = "none")

  if (base::isTRUE(show_counts)) {
    p <- p + ggplot2::geom_text(
      data = summ,
      ggplot2::aes(x = .data$Category, y = .data$y_label, label = .data$label),
      size = 3, vjust = 1
    )
  }

  ggplot2::ggsave(filename = pdf_path, plot = p, device = "pdf", width = width, height = height)

  link_text <- base::sprintf("[%s](%s)\n\n", base::basename(pdf_path), pdf_path)
  mdcat(link_text)
  base::print(p)

  base::attr(p, "pdf_path") <- pdf_path
  base::invisible(p)
}
