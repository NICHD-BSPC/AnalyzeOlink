#' QQ plots of raw NPX by contrast, arm, and protein
#'
#' Builds a three-level tabset in the knitted report (**Contrast → Arm → Protein**),
#' drawing QQ plots of *raw NPX* per arm and protein, and writing a PDF per plot
#' to `config$output_paths$differential_expression$qc`.
#'
#' @param contrast_list List whose first element describes the contrast and
#'   includes fields `name`, `se_object_name`, `sampletable_column`,
#'   `treatment`, and `control`.
#' @param se_list Named list of `SummarizedExperiment` objects indexed by
#'   `se_object_name`.
#' @param config List with `output_paths$differential_expression$qc` (writable directory).
#' @param proteins Character vector of assay IDs to include (e.g., `c("CTSL")`).
#' @return `invisible(NULL)`; called for side effects (plots, PDFs, markdown).
#' @export
plot_qq_npx <- function(contrast_list, se_list, config, proteins) {
  # We are calling this function on one contrast at a time using
  # a contrast specific protein set
  con <- contrast_list[[1]]
  contrast_name <- con$name
  se_obj        <- se_list[[con$se_object_name]]
  group_col     <- con$sampletable_column
  arms          <- c(con$treatment, con$control)

  long_df <- S4Vectors::metadata(se_obj)$long_data %>%
    dplyr::filter(Assay %in% proteins)

  for (arm in arms) {
    mdcat(paste0("##### ", arm, " {.tabset}"))            # level-2 tabset
    arm_df <- long_df[long_df[[group_col]] == arm, , drop = FALSE]

    for (prot in unique(arm_df$Assay)) {
      mdcat(paste0("###### ", prot))                      # level-3 tab

      prot_df <- arm_df[arm_df$Assay == prot, , drop = FALSE]
      if (nrow(prot_df) < 3) {
        mdcat("*Not enough observations to draw QQ plot*")
        next
      }

      qqp <- ggplot2::ggplot(prot_df, ggplot2::aes(sample = NPX)) +
        ggplot2::stat_qq() +
        ggplot2::stat_qq_line() +
        ggplot2::labs(
          title = paste(contrast_name, arm, prot, sep = "  –  "),
          x = "Theoretical Quantiles",
          y = "Sample Quantiles"
        )

      pdf_path <- file.path(
        config$output_paths$differential_expression$qc,
        paste0(contrast_name, "_", arm, "_", prot, "_qq.pdf")
      )
      ggplot2::ggsave(pdf_path, qqp, width = 6, height = 5)
      mdcat("- [Download PDF](", pdf_path, ")")

      subchunkify(qqp)
    }
  }
  invisible(NULL)
}


#' Histogram of raw NPX per contrast
#'
#' For each contrast, extracts the NPX assay matrix from the referenced
#' `SummarizedExperiment`, flattens to numeric (ignoring `NA` and literal
#' `"NULL"` strings), draws a histogram with a Normal overlay, writes a PDF,
#' and previews it in the report.
#'
#' @param res_list Named list of contrast results; each element must include
#'   `params$se_object_name` identifying the NPX source.
#' @param se_list Named list of `SummarizedExperiment` objects containing `"NPX"`.
#' @param outdir String; output directory for PDFs (must exist).
#' @param breaks Integer; number of histogram bins (default `30`).
#' @return `invisible(NULL)`; side effects are plots, PDFs, and markdown.
#' @export
npx_hist <- function(res_list, se_list, outdir, breaks = 30) {
  for (contrast_name in names(res_list)) {

    # Locate the SummarizedExperiment holding the raw NPX measurements
    se_name <- res_list[[contrast_name]]$params$se_object_name
    if (is.null(se_name))
      stop(sprintf(
        "Could not find `se_object_name` for contrast '%s'", contrast_name))

    se_obj <- se_list[[se_name]]
    if (!inherits(se_obj, "SummarizedExperiment"))
      stop(sprintf("Object '%s' is not a SummarizedExperiment", se_name))

    # Pull out ALL raw NPX values and flatten to a numeric vector
    npx <- SummarizedExperiment::assay(se_obj, "NPX")
    npx <- as.character(npx)
    npx <- npx[!is.na(npx) & npx != "NULL"]
    npx <- as.numeric(npx)

    # Summary label
    lab_stats <- sprintf(
      "n = %d,  mean = %.2f  (sd = %.2f)",
      length(npx), mean(npx), stats::sd(npx))

    # Set up paths
    pdf_path <- file.path(
      outdir,
      paste0(contrast_name, "_NPX_hist.pdf")
    )

    # PDF version
    grDevices::pdf(pdf_path)

    ## (a) Build *silent* histogram to capture bin info
    h <- graphics::hist(npx, breaks = breaks, plot = FALSE)

    ## (b) Draw the bars
    graphics::plot(
      h, col = "black", border = "white",
      main = paste("Raw NPX values:", contrast_name),
      xlab = "NPX"
    )

    ## (c) Overlay the theoretical Normal (same μ and σ as your data)
    m <- mean(npx);  s <- stats::sd(npx)
    xfit <- seq(min(npx), max(npx), length.out = 200)
    yfit <- stats::dnorm(xfit, m, s) * diff(h$breaks[1:2]) * length(npx)
    graphics::lines(xfit, yfit, lwd = 2, col = "red")

    graphics::mtext(lab_stats, side = 3, line = 0.5, cex = 0.8)
    grDevices::dev.off()

    # Markdown preview
    mdcat(paste0("#### ", contrast_name), "\n")
    mdcat(lab_stats, "\n")
    mdcat("- [Download PDF](", pdf_path, ")\n\n")

    ## 7. Inline preview in the Rmd report
    subchunkify({
      h <- graphics::hist(npx, breaks = breaks, plot = FALSE)
      graphics::plot(
        h, col = "black", border = "white",
        main = paste("Raw NPX values:", contrast_name),
        xlab = "NPX"
      )
      m <- mean(npx);  s <- stats::sd(npx)
      xfit <- seq(min(npx), max(npx), length.out = 200)
      yfit <- stats::dnorm(xfit, m, s) * diff(h$breaks[1:2]) * length(npx)
      graphics::lines(xfit, yfit, lwd = 2, col = "red")
      graphics::mtext(lab_stats, side = 3, line = 0.5, cex = 0.8)
    })
  }

  invisible(NULL)
}


#' Normality-test p-value histograms by contrast
#'
#' Filters each contrast to its main effect, checks required columns, plots the
#' histogram of test p-values, writes a PDF, and prints a concise summary of
#' (i) tests with FDR < `alpha` and (ii) overlap with ANOVA-significant effects.
#'
#' @param res_list Named list of contrast results; each element contains
#'   `$anova` and `$params$column` (the main effect name).
#' @param test One of `"SW"`, `"AD"`, or `"FT"` (case-insensitive).
#' @param alpha Numeric in (0,1); FDR threshold used in the summary.
#' @param outdir String; output directory for PDFs (must exist).
#' @param breaks Integer; number of histogram bins (default `10`).
#' @return A tibble with columns `contrast`, `test_type`,
#'   `n_test_significant`, and `n_effect_and_test_significant`.
#' @export
norm_hist <- function(res_list, test, alpha, outdir, breaks = 10) {
  test <- toupper(test)
  if (test == "SW") {
    pval_col   <- "SW_pval"
    padj_col   <- "SW_padj"
    test_label <- "Shapiro–Wilk"
  } else if (test == "AD") {
    pval_col   <- "AD_pval"
    padj_col   <- "AD_padj"
    test_label <- "Anderson–Darling"
  } else if (test == "FT") {
    pval_col   <- "FT_pval"
    padj_col   <- "FT_padj"
    test_label <- "F-test"
  } else {
    stop('`test` must be "SW", "AD", or "FT"')
  }

  summary_list <- vector("list", length(res_list))
  i <- 1

  for (contrast_name in names(res_list)) {
    # 1) grab & filter to main effect
    tbl       <- res_list[[contrast_name]]$anova
    main_eff  <- res_list[[contrast_name]]$params$column
    tbl_main  <- tbl %>% dplyr::filter(Effect == main_eff)

    # 2) check required columns
    required <- c(pval_col, padj_col, "Threshold")
    missing  <- setdiff(required, names(tbl_main))
    if (length(missing)) {
      stop(sprintf(
        "Contrast '%s' is missing column(s): %s",
        contrast_name, paste(missing, collapse = ", ")
      ))
    }

    # 3) compute counts
    pvals    <- tbl_main[[pval_col]]
    padjs    <- tbl_main[[padj_col]]
    thres    <- tbl_main$Threshold

    n_total                  <- length(pvals)
    n_test_significant       <- sum(padjs < alpha, na.rm = TRUE)
    n_effect_and_test_sig    <- sum((padjs < alpha) & (thres == "Significant"), na.rm = TRUE)
    total_anova_sig          <- sum(thres == "Significant", na.rm = TRUE)

    pct_test    <- n_test_significant    / n_total       * 100
    pct_inter   <- if (total_anova_sig>0)
                     n_effect_and_test_sig / total_anova_sig * 100
                   else
                     NA_real_

    # 4) labels
    label1 <- sprintf(
      "Test FDR < %.3g: %d (%.1f%% of %d)",
      alpha,
      n_test_significant,
      pct_test,
      n_total
    )
    label2 <- sprintf(
      "ANOVA sig. & non-normal: %d (%.1f%% of %d)",
      n_effect_and_test_sig,
      pct_inter,
      total_anova_sig
    )

    # 5) print markdown + histogram
    mdcat(paste0("#### ", contrast_name), "\n")
    mdcat(label1, "\n")
    mdcat(label2, "\n\n")

    pdf_path <- file.path(
      outdir,
      paste0(contrast_name, "_", test, "_pval_hist.pdf")
    )
    grDevices::pdf(pdf_path)
    graphics::hist(
      pvals, breaks = breaks,
      main = paste(test_label, "p-values:", contrast_name),
      xlab = paste(test_label, "p-value"),
      col  = "gray", border = "white"
    )
    graphics::mtext(label1, side = 3, line = 0.5, cex = 0.8)
    graphics::mtext(label2, side = 3, line = -0.5, cex = 0.8)
    grDevices::dev.off()
    mdcat("- [Download PDF](", pdf_path, ")\n")

    subchunkify({
      graphics::hist(
        pvals, breaks = breaks,
        main = paste(test_label, "p-values:", contrast_name),
        xlab = paste(test_label, "p-value"),
        col  = "gray", border = "white"
      )
      graphics::mtext(label1, side = 3, line = 0.5, cex = 0.8)
      graphics::mtext(label2, side = 3, line = -0.5, cex = 0.8)
    })

    # 6) collect summary
    summary_list[[i]] <- tibble::tibble(
      contrast                      = contrast_name,
      test_type                     = test,
      n_test_significant            = n_test_significant,
      n_effect_and_test_significant = n_effect_and_test_sig
    )
    i <- i + 1
  }

  dplyr::bind_rows(summary_list)
}


#' QQ plots of residuals by p-value bin (per protein)
#'
#' Builds a tabset per contrast, subdivides proteins into p-value bins (by
#' `SW_pval`), and for up to `num_proteins` per bin draws QQ plots of model
#' residuals. Optional labels can be added via `label_map`.
#'
#' @param res_list Named list from `run_tests()`; each element has `$anova`
#'   with columns `Assay`, `Panel`, `SW_pval`, and list-column `resids`.
#' @param outdir String; output directory for PDFs (must exist).
#' @param label_map Optional list; `label_map[[contrast]]` is a named list whose
#'   names are `SampleID`s and whose values are character vectors of
#'   `"Assay_Panel"` to label in the QQ plot.
#' @param num_proteins Integer; max proteins per bin (default `10`).
#' @param bins Numeric vector of cut points in (0,1] (default `c(0,0.1,0.25,0.50,0.75,1)`).
#' @param seed Integer; RNG seed for reproducibility.
#' @return `invisible(NULL)`; side effects are plots, PDFs, and markdown.
#' @export
plot_qq_protein <- function(res_list,
                            outdir,
                            label_map    = NULL,
                            num_proteins = 10,
                            bins         = c(0, 0.1, 0.25, 0.50, 0.75, 1),
                            seed         = 1) {

  # Checks
  stopifnot(
    dir.exists(outdir),
    is.list(res_list),
    is.list(label_map) | is.null(label_map)
  )

  for (contrast_name in names(res_list)) {

    # Gather ANOVA results
    term <- res_list[[contrast_name]]$params$column
    anova_tbl <- res_list[[contrast_name]]$anova %>%
      dplyr::filter(Effect == term)

    if (!all(c("Assay", "Panel", "SW_pval", "resids") %in% names(anova_tbl))) {
      warning("Contrast ", contrast_name, " lacks expected columns; skipped.")
      next
    }

    mdcat(paste0("#### ", contrast_name, " {.tabset}"))

    # Build p-value bins
    cuts   <- bins
    labels <- paste0(
      "Bin ", seq_len(length(cuts) - 1),
      " (", head(cuts, -1), ", ", tail(cuts, -1), "]"
    )

    bin_fac <- cut(anova_tbl$SW_pval,
                   breaks = cuts,
                   include.lowest = TRUE,
                   labels = labels)

    bin_tabs <- lapply(labels, function(lb) {
      anova_tbl[bin_fac == lb, , drop = FALSE] %>%
        dplyr::arrange(SW_pval) %>%
        dplyr::slice_head(n = num_proteins)
    })
    names(bin_tabs) <- labels

    # Iterate over bins
    for (bin_name in names(bin_tabs)) {
      mdcat(paste0("##### ", bin_name, " {.tabset}"))
      bin_df <- bin_tabs[[bin_name]]

      if (nrow(bin_df) == 0) {
        mdcat("*No proteins fall into this bin*")
        next
      }

      # Iterate over selected proteins
      for (i in seq_len(nrow(bin_df))) {

        prot        <- bin_df$Assay[i]
        panel       <- bin_df$Panel[i]
        prot_panel  <- paste0(prot, "_", panel)

        mdcat(paste0("###### ", prot, " (", panel, ")"))

        resid_vec <- unlist(bin_df$resids[[i]])
        if (length(resid_vec) < 3) {
          mdcat("*Not enough residuals to draw a QQ-plot*")
          next
        }

        qq_stats <- stats::qqnorm(resid_vec, plot.it = FALSE)
        plot_df  <- tibble::tibble(
          SampleID    = names(resid_vec),
          resid       = resid_vec,
          theor_quant = qq_stats$x,
          samp_quant  = qq_stats$y
        )

        # Pick labels: which SampleIDs contain prot_panel in their vector?
        lbls <- if (is.list(label_map)) {
          contr_map <- label_map[[contrast_name]]
          if (is.null(contr_map)) character(0) else
            names(contr_map)[vapply(
              contr_map,
              function(v) prot_panel %in% v,
              logical(1)
            )]
        } else {
          label_map
        }

        # Build QQ plot
        qqp <- ggplot2::ggplot(plot_df,
                               ggplot2::aes(sample = resid)) +
          ggplot2::stat_qq() +
          ggplot2::stat_qq_line() +
          ggplot2::labs(
            title = paste(contrast_name, "\n", prot_panel),
            x = "Theoretical Quantiles",
            y = "Sample Quantiles"
          )

        if (length(lbls) > 0) {
          qqp <- qqp +
            ggrepel::geom_label_repel(
              data  = dplyr::filter(plot_df, SampleID %in% lbls),
              ggplot2::aes(x = theor_quant, y = samp_quant, label = SampleID),
              inherit.aes       = FALSE,
              color             = "black",
              fill              = "white",
              segment.color     = "black",
              treatment_1.segment.length = 0,
              label.size        = 0.25,
              box.padding       = 0.7,
              point.padding     = 0.25,
              max.overlaps      = nrow(plot_df),
              show.legend       = FALSE
            )
        }

        # Save PDF and embed
        pdf_path <- file.path(
          outdir,
          paste0(contrast_name, "_", prot_panel, "_qq_residuals.pdf")
        )
        ggplot2::ggsave(filename = pdf_path, plot = qqp)

        mdcat("- [Download PDF](", pdf_path, ")")
        subchunkify(qqp)
      }
    }
  }

  invisible(NULL)
}


#' Residual outlier frequency tables (per contrast)
#'
#' Selects proteins whose Shapiro–Wilk p-values fall in `pval_bin`, flags
#' residual outliers per protein via a Tukey IQR rule (`k` × IQR), and returns
#' a nested list of outliers. If `print_plot = TRUE`, also writes a per-contrast
#' frequency table to Excel plus an optional bar plot, and returns the list of
#' frequency tables.
#'
#' @param res_list Named list from `run_tests()`; each element has `$anova`
#'   with `Assay`, `Panel`, `SW_pval`, and list-column `resids` (names are `SampleID`).
#' @param outdir String; output directory for workbooks/plots (must exist).
#' @param pval_bin Numeric length-2 vector giving the open/closed interval
#'   `(low, high]` for selecting proteins (default `c(0, 0.1)`).
#' @param k Numeric selectivity constant for Tukey fences (default `1.5`).
#' @param print_plot Logical; if `TRUE`, write Excel + plots and return frequency tables.
#' @return If `print_plot = FALSE`, returns (invisibly) `outlier_list` where
#'   `outlier_list[[contrast]][[SampleID]]` is a character vector of `"Assay_Panel"`.
#'   If `print_plot = TRUE`, returns (invisibly) a named list of frequency tables.
#' @export
write_outlier_frequency_table <- function(res_list,
                                          outdir,
                                          pval_bin   = c(0, 0.1),
                                          k          = 1.5,
                                          print_plot = FALSE
) {

  # Checks
  stopifnot(
    dir.exists(outdir),
    is.list(res_list),
    length(pval_bin) == 2,
    pval_bin[1] < pval_bin[2]
  )

  workbook_paths <- character()
  outlier_list   <- list()   # NESTED: [[contrast]][[SampleID]][[Protein_Panel]]
  freq_tbl_list  <- list()   # NESTED: [[contrast]]

  for (contrast in names(res_list)) {

    # Gather p-values and residuals
    term <- res_list[[contrast]]$params$column
    anova_tbl <- res_list[[contrast]]$anova %>%
      dplyr::filter(Effect == term) %>%
      dplyr::arrange(SW_pval)

    if (!all(c("Assay", "Panel", "SW_pval", "resids") %in% names(anova_tbl)))
      next

    sel <- anova_tbl$SW_pval  > pval_bin[1] &
           anova_tbl$SW_pval <= pval_bin[2]

    # Containers
    sample_outliers <- list()   # <- build this new structure
    all_outlier_ids <- character(0)

    for (idx in which(sel)) {

      prot_panel <- paste0(anova_tbl$Assay[idx], "_", anova_tbl$Panel[idx])
      resid_vec  <- unlist(anova_tbl$resids[[idx]])
      if (!length(resid_vec)) next

      # IQR outlier rule
      q1  <- stats::quantile(resid_vec, 0.25, na.rm = TRUE)
      q3  <- stats::quantile(resid_vec, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      low <- q1 - k * iqr
      upp <- q3 + k * iqr
      is_out <- resid_vec < low | resid_vec > upp

      ids <- names(resid_vec)[is_out]

      # Fill the nested list: [[SampleID]][[Protein_Panel]]
      for (sid in ids) {
        sample_outliers[[sid]][[prot_panel]] <- TRUE
      }

      all_outlier_ids <- c(all_outlier_ids, ids)
    }

    outlier_list[[contrast]] <- lapply(sample_outliers, names)

    if (print_plot) {

      # Build frequency table for the bar-plot
      freq_tbl <- if (length(all_outlier_ids)) {
        out <- as.data.frame(table(all_outlier_ids),
                             stringsAsFactors = FALSE,
                             responseName     = "Frequency")
        names(out)[1] <- "SampleID"
        out[order(-out$Frequency, out$SampleID), ]
      } else {
        data.frame(SampleID = character(), Frequency = integer())
      }

      freq_tbl_list[[contrast]] <- freq_tbl

      # Save workbook and add to report
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "outlier_frequency")
      openxlsx::writeData(wb, 1, freq_tbl, withFilter = TRUE)

      file_name <- paste0(contrast, "_residual_outliers.xlsx")
      xlsx_path <- file.path(outdir, file_name)
      openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
      workbook_paths[contrast] <- xlsx_path

      mdcat(paste0("#### ", contrast))
      mdcat(paste0("- [Download frequency table](", xlsx_path, ")"))

      if (nrow(freq_tbl)) {
        # PDF bar-plot
        pdf_path <- file.path(outdir, paste0(contrast, "_freq_barplot.pdf"))
        grDevices::pdf(pdf_path, width = 8, height = 5)
        graphics::barplot(
          height    = freq_tbl$Frequency,
          names.arg = freq_tbl$SampleID,
          las       = 2,
          col       = "gray",
          border    = "white",
          main      = paste("Residual outliers –", contrast),
          xlab      = "SampleID",
          ylab      = "Frequency"
        )
        grDevices::dev.off()
        mdcat(paste0("- [Download bar-plot](", pdf_path, ")"))

        # Inline copy
        subchunkify({
          old_mar <- par("mar"); old_oma <- par("oma"); old_mgp <- par("mgp")
          par(mar = c(4,4,4,2)+0.1, oma = c(3,0,0,0), mgp = c(3,1,0))
          graphics::barplot(
            height    = freq_tbl$Frequency,
            names.arg = freq_tbl$SampleID,
            las       = 2,
            col       = "gray",
            border    = "white",
            main      = paste("Residual outliers –", contrast),
            xlab      = NA,
            ylab      = "Frequency"
          )
          mtext("SampleID", side = 1, line = 1, outer = TRUE)
          par(mar = old_mar, oma = old_oma, mgp = old_mgp)
        })

      } else {
        mdcat("No residual outliers found in this p-value bin.")
      }
    }
  }

  if (print_plot) {
    invisible(freq_tbl_list)
  } else {
    invisible(outlier_list)
  }
}


#' Residual QQ plots by protein or by sample
#'
#' Wrapper that either (i) delegates to `plot_qq_protein()` for per-protein,
#' p-value-binned QQ plots, or (ii) for `type = "sample"`, selects the top
#' `num_samples` by outlier frequency and, for each, draws QQ plots labeling
#' all sample outliers for selected proteins.
#'
#' @param res_list Named list from `run_tests()`; each element has `$anova`
#'   with residuals.
#' @param outdir String; output directory for PDFs (must exist).
#' @param label_map Optional list of outliers:
#'   `label_map[[contrast]][[SampleID]] = c("Assay_Panel", ...)`.
#' @param num_proteins Integer; number of proteins to show per bin/sample (default `10`).
#' @param num_samples Integer; number of samples to display (default `10`).
#' @param type `"protein"` or `"sample"` (default).
#' @param freq_table When `type = "sample"`, a data frame or list of data frames
#'   with columns `SampleID` and `Frequency` (and optionally `Contrast`).
#' @param bins Numeric vector of p-value cut points (used only for `type="protein"`).
#' @param seed Integer RNG seed.
#' @return `invisible(NULL)`; side effects are plots, PDFs, and markdown.
#' @export
plot_qq <- function(res_list,
                    outdir,
                    label_map    = NULL,
                    num_proteins = 10,
                    num_samples  = 10,
                    type         = c("protein", "sample"),
                    freq_table   = NULL,
                    bins         = c(0, 0.1, 0.25, 0.50, 0.75, 1),
                    seed         = 1) {

  type   <- match.arg(type)

  # Checks
  stopifnot(
    dir.exists(outdir),
    is.list(res_list),
    is.list(label_map) | is.null(label_map)
  )

  # Delegate to per‐protein version
  if (type == "protein") {
    plot_qq_protein(
      res_list,
      outdir,
      label_map    = label_map,
      num_proteins = num_proteins,
      bins         = bins,
      seed         = seed
    )
    return(invisible(NULL))
  }

  # type = "sample" needs freq_table
  if (is.null(freq_table))
    stop("When type = 'sample' you must supply freq_table.")

  set.seed(seed)

  # Loop over each contrast
  for (contrast_name in names(res_list)) {

    # Extract this contrast’s freq_table
    if (is.list(freq_table)) {
      ft_contrast <- freq_table[[contrast_name]]
    } else {
      ft_contrast <- subset(freq_table, Contrast == contrast_name)
    }
    if (is.null(ft_contrast) ||
        !all(c("SampleID", "Frequency") %in% names(ft_contrast))) {
      warning("Missing or malformed freq_table for contrast ", contrast_name)
      next
    }

    # Top-n samples for this contrast
    samp_order <- ft_contrast[order(-ft_contrast$Frequency,
                                    ft_contrast$SampleID), ]
    samp_ids   <- head(samp_order$SampleID, num_samples)

    # Gather ANOVA residuals
    term <- res_list[[contrast_name]]$params$column
    anova_tbl <- res_list[[contrast_name]]$anova %>%
      dplyr::filter(Effect == term)

    if (!all(c("Assay","Panel","SW_pval","resids") %in% names(anova_tbl))) {
      warning("Contrast ", contrast_name, " lacks expected columns; skipped.")
      next
    }

    mdcat(paste0("#### ", contrast_name, " {.tabset}"))

    # For each sample
    for (sid in samp_ids) {

      mdcat(paste0("##### ", sid, " {.tabset}"))

      # Proteins where this sample was flagged as outlier
      prot_panels <- if (is.list(label_map)) {
        label_map[[contrast_name]][[sid]]
      } else character(0)

      if (length(prot_panels) == 0) next

      prot_panels <- unique(prot_panels)
      if (length(prot_panels) > num_proteins)
        prot_panels <- sample(prot_panels, num_proteins)

      # Loop proteins
      for (pp in prot_panels) {

        assay <- sub("_(.*)$", "",  pp)
        panel <- sub("^[^_]+_",  "", pp)

        mdcat(paste0("###### ", assay, " (", panel, ")"))

        resid_row <- anova_tbl[
          anova_tbl$Assay == assay & anova_tbl$Panel == panel, ]
        if (nrow(resid_row) == 0) {
          mdcat("*Protein-panel not found in ANOVA table*")
          next
        }

        resid_vec <- unlist(resid_row$resids[[1]])
        if (length(resid_vec) < 3) {
          mdcat("*Not enough residuals to draw a QQ-plot*")
          next
        }

        qq_stats <- stats::qqnorm(resid_vec, plot.it = FALSE)
        plot_df  <- tibble::tibble(
          SampleID    = names(resid_vec),
          resid       = resid_vec,
          theor_quant = qq_stats$x,
          samp_quant  = qq_stats$y
        )

        # Label **all** sample outliers for this protein
        out_sids <- if (is.list(label_map)) {
          smap <- label_map[[contrast_name]]
          names(smap)[vapply(smap, function(v) pp %in% v, logical(1))]
        } else character(0)
        lbls <- intersect(out_sids, plot_df$SampleID)

        # Build QQ‐plot
        qqp <- ggplot2::ggplot(plot_df,
                               ggplot2::aes(sample = resid)) +
          ggplot2::stat_qq() +
          ggplot2::stat_qq_line() +
          ggplot2::labs(
            title = paste(
              contrast_name,
              "\nSample:",  sid,
              "\nProtein:", assay,
              "\nPanel:",   panel
            ),
            x = "Theoretical Quantiles",
            y = "Sample Quantiles"
          )

        if (length(lbls) > 0) {
          qqp <- qqp +
            ggrepel::geom_label_repel(
              data  = dplyr::filter(plot_df, SampleID %in% lbls),
              ggplot2::aes(x = theor_quant,
                           y = samp_quant,
                           label = SampleID),
              inherit.aes        = FALSE,
              color              = "black",
              fill               = "white",
              segment.color      = "black",
              treatment_1.segment.length = 0,
              label.size         = 0.25,
              box.padding        = 0.7,
              point.padding      = 0.25,
              max.overlaps       = nrow(plot_df),
              show.legend        = FALSE
            )
        }

        # Save & embed
        pdf_path <- file.path(
          outdir,
          paste0(contrast_name, "_", sid, "_", pp, "_qq_residuals.pdf")
        )
        ggplot2::ggsave(filename = pdf_path, plot = qqp)

        mdcat("- [Download PDF](", pdf_path, ")")
        subchunkify(qqp)
      }
    }
  }

  invisible(NULL)
}
