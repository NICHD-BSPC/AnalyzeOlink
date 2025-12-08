#' Append Human Protein Atlas (HPA) brain annotations to results
#'
#' @description
#' Downloads and processes HPA Normal Tissue IHC data and tissue–organ mapping,
#' computes brain-focused annotation metrics per gene, and appends these to each
#' contrast table in `res`. Returns the updated `res` alongside the raw HPA IHC
#' table and a per-gene summary used for joins.
#'
#' @details
#' Key outputs added per gene:
#' - **BrainDetected**: binary; any brain tissue Level \> 0.
#' - **BrainSpecific**: binary; top organ is Brain **and** organ-level tau ≥ 0.9.
#' - **TauBrainSpecificity / TauOrganSpecificity**: organ-level tau (0–1), with a
#'   brain-specific copy set to 0 when Brain is not the top organ.
#' - **Regional grades** *(0–3)*: Cerebellum, Caudate, CerebralCortex, Hippocampus,
#'   and **Purkinje** (cell-type in Cerebellum).
#' - **BrainPurkinjeSpecific**: Purkinje grade only when **BrainSpecific == 1**, else 0.
#'
#' Data sources (HPA TSVs) are cached under `data_dir`:
#' - `normal_ihc_data.tsv.zip`, `normal_ihc_tissues.tsv.zip`
#'
#' @param res Named list of per-contrast results; each element contains an `$anova` table
#'   with a `GeneID` column to join on.
#' @param data_dir Directory to cache HPA downloads and extracted TSVs (created if missing).
#' @return A list with elements:
#'   - `res`: updated `res` where each `$anova` gained the HPA columns.
#'   - `hpa`: raw (filtered) HPA IHC table.
#'   - `hpa_summary`: per-gene summary used for annotation joins.
#' @seealso [add_extra_hpa_data()]
#' @export
append_brain_annotations_tbl <- function(res, data_dir) {
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # Download the HPA IHC dataset and organ-tissue mapping files
  ihc_url    <- "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"
  tissue_url <- "https://www.proteinatlas.org/download/tsv/normal_ihc_tissues.tsv.zip"
  ihc_zip    <- file.path(data_dir, basename(ihc_url))
  tissue_zip <- file.path(data_dir, basename(tissue_url))

  if (!file.exists(ihc_zip))    download.file(ihc_url,    destfile = ihc_zip,    mode = "wb")
  if (!file.exists(tissue_zip)) download.file(tissue_url, destfile = tissue_zip, mode = "wb")

  ihc_files    <- utils::unzip(ihc_zip,    exdir = data_dir)
  tissue_files <- utils::unzip(tissue_zip, exdir = data_dir)

  ihc_tsv    <- file.path(data_dir, "normal_ihc_data.tsv")
  tissue_tsv <- file.path(data_dir, "normal_ihc_tissues.tsv")

  ihc    <- readr::read_tsv(ihc_tsv,    col_types = readr::cols())
  tissue <- readr::read_tsv(tissue_tsv, col_types = readr::cols())

  ihc <- dplyr::rename(ihc, GeneID = Gene, symbol = `Gene name`) %>%
    dplyr::filter(!Level %in% c("Ascending", "Descending", "N/A", "Not representative"))

  # Add Organ via tissue mapping
  hpa <- ihc %>% dplyr::left_join(tissue, by = "Tissue")

  # Map Level to integers for aggregation
  hpa_num <- hpa %>%
    dplyr::mutate(Level_num = dplyr::case_when(
      Level == "Not detected" ~ 0L,
      Level == "Low"          ~ 1L,
      Level == "Medium"       ~ 2L,
      Level == "High"         ~ 3L,
      TRUE                    ~ NA_integer_
    ))

  # Helper: max with 0 fallback (keeps 0–3 scale even if all-NA)
  safe_max <- function(x) {
    x <- x[!is.na(x)]
    if (length(x)) max(x) else 0L
  }

  # 1) Binary: BrainDetected
  hpa_1 <- hpa_num %>%
    dplyr::filter(Organ == "Brain") %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarize(BrainDetected = as.integer(any(Level_num > 0)), .groups = "drop")

  # 2) Brain specificity (organ-level tau + suppression outside brain)
  tau_thr <- 0.90

  organ_means <- hpa_num %>%
    dplyr::group_by(GeneID, Organ) %>%
    dplyr::summarize(org_mean = mean(Level_num, na.rm = TRUE), .groups = "drop")

  tau_tbl <- organ_means %>%
    dplyr::group_by(GeneID) %>%
    dplyr::arrange(dplyr::desc(org_mean), .by_group = TRUE) %>%
    dplyr::summarize(
      top_organ = dplyr::first(Organ),
      max_mean  = dplyr::first(org_mean),
      k         = dplyr::n(),
      tau = dplyr::if_else(max_mean > 0 & k > 1,
                           sum(1 - (org_mean / max_mean)) / (k - 1),
                           0),
      .groups = "drop"
    )

  brain_prop_tbl <- hpa_num %>%
    dplyr::filter(Organ == "Brain") %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarize(
      brain_prop_high = mean(Level_num >= 2L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(brain_prop_high = ifelse(is.nan(brain_prop_high), 0, brain_prop_high))

  extra_tbl <- hpa_num %>%
    dplyr::filter(Organ != "Brain") %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarize(
      extra_max      = max(Level_num, na.rm = TRUE),
      extra_hi_count = sum(Level_num >= 2L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      extra_max      = ifelse(is.infinite(extra_max), 0, extra_max),
      extra_hi_count = ifelse(is.na(extra_hi_count), 0L, extra_hi_count)
    )

  hpa_2 <- tau_tbl %>%
    dplyr::left_join(brain_prop_tbl, by = "GeneID") %>%
    dplyr::left_join(extra_tbl,       by = "GeneID") %>%
    dplyr::mutate(
      brain_prop_high = dplyr::coalesce(brain_prop_high, 0),
      extra_max       = dplyr::coalesce(extra_max, 0),
      extra_hi_count  = dplyr::coalesce(extra_hi_count, 0L),

      BrainSpecific   = as.integer(top_organ == "Brain" & tau >= tau_thr),

      TauOrganSpecificity = as.numeric(tau),

      TauBrainSpecificity = dplyr::if_else(top_organ == "Brain",
                                           TauOrganSpecificity,
                                           0.0),

      TopOrgan = top_organ
    ) %>%
    dplyr::mutate(
      TauOrganSpecificity = round(TauOrganSpecificity, 3),
      TauBrainSpecificity = round(TauBrainSpecificity, 3)
    ) %>%
    dplyr::select(GeneID, BrainSpecific, TauBrainSpecificity, TauOrganSpecificity, TopOrgan)

  # 3) Regional (0–3) brain tissue grades + Purkinje
  brain_levels <- hpa_num %>%
    dplyr::filter(Organ == "Brain") %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarize(
      Cerebellum     = safe_max(Level_num[Tissue == "Cerebellum"]),
      Caudate        = safe_max(Level_num[Tissue == "Caudate"]),
      CerebralCortex = safe_max(Level_num[Tissue == "Cerebral cortex"]),
      Hippocampus    = safe_max(Level_num[Tissue == "Hippocampus"]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(-GeneID, as.integer))

  purkinje_levels <- hpa_num %>%
    dplyr::filter(Organ == "Brain", Tissue == "Cerebellum") %>%
    dplyr::mutate(
      .has_celltype = "Cell type" %in% names(.),
      .is_purkinje  = ifelse(.has_celltype, grepl("Purkinje", .[["Cell type"]], ignore.case = TRUE), FALSE)
    ) %>%
    dplyr::filter(.is_purkinje) %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarize(Purkinje = safe_max(Level_num), .groups = "drop") %>%
    dplyr::mutate(Purkinje = as.integer(Purkinje))

  hpa_region_grades <- brain_levels %>%
    dplyr::left_join(purkinje_levels, by = "GeneID") %>%
    dplyr::mutate(Purkinje = dplyr::coalesce(Purkinje, 0L))

  # 4) Combine and coerce types
  hpa_summary <- hpa_1 %>%
    dplyr::left_join(hpa_2,             by = "GeneID") %>%
    dplyr::left_join(hpa_region_grades, by = "GeneID") %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("Purkinje","Cerebellum","Caudate","CerebralCortex","Hippocampus")),
        ~ as.integer(dplyr::coalesce(., 0L))
      ),
      BrainDetected = as.integer(dplyr::coalesce(BrainDetected, 0L)),
      BrainSpecific = as.integer(dplyr::coalesce(BrainSpecific, 0L))
    ) %>%
    dplyr::mutate(
      BrainPurkinjeSpecific = as.integer(dplyr::if_else(BrainSpecific == 1L & Purkinje > 0L, Purkinje, 0L))
    ) %>% add_extra_hpa_data(data_dir)

  # 5) Append to each contrast's ANOVA
  for (nm in names(res)) {
    annotate_cols_int <- c("BrainDetected","BrainSpecific",
                           "Purkinje","Cerebellum","Caudate",
                           "CerebralCortex","Hippocampus",
                           "BrainPurkinjeSpecific")

    res[[nm]]$anova <- res[[nm]]$anova %>%
      dplyr::left_join(hpa_summary, by = "GeneID") %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(annotate_cols_int),
          ~ as.integer(dplyr::coalesce(., 0L))
        ),
        TauBrainSpecificity = dplyr::coalesce(TauBrainSpecificity, 0),
        TauOrganSpecificity = dplyr::coalesce(TauOrganSpecificity, 0)
      )
  }

  res <- sort_all_results(res)
  return(list(res = res, hpa_data = hpa, hpa_summary = hpa_summary))
}


#' Add extra HPA-derived annotations (subcellular location, blood PEA DE, interactions)
#'
#' @description
#' Augments a per-gene HPA summary with:
#' 1) **SubCellLocation** (collapsed main/additional/extracellular locations),
#' 2) **Top10BloodDisease** (Healthy-control blood PEA DE: top 10 diseases by adjusted p-value),
#' 3) **Interactions** (consensus protein–protein partners mapped to gene symbols).
#'
#' @param hpa_summary Tibble/data frame with a `GeneID` column (Ensembl).
#' @param data_dir Directory used to cache HPA downloads (created if missing).
#' @return The input `hpa_summary` with new columns appended (where available).
#' @note This function downloads public TSVs from the Human Protein Atlas.
#' @seealso [append_brain_annotations_tbl()]
#' @export
add_extra_hpa_data <- function(hpa_summary, data_dir) {
  stopifnot(is.data.frame(hpa_summary), "GeneID" %in% names(hpa_summary))
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # helpers
  .safe_download <- function(url, dest) {
    if (!file.exists(dest)) utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    dest
  }
  .unzip_pick_tsv <- function(zip_path, pattern_hint = NULL) {
    extracted <- utils::unzip(zip_path, exdir = data_dir)
    tsvs <- extracted[grepl("\\.tsv$", extracted, ignore.case = TRUE)]
    if (!is.null(pattern_hint)) {
      hit <- tsvs[grepl(pattern_hint, basename(tsvs), ignore.case = TRUE)]
      if (length(hit)) return(hit[[1]])
    }
    if (length(tsvs)) return(tsvs[[1]])
    stop("No .tsv found in: ", zip_path)
  }
  .collapse_unique <- function(x, in_sep = ";", out_sep = ", ") {
    if (length(x) == 0 || all(is.na(x))) return(NA_character_)
    toks <- unlist(strsplit(paste(stats::na.omit(x), collapse = in_sep), in_sep, fixed = TRUE))
    toks <- trimws(toks); toks <- toks[nzchar(toks)]; toks <- unique(toks)
    if (!length(toks)) return(NA_character_)
    paste(toks, collapse = out_sep)
  }

  # 1) Subcellular location
  sub_url <- "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip"
  sub_zip <- file.path(data_dir, basename(sub_url))
  .safe_download(sub_url, sub_zip)
  sub_tsv <- .unzip_pick_tsv(sub_zip, pattern_hint = "subcellular")

  sub_raw <- readr::read_tsv(sub_tsv, col_types = readr::cols())
  req_sub <- c("Gene", "Gene name", "Main location", "Additional location", "Extracellular location")
  miss_sub <- setdiff(req_sub, names(sub_raw))
  if (length(miss_sub)) stop("Subcellular TSV missing: ", paste(miss_sub, collapse = ", "))

  sub_map <- sub_raw %>%
    dplyr::transmute(GeneID = .data$Gene, Symbol = .data$`Gene name`)

  sub_loc <- sub_raw %>%
    dplyr::rename(GeneID = .data$Gene) %>%
    dplyr::mutate(
      `Main location`          = dplyr::na_if(`Main location`, ""),
      `Additional location`    = dplyr::na_if(`Additional location`, ""),
      `Extracellular location` = dplyr::na_if(`Extracellular location`, "")
    ) %>%
    tidyr::unite("SubCellLocation",
                 c("Main location","Additional location","Extracellular location"),
                 sep = ";", remove = TRUE, na.rm = TRUE) %>%
    dplyr::group_by(.data$GeneID) %>%
    dplyr::summarize(
      SubCellLocation = .collapse_unique(.data$SubCellLocation, in_sep = ";", out_sep = ", "),
      .groups = "drop"
    )

  # 2) Blood PEA DE (Healthy subset; top 10 per gene)
  pea_url <- "https://www.proteinatlas.org/download/tsv/blood_pea_disease_de.tsv.zip"
  pea_zip <- file.path(data_dir, basename(pea_url))
  .safe_download(pea_url, pea_zip)
  pea_tsv <- .unzip_pick_tsv(pea_zip, pattern_hint = "blood")

  pea_raw <- readr::read_tsv(pea_tsv, col_types = readr::cols())
  req_pea <- c("Gene", "ENSG ID", "Disease", "Class", "Control", "p-value adjusted")
  miss_pea <- setdiff(req_pea, names(pea_raw))
  if (length(miss_pea)) stop("Blood PEA TSV missing: ", paste(miss_pea, collapse = ", "))

  pea_map <- pea_raw %>%
    dplyr::transmute(GeneID = .data$`ENSG ID`, Symbol = .data$Gene)

  pea_top10 <- pea_raw %>%
    dplyr::filter(.data$Control == "Healthy") %>%
    dplyr::rename(Assay = .data$Gene,
                  GeneID = .data$`ENSG ID`,
                  padj   = .data$`p-value adjusted`) %>%
    dplyr::mutate(padj = suppressWarnings(as.numeric(.data$padj))) %>%
    dplyr::group_by(.data$GeneID) %>%
    dplyr::arrange(.data$padj, .by_group = TRUE) %>%
    dplyr::distinct(.data$Disease, .data$Class, .keep_all = TRUE) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::summarize(
      Top10BloodDisease = paste0(.data$Disease, " (", .data$Class, ")", collapse = ", "),
      .groups = "drop"
    )

  # 3) Interaction consensus → gene-symbol partners
  ensg_to_assay <- dplyr::bind_rows(
      pea_map %>% dplyr::mutate(.priority = 1L),
      sub_map %>% dplyr::mutate(.priority = 2L)
    ) %>%
    dplyr::filter(!is.na(.data$GeneID), nzchar(.data$GeneID)) %>%
    dplyr::arrange(.data$GeneID, .data$.priority) %>%
    dplyr::group_by(.data$GeneID) %>%
    dplyr::summarize(Assay = dplyr::first(.data$Symbol[nzchar(.data$Symbol)]), .groups = "drop")

  int_url <- "https://www.proteinatlas.org/download/tsv/interaction_consensus.tsv.zip"
  int_zip <- file.path(data_dir, basename(int_url))
  .safe_download(int_url, int_zip)
  int_tsv <- .unzip_pick_tsv(int_zip, pattern_hint = "interaction")

  int_raw <- readr::read_tsv(int_tsv, col_types = readr::cols())
  req_int <- c("ensembl_gene_id_1", "ensembl_gene_id_2")
  miss_int <- setdiff(req_int, names(int_raw))
  if (length(miss_int)) stop("Interaction consensus TSV missing: ", paste(miss_int, collapse = ", "))

  pairs <- int_raw %>%
    dplyr::select(a = .data$ensembl_gene_id_1, b = .data$ensembl_gene_id_2) %>%
    dplyr::filter(!is.na(.data$a), !is.na(.data$b), nzchar(.data$a), nzchar(.data$b))

  int_sym <- dplyr::bind_rows(
      pairs %>% dplyr::rename(GeneID = .data$a, partner = .data$b),
      pairs %>% dplyr::rename(GeneID = .data$b, partner = .data$a)
    ) %>%
    dplyr::left_join(ensg_to_assay, by = c("partner" = "GeneID")) %>%
    dplyr::mutate(partner_assay = dplyr::coalesce(.data$Assay, .data$partner)) %>%
    dplyr::group_by(.data$GeneID) %>%
    dplyr::summarize(
      Interactions = .collapse_unique(.data$partner_assay, in_sep = ";", out_sep = ", "),
      .groups = "drop"
    )

  # Join extras to hpa_summary
  out <- hpa_summary %>%
    dplyr::left_join(sub_loc,   by = "GeneID") %>%
    dplyr::left_join(pea_top10, by = "GeneID") %>%
    dplyr::left_join(int_sym,   by = "GeneID")

  tibble::as_tibble(out)
}


#' Export raw HPA IHC table to an Excel workbook
#'
#' @description
#' Writes the (renamed) HPA IHC table to one or more Excel sheets, splitting
#' into chunks to avoid Excel row limits.
#'
#' @param hpa HPA IHC data frame (as returned from `append_brain_annotations_tbl()`).
#' @param outdir Output directory for the Excel file (created if missing).
#' @return Invisibly, the path to the written workbook.
#' @export
export_hpa_ihc_data <- function(hpa, outdir) {
  # ensure output directory exists
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # rename columns
  hpa_renamed <- hpa %>%
    dplyr::rename(
      IHCTissue = `IHC tissue name`,
      CellType  = `Cell type`
    )

  max_rows   <- 250000L
  total_rows <- nrow(hpa_renamed)
  n_sheets   <- ceiling(total_rows / max_rows)

  wb <- openxlsx::createWorkbook()
  for (i in seq_len(n_sheets)) {
    start_row <- (i - 1) * max_rows + 1L
    end_row   <- min(i * max_rows, total_rows)
    sheet_name <- if (n_sheets == 1) "HPA_IHC" else paste0("HPA_IHC_Part", i)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x     = hpa_renamed[start_row:end_row, , drop = FALSE],
      withFilter = TRUE
    )
  }

  out_file <- file.path(outdir, "hpa_tissue_ihc_data.xlsx")
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)

  message(sprintf("HPA IHC data exported to %s", out_file))
  invisible(out_file)
}


#' Format annotated results tables for export
#'
#' @description
#' Produces a cleaned, presentation-ready table from an analysis results table,
#' with optional handling for therapeutic-reversal (two-contrast) tables when a
#' `reversal` column and `estimate_*` columns are present.
#'
#' @param df Data frame of results (single-contrast or therapeutic-reversal).
#' @return A formatted data frame with selected/rounded columns and annotations
#'   (if present).
#' @note This function currently accepts a single argument `df`. Ensure callers
#'   pass only one argument (some code paths may still call with an extra `hpa`).
#' @export
format_annotations_results <- function(df) {
  # Drop columns that are entirely NA
  df <- df %>% dplyr::select(tidyselect::where(~ !all(is.na(.x))))

  # Detect therapeutic-reversal (joined two-contrast) tables
  is_tr_table <- (
    "reversal" %in% names(df) &&
    any(grepl("^estimate_", names(df)))
  )

  if (is_tr_table) {
    # Therapeutic reversal formatting
    est_cols   <- grep("^estimate_",       names(df), value = TRUE)
    sig_cols   <- grep("^sig_",            names(df), value = TRUE)
    adj_p_cols <- grep("^Adjusted_pval_",  names(df), value = TRUE)
    p_cols     <- grep("^p\\.value_",      names(df), value = TRUE)

    df %>%
      dplyr::select(
        -dplyr::matches("(padj)$"),
        -dplyr::any_of(c(
          "resids","AIC","BIC","AICc","symbol","AD_pval","LT_pval",
          "y_plot","shape","Assay_Panel","Threshold","status"
        ))
      ) %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::any_of(c(est_cols, "pes","t.ratio","SW_pval","FT_pval",
            "TauBrainSpecificity", "TauOrganSpecificity")),
          ~ round(.x, 3)
        ),
        dplyr::across(
          dplyr::any_of(c("p.value","Adjusted_pval", adj_p_cols, p_cols)),
          ~ formatC(as.numeric(.x), format = "e", digits = 2)
        )
      ) %>%
      dplyr::select(
        dplyr::any_of(c(
          "Assay","GeneID","OlinkID","UniProt","Panel",
          est_cols, "rank_abs", p_cols, adj_p_cols, sig_cols, "reversal_sig_both",
          "significant_in_both","reversal",
          "BrainDetected","BrainSpecific","TauBrainSpecificity",
          "TopOrgan","TauOrganSpecificity",
          "BrainPurkinjeSpecific","Cerebellum","Caudate","CerebralCortex","Hippocampus",
          "Purkinje","SubCellLocation","Top10BloodDisease","Interactions"
        )),
        dplyr::everything()
      ) %>%
      dplyr::arrange(dplyr::desc(reversal_sig_both), dplyr::desc(rank_abs))

  } else {
    # Standard single-contrast tables
    df %>%
      dplyr::select(
        -dplyr::matches("(padj)$"),
        -dplyr::any_of(c("resids","AIC","BIC","AICc","symbol","AD_pval","LT_pval"))
      ) %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::any_of(c("p.value","Adjusted_pval")),
          ~ formatC(as.numeric(.x), format = "e", digits = 2)
        ),
        dplyr::across(
          dplyr::any_of(c("pes","log2FC","estimate","t.ratio","SW_pval","FT_pval",
            "TauBrainSpecificity", "TauOrganSpecificity")),
          ~ round(.x, 3)
        )
      ) %>%
      dplyr::select(
        dplyr::any_of(c(
          "Assay","GeneID","OlinkID","UniProt","Panel","Effect","p.value",
          "Adjusted_pval","pes","log2FC","estimate","Threshold",
          "formula","df","MSE","F","t.ratio","SW_pval","FT_pval",
          "BrainDetected","BrainSpecific","TauBrainSpecificity",
          "TopOrgan","TauOrganSpecificity",
          "BrainPurkinjeSpecific","Cerebellum","Caudate","CerebralCortex","Hippocampus",
          "Purkinje","SubCellLocation","Top10BloodDisease","Interactions"
        )),
        dplyr::everything()
      )
  }
}


#' Write per-contrast annotation results to Excel workbooks
#'
#' @description
#' For each contrast in `res`, formats each results table via
#' `format_annotations_results()` and writes them to an Excel workbook with one
#' worksheet per analysis/test type.
#'
#' @param res Named list of contrast results; each element has nested data frames
#'   for test types (excluding `"emm_int"` and `"params"`).
#' @param hpa (Unused) Reserved for future use; currently ignored.
#' @param config Configuration list with `output_paths$annotations$results`.
#' @return A named character vector of workbook paths (one per contrast).
#' @export
export_annotations <- function(res, hpa, config) {
  outdir  <- config$output_paths$annotations$results

  workbook_paths <- character()

  for (contrast in names(res)) {
    wb   <- openxlsx::createWorkbook()
    keep <- FALSE

    for (test_type in base::setdiff(names(res[[contrast]]), c("emm_int", "params"))) {
      df <- res[[contrast]][[test_type]]
      if (is.null(df) || nrow(df) == 0) next

      sheet_name <- get_test_name(test_type)
      df_fmt     <- format_annotations_results(df)  # NOTE: format_* currently accepts only `df`

      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = df_fmt, withFilter = TRUE)
      keep <- TRUE
    }

    if (keep) {
      file_name <- paste0(contrast, "_annotations_results.xlsx")
      path      <- file.path(outdir, file_name)
      openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
      workbook_paths[contrast] <- path
    }
  }

  return(workbook_paths)
}


#' Parallel export of annotation results (one workbook per contrast)
#'
#' @description
#' Parallelized variant of [export_annotations()] using `parallel::mclapply`.
#' Each contrast is written to its own workbook if it contains at least one
#' non-empty results table.
#'
#' @param res Named list of contrast results.
#' @param outdir Output directory for the workbooks.
#' @return A named character vector of workbook paths (one per contrast),
#'   omitting contrasts with no written sheets.
#' @export
parallel_export_annotations <- function(res, outdir) {
  contrasts <- names(res)
  # Run one contrast per core
  paths_list <- parallel::mclapply(contrasts, function(contrast) {
    wb   <- openxlsx::createWorkbook()
    keep <- FALSE

    # Loop over test types (skip emm_int and params)
    for (test_type in base::setdiff(names(res[[contrast]]), c("emm_int", "params"))) {
      df <- res[[contrast]][[test_type]]
      if (is.null(df) || nrow(df) == 0) next

      sheet_name <- get_test_name(test_type)
      df_fmt     <- format_annotations_results(df)

      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = df_fmt, withFilter = TRUE)
      keep <- TRUE
    }

    # Save workbook if any sheets were written
    if (keep) {
      file_name <- paste0(contrast, "_annotations_results.xlsx")
      path      <- file.path(outdir, file_name)
      openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
      return(path)
    } else {
      return(NULL)
    }
  }, mc.cores = length(contrasts))

  # Flatten and name the result vector
  paths <- unlist(paths_list, use.names = FALSE)
  names(paths) <- contrasts[!vapply(paths_list, is.null, logical(1))]

  return(paths)
}


#' Build an annotated volcano plot (PDF + interactive Plotly)
#'
#' @description
#' Creates a volcano plot for a single contrast + test type, colored by an
#' annotation column (binary 0/1 or graded 0–3). Writes a PDF and returns the
#' ggplot object, a Plotly widget, and output paths.
#'
#' @param res Results data frame containing at least `Adjusted_pval` and either
#'   `estimate` (preferred) or `log2FC`.
#' @param se_list Named list of `SummarizedExperiment` objects (unused here, reserved).
#' @param contrast Contrast name (string) used in titles/filenames.
#' @param test_type Test name key (used for sheet titles and filenames).
#' @param outdir Output directory (created if needed).
#' @param parameters List with `treatment`, `control`, and `alpha`.
#' @param color_by Annotation column to map to color (supports 0/1 or 0–3 integers).
#' @param export_2nd Logical; also write a zoomed PDF limited by `export_2nd_y_lim`.
#' @param export_2nd_y_lim Numeric; y-axis (–log10 p) limit for the zoomed plot.
#' @param vert_line Logical; draw a vertical line at x = 0.
#' @return Invisibly, a list with `ggplot`, `p_plotly`, `pdf_primary`, `pdf_secondary`.
#' @export
annotations_make_volcano <- function(res,
                                    se_list,
                                    contrast,
                                    test_type,
                                    outdir,
                                    parameters,
                                    color_by,
                                    export_2nd       = FALSE,
                                    export_2nd_y_lim = 5,
                                    vert_line        = TRUE) {
  # Prepare x‐ and y‐axes
  if ("estimate" %in% names(res)) {
    x_var <- "estimate"; x_lab <- "estimate (adj. log2FC)"
  } else if ("log2FC" %in% names(res)) {
    x_var <- "log2FC";    x_lab <- "log2FC"
  } else {
    stop("Need either ‘estimate’ or ‘log2FC’ in results")
  }
  # Check color_by exists
  if (!color_by %in% names(res)) {
    stop("Column ‘", color_by, "’ not found in results")
  }

  # Determine whether the color_by column is binary (0/1) or graded (0–3)
  vals <- sort(unique(res[[color_by]]))
  is_binary <- all(vals %in% c(0L, 1L))
  color_levels <- if (is_binary) c(0L, 1L) else c(0L, 1L, 2L, 3L)

  # Compute –log10 p‐value and create grouping factor (supports 0/1 and 0–3)
  res <- dplyr::mutate(
    res,
    y           = -log10(Adjusted_pval),
    color_group = factor(as.integer(.data[[color_by]]), levels = color_levels)
  )
  present_levels <- levels(res$color_group)

  # Palette: binary uses two colors; graded uses four
  pal_all <- if (is_binary) {
    c("0" = "steelblue", "1" = "indianred")
  } else {
    c("0" = "grey70", "1" = "steelblue", "2" = "goldenrod", "3" = "indianred")
  }
  pal_use <- pal_all[present_levels]

  # Build ggplot
  p <- ggplot2::ggplot(
    res,
    ggplot2::aes(x = .data[[x_var]], y = y, color = color_group)
  )

  # Draw points in explicit z-order as: grey -> blue -> yellow -> red on top
  # present_levels already defined above; e.g., c("0","1","2","3") or c("0","1")
  layer_order <- intersect(c("0","1","2","3"), present_levels)
  for (lvl in layer_order) {
    p <- p + ggplot2::geom_point(
      data  = res[res$color_group == lvl, , drop = FALSE],
      na.rm = TRUE
    )
  }

  # Manual colors, legend order, titles, etc.
  p <- p +
    ggplot2::scale_color_manual(
      values = pal_use,
      breaks = present_levels,
      name   = color_by
    ) +
    ggplot2::labs(
      x        = x_lab,
      y        = "-log10(adj. p-value)",
      title    = paste("Volcano:", test_type, contrast),
      subtitle = sprintf(
        "%s vs %s (alpha = %g)",
        parameters$treatment, parameters$control, parameters$alpha
      )
    )
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    if (vert_line) {
      p <- p + ggplot2::geom_vline(
        xintercept = 0,
        linetype   = "dashed"
      )
    }
  p <- p +
    ggplot2::geom_hline(
      yintercept = -log10(parameters$alpha),
      linetype   = "dashed"
    ) +
    ggplot2::annotate(
      "text",
      x     = Inf,
      y     = -log10(parameters$alpha),
      label = paste0("-log10(", parameters$alpha, ")"),
      hjust = 1.1, vjust = -0.5, size = 3
    )

  # Interactive Plotly with tooltip and link-out to HPA
  p_plotly <- plotly::ggplotly(
    p +
      ggplot2::aes(
        text = paste0(
          "Protein: ", Assay,
          "<br>GeneID: ", GeneID,
          "<br>Panel: ", Panel,
          "<br>", x_lab, ": ", round(.data[[x_var]], 3),
          "<br>Adjusted p-value: ", sprintf("%.2e", Adjusted_pval),
          "<br>Sub Cellular Location: ", SubCellLocation,
          "<br>Top Diseases Causing Increase in Blood: ",
            stringr::str_trunc(Top10BloodDisease, width = 50, side = "right", ellipsis = "…"),
          "<br>Protein-Protein Interactions With: ",
            stringr::str_trunc(Interactions, width = 50, side = "right", ellipsis = "…")
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
          y0    = -log10(parameters$alpha),
          y1    = -log10(parameters$alpha),
          line  = list(dash = "dash", color = "black")
        )
      ),
      annotations = list(
        list(
          x         = 1,
          y         = -log10(parameters$alpha),
          text      = paste0("-log10(", parameters$alpha, ")"),
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
          var url   = 'https://www.proteinatlas.org/' + gene + '-' + assay;
          window.open(url, '_blank');
        });
      }"
    )

  # save PDF(s)
  pdf_primary <- file.path(outdir, paste0(test_type, "_", contrast, "_", color_by, "_annotations_volcano.pdf"))
  ggplot2::ggsave(pdf_primary, p, width = 7, height = 5)

  pdf_secondary <- NULL
  if (export_2nd) {
    p_zoom <- p + ggplot2::coord_cartesian(ylim = c(0, export_2nd_y_lim + 0.1))
    pdf_secondary <- file.path(
      outdir,
      paste0(test_type, "_", contrast, "_", color_by, "_annotations_volcano_zoom.pdf")
    )
    ggplot2::ggsave(pdf_secondary, p_zoom, width = 7, height = 5)
  }

  invisible(list(
    ggplot        = p,
    p_plotly      = p_plotly,
    pdf_primary   = pdf_primary,
    pdf_secondary = pdf_secondary
  ))
}


#' Knit volcano plots for all contrasts/test types
#'
#' @description
#' For each contrast and test type in `res_list`, renders a Plotly volcano plot
#' colored by `color_by` within the R Markdown report and links to the exported
#' Excel workbook and PDFs.
#'
#' @param res_list Named list of contrast results; each element includes results
#'   tables and a `params` list with `name`, `se_object_name`, `column`,
#'   `treatment`, `control`, and `alpha`.
#' @param se_list Named list of `SummarizedExperiment` objects (passed through).
#' @param outdir Directory for PDF outputs.
#' @param excel_paths Named vector of workbook paths per contrast.
#' @param color_by Column name used for color grouping in volcano plots.
#' @return `invisible(NULL)`; called for its side effects (markdown + files).
#' @export
annotations_volcano <- function(res_list, se_list, outdir, excel_paths, color_by) {
  for (contrast in names(res_list)) {
    mdcat(paste("###", contrast, "{.tabset}"))

    for (test_type in base::setdiff(names(res_list[[contrast]]), c("emm_int", "params"))) {
      test_header <- get_test_name(test_type)
      mdcat(paste("####", test_header, "Volcano Plot"))
      res <- res_list[[contrast]][[test_type]]
      parameters <- res_list[[contrast]]$params

      stopifnot(all(c("Assay", "estimate", "Adjusted_pval") %in% names(res)))
      required_params <- c("name", "se_object_name", "column", "treatment", "control", "alpha")
      if (!all(required_params %in% names(parameters))) {
        stop("contrast_list[['", contrast, "']] must contain: ", paste(required_params, collapse = ", "))
      }

      plot_n_path <- annotations_make_volcano(res, se_list, contrast, test_type, outdir, parameters, color_by)
      # Excel
      mdcat(paste0("- [", test_header, " ", contrast, " Excel data](", excel_paths[[contrast]], ")"))
      # PDF volcano
      mdcat(paste0("- [", test_header, " ", contrast, " Volcano Plot PDF](", plot_n_path$pdf_primary, ")"))
      if (contrast == "disease_vs_control") {
        mdcat(paste0("- [", test_header, " ", contrast, " Zoomed-in Volcano Plot PDF](", plot_n_path$pdf_secondary, ")"))
      }
      subchunkify(plot_n_path$p_plotly)
    }
  }
}


#' Two-contrast estimate scatter with brain-annotation coloring
#'
#' @description
#' Compares `estimate` across two selected contrasts for proteins significant in
#' both, colored by a brain-annotation column (supports 0/1 or 0–3). Writes a
#' PDF and optionally renders an interactive version.
#'
#' @param res Named list of contrast results (at least two).
#' @param outdir Output directory for the PDF.
#' @param contrasts Character length-2: the two contrast names to compare.
#' @param label_concord Logical; label top proteins by |x| + |y|.
#' @param label_sig_both Logical; (mutually exclusive with `label_concord`).
#' @param plotly Logical; also render an interactive Plotly widget.
#' @param color_by Annotation column used for color grouping.
#' @param y.lim Numeric length-2; y-axis limits (clipping is indicated).
#' @param label_n Integer; number of labels when `label_concord = TRUE`.
#' @return `invisible(NULL)`; side effects are a saved PDF and optional widget.
#' @export
annotations_estimate_plot <- function(res,
                                      outdir,
                                      contrasts,
                                      label_concord  = TRUE,
                                      label_sig_both = FALSE,
                                      plotly         = FALSE,
                                      color_by,
                                      y.lim = c(-2, 2),
                                      label_n = 10
) {

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

  pull_estimate <- function(ctr) {
    an_tbl      <- res[[ctr]]$anova
    effect_name <- res[[ctr]]$params$column
    an_sub      <- dplyr::filter(an_tbl, Effect == effect_name)

    if (!"estimate" %in% names(an_sub)) {
      stop("Contrast '", ctr, "' ANOVA does not have an 'estimate' column.")
    }
    if (!"Threshold" %in% names(an_sub)) {
      stop("Contrast '", ctr, "' ANOVA does not have a 'Threshold' column.")
    }
    if (!color_by %in% names(an_sub)) {
      stop("Contrast '", ctr, "' ANOVA does not have the required annotation column: ", color_by)
    }

    an_sub %>%
      dplyr::select(Assay, GeneID, Panel, estimate, Threshold,
                    SubCellLocation, Top10BloodDisease, Interactions,
                    !!rlang::sym(color_by)) %>%
      dplyr::rename(
        value = estimate,
        thr   = Threshold,
        ann   = !!rlang::sym(color_by)
      ) %>%
      dplyr::mutate(
        sig      = (thr == "Significant"),
        Contrast = ctr
      )
  }

  df1 <- pull_estimate(contrasts[1])
  df2 <- pull_estimate(contrasts[2])

  joined <- dplyr::full_join(
    df1, df2,
    by     = c("Assay", "Panel", "GeneID"),
    suffix = c(".A", ".B")
  ) %>%
    dplyr::filter(sig.A & sig.B, !is.na(ann.A), !is.na(ann.B)) %>%
    dplyr::mutate(
      color_group = {
        vals <- sort(unique(as.integer(ann.A)))
        is_binary <- all(vals %in% c(0L, 1L))
        lvls <- if (is_binary) c(0L, 1L) else c(0L, 1L, 2L, 3L)
        factor(as.integer(ann.A), levels = lvls)
      },
      y_plot = dplyr::case_when(
        value.B < y.lim[1] ~ y.lim[1],
        value.B > y.lim[2] ~ y.lim[2],
        TRUE               ~ value.B
      )
    )

  # palettes for binary vs graded columns
  present_levels <- levels(joined$color_group)
  is_binary_plot <- identical(present_levels, c("0", "1"))
  pal_all <- if (is_binary_plot) {
    c("0" = "steelblue", "1" = "indianred")
  } else {
    c("0" = "grey70", "1" = "steelblue", "2" = "goldenrod", "3" = "indianred")
  }
  pal_use <- pal_all[present_levels]

  gg <- ggplot2::ggplot(
    joined,
    ggplot2::aes(x = value.A, y = y_plot, color = color_group)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  # Draw points in explicit z-order
  # present_levels are defined above as levels(joined$color_group)
  layer_order <- intersect(c("0","1","2","3"), present_levels)
  for (lvl in layer_order) {
    gg <- gg + ggplot2::geom_point(
      data  = joined[joined$color_group == lvl, , drop = FALSE],
      size  = 2,
      alpha = 0.7,
      na.rm = TRUE
    )
  }

  gg <- gg +
    ggplot2::scale_color_manual(values = pal_use, breaks = present_levels, name = color_by) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::scale_y_continuous(limits = y.lim) +
    ggplot2::labs(
      x     = paste0("estimate (", contrasts[1], ")"),
      y     = paste0("estimate (", contrasts[2], ")"),
      title = paste0("estimate: ", contrasts[1], " vs ", contrasts[2])
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (isTRUE(label_concord)) {
    to_label <- joined %>%
      dplyr::mutate(rank_abs = abs(value.A) + abs(value.B)) %>%
      dplyr::slice_max(rank_abs, n = label_n, with_ties = FALSE)

    gg <- gg +
      ggrepel::geom_label_repel(
        data          = to_label,
        mapping       = ggplot2::aes(label = Assay),
        inherit.aes   = TRUE,
        color         = "black",
        fill          = "white",
        segment.color = "black",
        label.size    = 0.25,
        box.padding   = 0.3,
        point.padding = 0.25,
        max.overlaps  = nrow(joined),
        show.legend   = FALSE
      )
  }

  out_file <- file.path(
    outdir,
    paste0(paste(contrasts, collapse = "_"), "_", color_by, "_annotations_estimate_plot.pdf")
  )
  ggplot2::ggsave(filename = out_file, plot = gg)

  mdcat(paste0(
    "- [estimate Plot: ", contrasts[1], " vs ", contrasts[2], "](", out_file, ")"
  ))

  if (isTRUE(plotly)) {
    ply <- plotly::ggplotly(
      gg +
        ggplot2::aes(
          text = paste0(
            "Protein: ", Assay,
            "<br>GeneID: ", GeneID,
            "<br>Panel: ", Panel,
            "<br>", color_by, ": ", color_group,
            "<br>estimate(", contrasts[1], "): ", round(value.A, 3),
            "<br>estimate(", contrasts[2], "): ", round(value.B, 3),
            "<br>Sub Cellular Location: ", SubCellLocation.A,
            "<br>Top Diseases Causing Increase in Blood: ", Top10BloodDisease.A,
            "<br>Protein-Protein Interactions With: ", Interactions.A
          ),
          customdata = paste(GeneID, Assay, sep = "|")
        ),
      tooltip = "text"
    ) %>%
      htmlwidgets::onRender(
        "function(el, x) {
           el.on('plotly_click', function(data) {
             var pt    = data.points[0];
             var parts = pt.customdata.split('|');
             var gene  = parts[0];
             var assay = parts[1];
             var url   = 'https://www.proteinatlas.org/' + gene + '-' + assay;
             window.open(url, '_blank');
           });
         }"
      )
    subchunkify(ply, fig.height = 10, fig.width = 10)
  } else {
    print(gg)
  }

  invisible(NULL)
}


#' Internal: reduce duplicated assay flags to any(==1)
#' @keywords internal
#' @noRd
.any01 <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(any(x == 1, na.rm = TRUE))
}

#' Internal: direction call from padj/estimate
#' @param padj Numeric or factor coercible to numeric adjusted p-value.
#' @param est Numeric or factor coercible to numeric effect estimate.
#' @param thr Numeric FDR cutoff (default 0.10).
#' @return Factor with levels `c("up","down","no_change")`.
#' @keywords internal
#' @noRd
.dir_from <- function(padj, est, thr = 0.10) {
  # Coerce factors safely
  if (is.factor(padj)) padj <- as.numeric(as.character(padj))
  if (is.factor(est))  est  <- as.numeric(as.character(est))

  out <- ifelse(!is.na(padj) & padj <= thr & !is.na(est) & est > 0, "up",
         ifelse(!is.na(padj) & padj <= thr & !is.na(est) & est < 0, "down",
                "no_change"))
  base::factor(out, levels = c("up","down","no_change"))
}

#' Brain-annotation heatmap (regions × assays)
#'
#' @description
#' Builds an interactive (HTML) and static (PDF) heatmap summarizing brain-region
#' presence flags per assay, with row annotations for panel, disease direction,
#' optional treatment_1 direction, and brain detection status.
#'
#' @param res Results data frame for a single contrast containing region columns
#'   (e.g., `HigherInCaudate`, `HigherInHippocampus`, ...).
#' @param outdir Output directory for HTML/PDF (created if missing).
#' @param treatment_1_res Optional results table for a second contrast to compute
#'   treatment_1lustat direction per assay.
#' @param region_cols Character vector of region column names (defaults shown).
#' @param assay_col, panel_col Column names for assay and panel.
#' @param padj_col, est_col Column names for adjusted p-value and estimate.
#' @param brain_detect_col Column indicating brain detection (0/1).
#' @param padj_thr Numeric FDR threshold for direction calls (default 0.10).
#' @param file_prefix File name prefix (without extension).
#' @return Invisible Plotly widget (the interactive heatmap); files are written
#'   to `outdir`.
#' @export
annotations_heatmap <- function(
  res,
  outdir,
  treatment_1_res = NULL,
  region_cols = c("HigherInCaudate","HigherInHippocampus","HigherInCerebralCortex","HigherInCerebellum"),
  assay_col = "Assay",
  panel_col = "Panel",
  padj_col = "Adjusted_pval",
  est_col  = "estimate",
  brain_detect_col = "BrainDetected",
  padj_thr = 0.10,
  file_prefix = "brain_annotations_heatmap"
) {
  stopifnot(all(region_cols %in% names(res)))
  if (!base::dir.exists(outdir)) base::dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Prepare Treatment directions on the provided 'res'
  res_disease <- res %>%
    dplyr::mutate(
      disease_dir = .dir_from(!!rlang::sym(padj_col), !!rlang::sym(est_col), thr = padj_thr),
      brain_det = dplyr::case_when(
        !!rlang::sym(brain_detect_col) %in% 1 ~ "detected",
        !!rlang::sym(brain_detect_col) %in% 0 ~ "not_detected",
        TRUE ~ "NA"
      )
    )

  # Optional: compute treatment_1 direction by Assay if treatment_1_res provided
  treatment_1_dirs <- NULL
  if (!base::is.null(treatment_1_res)) {
    treatment_1_dirs <- treatment_1_res %>%
      dplyr::select(dplyr::all_of(c(assay_col, padj_col, est_col))) %>%
      dplyr::mutate(treatment_1_dir = .dir_from(!!rlang::sym(padj_col), !!rlang::sym(est_col), thr = padj_thr)) %>%
      dplyr::select(dplyr::all_of(c(assay_col, "treatment_1_dir"))) %>%
      dplyr::distinct()
  }

  # Collapse to unique assays, aggregating region flags across duplicates
  agg <- res_disease %>%
    dplyr::group_by(!!rlang::sym(assay_col)) %>%
    dplyr::summarize(
      !!rlang::sym(panel_col) := dplyr::first(!!rlang::sym(panel_col)),
      disease_dir = dplyr::first(disease_dir),
      brain_det = dplyr::first(base::factor(brain_det, levels = c("detected","not_detected","NA"))),
      dplyr::across(
        dplyr::all_of(region_cols),
        ~ as.integer(any(.x %in% c(1L, TRUE), na.rm = TRUE))
      ),
      .groups = "drop"
    ) %>%
    {
      if (!base::is.null(treatment_1_dirs)) dplyr::left_join(., treatment_1_dirs, by = assay_col)
      else dplyr::mutate(., treatment_1_dir = base::factor("no_change", levels = c("up","down","no_change")))
    }

  # Matrix for heatmap: rows = Assay, cols = 4 regions
  mat <- agg %>%
    dplyr::select(dplyr::all_of(region_cols)) %>%
    base::as.matrix()
  base::rownames(mat) <- agg[[assay_col]]

  # Row/Col annotations
  row_anno <- base::data.frame(
    panel = agg[[panel_col]],
    disease = agg$disease_dir,
    treatment_1 = agg$treatment_1_dir,
    brain = agg$brain_det,
    row.names = base::rownames(mat),
    check.names = FALSE
  )
  col_anno <- base::data.frame(
    region = base::factor(base::colnames(mat), levels = region_cols),
    row.names = base::colnames(mat),
    check.names = FALSE
  )

  # Palettes for annotations
  panel_lvls <- base::unique(row_anno$panel)
  pan_cols <- stats::setNames(colorspace::qualitative_hcl(base::length(panel_lvls), palette = "Dark 3"), panel_lvls)
  dir_cols  <- c(up = "#E41A1C", down = "#377EB8", no_change = "#BDBDBD")
  brain_cols <- c(detected = "#4DAF4A", not_detected = "#F0F0F0", "NA" = "#CCCCCC")
  region_cols_pal <- stats::setNames(colorspace::qualitative_hcl(base::length(region_cols), palette = "Set2"), region_cols)

  # Interactive (heatmaply)
  widget <- heatmaply::heatmaply(
    x = mat,
    na.value = "grey85",
    colors = viridisLite::viridis,
    row_side_colors = row_anno,
    col_side_colors = col_anno,
    row_side_palette = base::list(
      panel = pan_cols,
      disease = dir_cols,
      treatment_1 = dir_cols,
      brain = brain_cols
    ),
    col_side_palette = base::list(
      region = region_cols_pal
    ),
    showticklabels = c(TRUE, TRUE),
    xlab = NULL, ylab = NULL,
    plot_method = "plotly",
    k_col = 2,
    k_row = NA
  )

  # Save HTML
  html_path <- base::file.path(outdir, base::paste0(file_prefix, ".html"))
  base::suppressWarnings(htmlwidgets::saveWidget(widget, file = html_path, selfcontained = TRUE))

  # Static (PDF via pheatmap)
  ann_colors <- base::list(
    panel = pan_cols,
    disease = dir_cols,
    treatment_1 = dir_cols,
    brain = brain_cols,
    region = region_cols_pal
  )
  pdf_path <- base::file.path(outdir, base::paste0(file_prefix, ".pdf"))
  grDevices::pdf(pdf_path, width = 10, height = 12)
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = viridisLite::viridis(256),
    na_col = "#E5E5E5",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_row = row_anno,
    annotation_col = col_anno,
    annotation_colors = ann_colors,
    fontsize = 8, border_color = NA
  )
  grDevices::dev.off()

  base::message("Saved: ", pdf_path, " and ", html_path)
  base::invisible(widget)
}
