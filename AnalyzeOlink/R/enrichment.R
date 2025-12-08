#' Run Functional Enrichment
#'
#' @description
#' Performs functional enrichment (ORA, GSEA) for each contrast, method, and
#' test type. Only runs enrichment if at least `cutoff` assays are marked
#' `"Significant"`. Results are split by ontology (CC, BP, MF, other).
#'
#' @param res_list List of test results from `run_tests()` and
#'   `add_estimate_column()`.
#' @param config List with `enrichment$methods`, `enrichment$organism`,
#'   `enrichment$ontology`, and `alpha`.
#' @param cutoff Minimum number of significant assays required (default 5).
#'
#' @return Nested list `enrich_list[[contrast]][[method]][[test]][[direction]]`,
#'   where each direction is a list of enrichment tibbles by ontology.
#' @export
enrich <- function(res_list, config, cutoff = 5) {
  requireNamespace("dplyr"); requireNamespace("stringr")
  methods  <- config$enrichment$methods
  organism <- config$enrichment$organism
  ontology <- config$enrichment$ontology

  if ("ORA" %in% methods && !requireNamespace("msigdbdf", quietly = TRUE)) {
    install.packages(
      "msigdbdf",
      repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org")
    )
  }

  enrich_list <- list()
  directions  <- c("up", "down", "both")
  ontos       <- c("CC", "BP", "MF", "other")

  for (contrast in names(res_list)) {
    rd <- res_list[[contrast]]
    if (!is.list(rd)) next
    long_data <- rd$params$df_work

    enrich_list[[contrast]] <- list()
    for (method in methods) {
      enrich_list[[contrast]][[method]] <- list()

      for (test_type in setdiff(names(rd), c("params","emm_int"))) {
        test_res <- rd[[test_type]]
        if (!"estimate" %in% names(test_res)) {
          stop("No estimate column, can't enrich results for: ", contrast)
        }

        enrich_list[[contrast]][[method]][[test_type]] <- list()

        for (dir in directions) {
          tr <- switch(dir,
            up   = dplyr::filter(test_res, Threshold=="Significant", estimate >  0),
            down = dplyr::filter(test_res, Threshold=="Significant", estimate <  0),
            both = dplyr::filter(test_res, Threshold=="Significant")
          )

          if (nrow(tr) < cutoff) {
            enrich_list[[contrast]][[method]][[test_type]][[dir]] <-
              setNames(vector("list", length(ontos)), ontos)
            next
          }

          enriched <- tryCatch(
            OlinkAnalyze::olink_pathway_enrichment(
              data          = long_data,
              test_results  = tr,
              method        = method,
              ontology      = ontology,
              organism      = organism,
              pvalue_cutoff = config$alpha
            ),
            error = function(e) {
              warning("Enrichment error for ", contrast, "/", test_type, "/", dir, ": ", e$message)
              NULL
            }
          )
          if (is.null(enriched) || nrow(enriched)==0) {
            enrich_list[[contrast]][[method]][[test_type]][[dir]] <-
              setNames(vector("list", length(ontos)), ontos)
            next
          }

          enriched <- enriched %>%
            dplyr::arrange(p.adjust) %>%
            dplyr::mutate(
              Description_short = dplyr::if_else(
                stringr::str_length(Description) <= 25,
                Description,
                paste0(stringr::str_sub(Description,1,25),"…")
              ),
              ID        = paste0(dplyr::row_number(),"_",Description_short),
              Threshold = if_else(p.adjust < config$alpha, "Significant","Non-significant")
            ) %>%
            dplyr::select(ID, dplyr::everything(), -Description_short)

          CC    <- dplyr::filter(enriched,
                    stringr::str_starts(Description,"GOCC_") %>%
                    stringr::str_detect(ID,"^\\d+_GOCC_")) %>%
                    tibble::as_tibble()
          BP    <- dplyr::filter(enriched,
                    stringr::str_starts(Description,"GOBP_") %>%
                    stringr::str_detect(ID,"^\\d+_GOBP_")) %>%
                    tibble::as_tibble()
          MF    <- dplyr::filter(enriched,
                    stringr::str_starts(Description,"GOMF_") %>%
                    stringr::str_detect(ID,"^\\d+_GOMF_")) %>%
                    tibble::as_tibble()
          other <- dplyr::filter(enriched,
                    ! (stringr::str_starts(Description,"GOCC_|GOBP_|GOMF_") %>%
                       stringr::str_detect(ID,"^\\d+_(GOCC|GOBP|GOMF)_"))) %>%
                       tibble::as_tibble()

          enrich_list[[contrast]][[method]][[test_type]][[dir]] <- list(
            CC    = CC,
            BP    = BP,
            MF    = MF,
            other = other
          )
        }
      }
    }
  }

  enrich_list
}


#' Run Functional Enrichment in Parallel
#'
#' @description
#' Runs `enrich()` across results in parallel if `config$parallel` is TRUE,
#' dispatching each element of `res_list` to its own core. Otherwise runs in series.
#'
#' @param res_list List of result objects from `run_tests()`.
#' @param config List with fields:
#'   * `parallel`: logical, run in parallel if TRUE.
#'   * `cores`: integer, number of cores to use.
#'   * other fields used by `enrich()`.
#'
#' @return Named list of enrichment results matching `names(res_list)`.
#' @export
parallel_enrich <- function(res_list, config) {
  if (isTRUE(config$parallel)) {
    ncores <- config$cores

    .run_one_enricher <- function(name) {
      single_input <- res_list[name]
      out <- enrich(single_input, config)
      out[[1]]
    }

    enr_list <- parallel::mclapply(
      X        = names(res_list),
      FUN      = .run_one_enricher,
      mc.cores = ncores
    )
    names(enr_list) <- names(res_list)
    return(enr_list)
  } else {
    return(enrich(res_list, config))
  }
}


#' Format Enrichment Results
#'
#' @description
#' Standardizes and formats enrichment results by:
#' * Dropping all-NA columns.
#' * Reordering core columns.
#' * Optionally dropping helper columns.
#' * Formatting numeric columns.
#'
#' @param df Enrichment results tibble.
#' @param drop_enrich_cols Logical; drop enrichment-specific helper columns.
#' @param scientific_notation Logical; format p-values in scientific notation.
#'
#' @return Tibble of formatted enrichment results.
#' @export
format_enrich_results <- function(df,
                                  drop_enrich_cols = FALSE,
                                  scientific_notation = FALSE) {
  df <- df %>% dplyr::select(where(~ !all(is.na(.x))))

  formatted <- as.data.frame(df) %>%
    dplyr::select(
      any_of(c(
        "ID", "Description", "GeneRatio", "BgRatio",
        "Count", "geneID", "leading_edge", "core_enrichment",
        "pvalue", "p.adjust", "qvalue", "Threshold"
      )),
      everything()
    )

  if (isTRUE(drop_enrich_cols)) {
    formatted <- formatted %>%
      dplyr::select(-any_of(c("ID", "leading_edge", "core_enrichment", "geneID")))
  }

  if (isTRUE(scientific_notation)) {
    formatted <- formatted %>%
      dplyr::mutate(
        dplyr::across(
          c(pvalue, p.adjust, qvalue),
          ~ sprintf("%.2e", .)
        )
      ) %>%
      dplyr::select(-any_of("rank", "qvalue"))
  } else {
    formatted <- formatted %>%
      dplyr::mutate(
        dplyr::across(
          where(is.numeric),
          ~ round(., 3)
        )
      )
  }

  tibble::as_tibble(formatted)
}
