#' @title Post-hoc emmeans handler
#' @description
#' Internal helper for handling emmeans post-hoc results. Supports both
#' simple main-effect contrasts and covariate interactions.
#'
#' @param anova_tbl ANOVA table tibble (already constructed).
#' @param fit A fitted model object (lmrob, rlmer, afex model, etc.).
#' @param emm A list of emmeans specifications (simple contrasts or interactions).
#' @param variable Name of the grouping variable (character).
#' @param alpha Numeric significance level for confidence intervals.
#' @param key A one-row tibble of identifying columns (Assay, OlinkID, UniProt, Panel).
#'
#' @return A list with components:
#' \describe{
#'   \item{anova_tbl}{Possibly updated ANOVA table tibble.}
#'   \item{emm_int}{Interaction tibble with trends and contrasts (or empty).}
#'   \item{repl_tbl}{Replacement tibble for simple main-effect rows (or empty).}
#' }
#' @noRd
posthoc_emm_handler <- function(anova_tbl, fit, emm, variable, alpha, key) {
  # emmeans post-hoc handling
  #   - simple main effect           (1 row -> anova_tbl)
  #   - interaction Int              (separate emm_int table)
  emm_int  <- tibble::tibble()
  repl_tbl <- tibble::tibble()

    ## make sure columns exist so coalesce() works later
    if (!"estimate" %in% names(anova_tbl))
      anova_tbl <- dplyr::mutate(anova_tbl, estimate = NA_real_)
    if (!"df" %in% names(anova_tbl))
      anova_tbl <- dplyr::mutate(anova_tbl, df = NA_real_)

    ## ---------- A) interaction branch (build emm_int) ---------------
    if ("Int" %in% names(emm)) {

      emm_int <- purrr::map_dfr(
        emm["Int"],                                  # keeps the key "Int"
        \(spec, nm = "Int") {
          res <- interaction_emm(fit, spec)

          tr <- tibble::as_tibble(res$trends) %>%
                  dplyr::mutate(name = nm, component = "trend")

          ct <- tibble::as_tibble(res$contrasts) %>%
                  dplyr::mutate(name = nm, component = "contrast",
                                Treatment = NA, Age_c.trend = NA_real_) %>%
                  dplyr::relocate(Treatment, Age_c.trend, .after = component)

          dplyr::bind_rows(tr, ct)
        }
      ) %>%
        dplyr::bind_cols(key) %>%
        dplyr::relocate(component, name, .before = Treatment)
    }

    ## ---------- B) simple main-effect contrast ----------------------
    simple_names <- setdiff(names(emm), "Int")   # usually one entry

    if (length(simple_names)) {

      nm   <- simple_names[1]                    # use the first (or only) one
      spec <- emm[[nm]]

      ## build contrast ----------------------------------------------
      if (is.numeric(spec)) {
        ctr <- emmeans::contrast(
                 emmeans::emmeans(fit, specs = variable),
                 method = setNames(list(spec), nm),
                 adjust = "none")
      } else {                                   # length-3 pairwise
        lvl1 <- spec[1]; lvl2 <- spec[2]
        ctr  <- emmeans::contrast(
                 emmeans::emmeans(fit, specs = nm),
                 method  = "pairwise",
                 levels  = c(lvl1, lvl2),
                 adjust  = "none")
      }

      s <- summary(ctr, infer = TRUE, level = 1 - alpha)
    }

    repl_tbl <- tibble::tibble(
        Effect   = nm,                # Treatment
        estimate = s$estimate,
        df       = s$df,
        p.value  = s$p.value,
        t.ratio  = s$t.ratio
      )

  # Return the updated anova_tbl, plus any derived tables
  list(
    anova_tbl = anova_tbl,
    emm_int   = emm_int,
    repl_tbl  = repl_tbl
  )
}


#' Simple Wald χ² ANOVA for robustlmm::rlmer models
#'
#' Perform a Wald‐type χ² analysis of variance on a `rlmerMod` object from the
#' robustlmm package.
#'
#' @param model An object of class \code{rlmerMod}, as returned by
#'   \code{robustlmm::rlmer()}. Must inherit from \code{"rlmerMod"}.
#'
#' @return A \code{tibble} with the following columns:
#' \itemize{
#'   \item \code{Effect}: Term names with trailing \dQuote{1} (or
#'     \dQuote{1} before a colon) removed.
#'   \item \code{df}: Degrees of freedom for each test (always 1).
#'   \item \code{Chisq}: Wald χ² statistic (square of the \code{t} value).
#'   \item \code{p.value}: P‐value from the χ² distribution with 1 df.
#'   \item \code{estimate}: Estimated coefficient for each term.
#' }
#'
#' @examples
#' \dontrun{
#' library(robustlmm)
#' fit <- rlmer(y ~ x + (1 | g), data = df)
#' wald_anova_rlmer(fit)
#' }
#'
#' @importFrom tibble tibble
#' @importFrom stats pchisq
#' @export
wald_anova_rlmer <- function(model) {
  stopifnot(inherits(model, "rlmerMod"))

  cm <- summary(model)$coefficients   # matrix with Estimate, SE, t value

  # remove the coding “1” that lme4/robustlmm appends to factor names
  # delete a trailing 1 ...or.. a 1 immediately in front of “:”
  clean_eff <- gsub("1(?=:|$)", "", rownames(cm), perl = TRUE)

  tibble::tibble(
    Effect   = clean_eff,
    df       = 1L,
    Chisq    = cm[, "t value"]^2,     # Wald statistic
    p.value  = stats::pchisq(Chisq, df = 1L, lower.tail = FALSE),
    estimate = cm[, "Estimate"]
  )
}


#' Robust ANCOVA (robustbase::lmrob) analysis with diagnostics and post-hoc tests
#'
#' Fits a robust ANCOVA model via \code{lmrob}, performs a Type-III ANOVA
#' with partial η², adds residual diagnostics, information criteria, and
#' optional raw log2-fold-change. Supports post-hoc contrasts via \pkg{emmeans}.
#'
#' @inheritParams one_mixed
#' @return A list with:
#' \itemize{
#'   \item \code{anova}: ANOVA table with diagnostics, \code{pes}, optional \code{log2FC}.
#'   \item \code{emm_int}: Interaction post-hoc table (empty if none).
#' }
#' @export
one_robustbase <- function(
  df_sub, key, formula, form_obj, observed,
  variable, alt_lfc_col, report_lfc, emm,
  alpha, outcome, treatment, control
) {

  # 1) Fit model & make a Type-III table
  # ----- robust ANCOVA ----
  fit <- robustbase::lmrob(
    formula = form_obj,
    data    = df_sub,
    setting = "KS2014"
  )

  anova_raw <- car::Anova(fit, type = 3, test.statistic = "F", singular.ok = TRUE)        # S4 object
  anova_src <- as.data.frame(as.matrix(anova_raw)) %>%
    tibble::rownames_to_column("Effect") %>%
    dplyr::rename(df = Df)

  pcol <- grep("^Pr\\(>", names(anova_src), value = TRUE)
  anova_src <- dplyr::rename(anova_src,
                             p.value = dplyr::all_of(pcol)) %>%
    dplyr::filter(Effect != "Residuals")

  # Add estimate for all df = 1 effects
  ## 1. get coefficients and the ‘assign’ attribute
  beta        <- stats::coef(fit)
  mm          <- stats::model.matrix(fit)
  term_id     <- attr(mm, "assign")           # one integer per model column
  coef_map    <- split(seq_along(beta), term_id)

  ## 2. pick terms that have exactly *one* column -> 1-df in the ANOVA table
  # keep only real terms (assign > 0), no intercept
  one_df_idx  <- names(coef_map)[vapply(coef_map, length, 1L) == 1L &
                              names(coef_map) != "0"]

  coef_tbl <- tibble::tibble(
    Effect   = attr(stats::terms(fit), "term.labels")[as.integer(one_df_idx)],
    estimate = beta[unlist(coef_map[one_df_idx])]
  )

  ## 3. join – every 1-df (two-level factor) Effect gets its coefficient, everything else stays NA
  anova_src <- anova_src %>%
                 dplyr::left_join(coef_tbl, by = "Effect") %>%
                 dplyr::mutate(
                   estimate = dplyr::coalesce(estimate, NA_real_)  # keep column consistent
                 )

  anova_src <- dplyr::filter(anova_src, Effect != "(Intercept)")

  # 2) Add partial η²
  anova_src <- add_pes_rob(anova_src, model = fit)

  # 3) Information criteria
  res_vec  <- stats::residuals(fit)

  anova_tbl <- anova_src |>
    dplyr::mutate(
      resids = list(setNames(res_vec, df_sub$SampleID)),
    )

  # 4) Normality / variance diagnostics
  tmp <- data.frame(
           resid = res_vec,
           grp   = df_sub[[variable]]
         ) |>
         dplyr::filter(grp %in% c(treatment, control)) |>
         droplevels()

  sw_p <- tryCatch(stats::shapiro.test(tmp$resid)$p.value,
                   error = \(e) NA_real_)
  ad_p <- tryCatch(nortest::ad.test(tmp$resid)$p.value,
                   error = \(e) NA_real_)
  lt_p <- tryCatch(car::leveneTest(resid ~ grp, data = tmp,
                                   center = median)[1, "Pr(>F)"] |>
                        as.numeric(),
                   error = \(e) NA_real_)
  ft_p <- tryCatch(stats::var.test(resid ~ grp, data = tmp)$p.value |>
                        as.numeric(),
                   error = \(e) NA_real_)

  anova_tbl <- anova_tbl |>
    dplyr::mutate(
      SW_pval = sw_p,
      AD_pval = ad_p,
      LT_pval = lt_p,
      FT_pval = ft_p
    ) |>
    dplyr::bind_cols(key)

  # 5) Optional raw log2-fold-change
  if (isTRUE(report_lfc)) {
    lfc_col <- if (!is.null(alt_lfc_col) && alt_lfc_col %in% names(df_sub))
                 alt_lfc_col else variable
    log2fc <- if (
      lfc_col %in% names(df_sub) &&
      is.factor(df_sub[[lfc_col]]) &&
      nlevels(df_sub[[lfc_col]]) == 2
    ) {
      lv <- levels(df_sub[[lfc_col]])
      mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[1]], na.rm = TRUE) -
        mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[2]], na.rm = TRUE)
    } else {
      NA_real_
    }
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = log2fc)
  } else {
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = NA_real_)
  }

  # emmeans post-hoc handling
  #   - simple main effect           (1 row -> anova_tbl)
  #   - interaction Int              (separate emm_int table)
  emm_int  <- tibble::tibble()
  repl_tbl <- tibble::tibble()

  if (is.list(emm)) {
    emm_res   <- posthoc_emm_handler(anova_tbl, fit, emm, variable, alpha, key)
    anova_tbl <- emm_res$anova_tbl
    emm_int   <- emm_res$emm_int
    repl_tbl  <- emm_res$repl_tbl
  }

  ##  splice-in the emmeans main-effect row (lmrob only)
  if (nrow(repl_tbl)) {                # <- repl_tbl only exists for lmrob

    anova_tbl <- anova_tbl %>%
      dplyr::left_join(repl_tbl, by = "Effect", suffix = c("", ".emm")) %>%

      # overwrite for that ONE row
      dplyr::mutate(
        estimate = dplyr::coalesce(estimate.emm, estimate),
        df       = dplyr::coalesce(df.emm,       df),
        p.value  = dplyr::coalesce(p.value.emm,  p.value),

        # keep the F column except on the overridden row
        F        = dplyr::case_when(
                     !is.na(p.value.emm) ~ NA_real_,   # we used emmeans here
                     TRUE                ~ F           # keep original F
                   ),

        # create t.ratio only where we have the emmeans row
        t.ratio  = dplyr::case_when(
                     !is.na(p.value.emm) ~ t.ratio,
                     TRUE                ~ NA_real_
                   )
      ) %>%
      dplyr::select(-dplyr::ends_with(".emm"))
  }

  ## -----------------------------------------------------------
  ##  result anova_tbl:
  ##  * lmrob: columns F (+ t.ratio only on the main effect row)
  ##  * rlmer:  column Chisq, never t.ratio
  ## -----------------------------------------------------------
  list(anova = anova_tbl, emm_int = emm_int)
}


#' Robust linear mixed model (robustlmm::rlmer) with diagnostics and post-hoc tests
#'
#' Fits a robust mixed-effects model via \code{rlmer}, performs a Type-III
#' Wald χ² ANOVA with partial η², adds residual diagnostics and optional
#' raw log2-fold-change, plus post-hoc contrasts.
#'
#' @inheritParams one_mixed
#' @return A list with:
#' \itemize{
#'   \item \code{anova}: ANOVA table with diagnostics, \code{pes}, optional \code{log2FC}.
#'   \item \code{emm_int}: Interaction post-hoc table.
#' }
#' @export
one_robustlmm <- function(
  df_sub, key, formula, form_obj, observed,
  variable, alt_lfc_col, report_lfc, emm,
  alpha, outcome, treatment, control
) {

  # 1) Fit model & make a Type-III table
  # ----- robust LME ------
  fit <- robustlmm::rlmer(
    formula = form_obj,
    data    = df_sub
  )

  ## -- Type-III Wald test (rlmer) --
  anova_src <- wald_anova_rlmer(fit)
  anova_src <- dplyr::filter(anova_src, Effect != "(Intercept)")

  # 2) Add partial η²
  anova_src <- add_pes_rob(anova_src, model = fit)

  # 3) Information criteria
  res_vec  <- stats::residuals(fit)

  anova_tbl <- anova_src |>
    dplyr::mutate(
      resids = list(setNames(res_vec, df_sub$SampleID)),
    )

  # 4) Normality / variance diagnostics
  tmp <- data.frame(
           resid = res_vec,
           grp   = df_sub[[variable]]
         ) |>
         dplyr::filter(grp %in% c(treatment, control)) |>
         droplevels()

  sw_p <- tryCatch(stats::shapiro.test(tmp$resid)$p.value,
                   error = \(e) NA_real_)
  ad_p <- tryCatch(nortest::ad.test(tmp$resid)$p.value,
                   error = \(e) NA_real_)
  lt_p <- tryCatch(car::leveneTest(resid ~ grp, data = tmp,
                                   center = median)[1, "Pr(>F)"] |>
                        as.numeric(),
                   error = \(e) NA_real_)
  ft_p <- tryCatch(stats::var.test(resid ~ grp, data = tmp)$p.value |>
                        as.numeric(),
                   error = \(e) NA_real_)

  anova_tbl <- anova_tbl |>
    dplyr::mutate(
      SW_pval = sw_p,
      AD_pval = ad_p,
      LT_pval = lt_p,
      FT_pval = ft_p
    ) |>
    dplyr::bind_cols(key)

  # 5) Optional raw log2-fold-change
  if (isTRUE(report_lfc)) {
    lfc_col <- if (!is.null(alt_lfc_col) && alt_lfc_col %in% names(df_sub))
                 alt_lfc_col else variable
    log2fc <- if (
      lfc_col %in% names(df_sub) &&
      is.factor(df_sub[[lfc_col]]) &&
      nlevels(df_sub[[lfc_col]]) == 2
    ) {
      lv <- levels(df_sub[[lfc_col]])
      mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[1]], na.rm = TRUE) -
        mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[2]], na.rm = TRUE)
    } else {
      NA_real_
    }
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = log2fc)
  } else {
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = NA_real_)
  }

  # emmeans post-hoc handling
  #   - simple main effect           (1 row -> anova_tbl)
  #   - interaction Int              (separate emm_int table)
  emm_int  <- tibble::tibble()
  repl_tbl <- tibble::tibble()

  if (is.list(emm)) {
    emm_res   <- posthoc_emm_handler(anova_tbl, fit, emm, variable, alpha, key)
    anova_tbl <- emm_res$anova_tbl
    emm_int   <- emm_res$emm_int
    repl_tbl  <- emm_res$repl_tbl
  }

  ##  splice-in the emmeans main-effect row (lmrob only)
  if (nrow(repl_tbl)) {                # <- repl_tbl only exists for lmrob

    anova_tbl <- anova_tbl %>%
      dplyr::left_join(repl_tbl, by = "Effect", suffix = c("", ".emm")) %>%

      # overwrite for that ONE row
      dplyr::mutate(
        estimate = dplyr::coalesce(estimate.emm, estimate),
        df       = dplyr::coalesce(df.emm,       df),
        p.value  = dplyr::coalesce(p.value.emm,  p.value),

        # keep the F column except on the overridden row
        F        = dplyr::case_when(
                     !is.na(p.value.emm) ~ NA_real_,   # we used emmeans here
                     TRUE                ~ F           # keep original F
                   ),

        # create t.ratio only where we have the emmeans row
        t.ratio  = dplyr::case_when(
                     !is.na(p.value.emm) ~ t.ratio,
                     TRUE                ~ NA_real_
                   )
      ) %>%
      dplyr::select(-dplyr::ends_with(".emm"))
  }

  ## -----------------------------------------------------------
  ##  result anova_tbl:
  ##  * lmrob: columns F (+ t.ratio only on the main effect row)
  ##  * rlmer:  column Chisq, never t.ratio
  ## -----------------------------------------------------------
  list(anova = anova_tbl, emm_int = emm_int)
}


#' Robust ANCOVA / robust LME analysis wrapper
#'
#' Run \code{one_robustbase()} or \code{one_robustlmm()} across all protein-panels,
#' apply multiple-testing correction, annotate significance thresholds, and
#' return combined ANOVA and post-hoc results.
#'
#' @inheritParams robust_test
#' @return A list with:
#' \itemize{
#'   \item \code{anova}: Combined ANOVA table.
#'   \item \code{emm_int}: Combined interaction table.
#' }
#' @export
robust_test <- function(
  long_df, formula, treatment, control,
  variable = NULL, observed = NULL, alpha = 0.10,
  report_lfc = FALSE, emm = NULL, alt_lfc_col = NULL,
  lmer = FALSE, rm_incomplete = TRUE
) {
  if (missing(long_df) || missing(formula)) {
    stop("Both 'long_df' and 'formula' must be supplied.")
  }

  form_obj <- stats::as.formula(formula)
  outcome  <- all.vars(form_obj)[1]
  vars_nom <- all.vars(form_obj)
  extras   <- c("SampleID", "Assay", "OlinkID", "UniProt", "Panel")
  data_sel <- dplyr::select(long_df,
                            dplyr::all_of(c(vars_nom, extras, alt_lfc_col)))

  ## 1) optional NA pruning
  if (isTRUE(rm_incomplete)) {
    preds    <- setdiff(vars_nom, outcome)
    cmpl_idx <- stats::complete.cases(data_sel[preds])
    if (any(!cmpl_idx)) {
      warning("Removed ", sum(!cmpl_idx), " incomplete rows.")
      data_sel <- data_sel[cmpl_idx, ]
    }
  }

  ## 2) per-protein analysis
  if (!isTRUE(lmer)) {
    results_list <- data_sel %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::group_map(
        \(df_sub, key)
          one_robustbase(
            df_sub, key, formula, form_obj, observed, variable,
            alt_lfc_col, report_lfc, emm, alpha, outcome,
            treatment, control
          ),
        .keep = TRUE
      )
  } else {
    ## parallel branch for robust LME
    grp  <- data_sel %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel)
    keys <- grp %>% dplyr::group_keys()
    dfs  <- grp %>% dplyr::group_split()

    params <- purrr::map2(
      dfs,
      split(keys, seq(nrow(keys))),
      \(df, key_row) list(df_sub = df, key = key_row)
    )

    # 100 proteins per chunk
    chunks <- split(params, ceiling(seq_along(params) / 400))
    # parallel::detectCores() will grab every core on your node restict the
    # number of cores used for every batch of 200 proteins to 8 Since 8 cores
    # * 400 proteins = 3200 proteins run simultaneously (>2888).
    ncores <- 8

    results_list <- parallel::mclapply(
      X        = chunks,
      FUN      = function(chunk) {
        lapply(chunk, function(p)
          one_robustlmm(
            df_sub      = p$df_sub,
            key         = p$key,
            formula     = formula,
            form_obj    = form_obj,
            observed    = observed,
            variable    = variable,
            alt_lfc_col = alt_lfc_col,
            report_lfc  = report_lfc,
            emm         = emm,
            alpha       = alpha,
            outcome     = outcome,
            treatment   = treatment,
            control     = control
          )
        )
      },
      mc.cores = ncores
    )
    # flatten nested lists
    results_list <- unlist(results_list, recursive = FALSE)
  }

  anova_raw   <- purrr::map(results_list, "anova")   %>% dplyr::bind_rows()
  emm_int_tbl <- purrr::map(results_list, "emm_int") %>%
                 purrr::compact() %>%
                 dplyr::bind_rows()

  ## 3) multiple-testing adjustment, append raw npx data, and sort
  anova_final <- anova_raw %>%
    dplyr::group_by(Effect) %>%
    dplyr::mutate(
      Adjusted_pval = stats::p.adjust(p.value, method = "BH"),
      SW_padj       = stats::p.adjust(SW_pval, method = "BH"),
      AD_padj       = stats::p.adjust(AD_pval, method = "BH"),
      LT_padj       = stats::p.adjust(LT_pval, method = "BH"),
      FT_padj       = stats::p.adjust(FT_pval, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Threshold = ifelse(Adjusted_pval < alpha,
                         "Significant", "Non-significant"),
      formula   = formula
    ) %>%
    add_raw_npx(long_df) %>%
    sort_results()

  list(anova = anova_final, emm_int = emm_int_tbl)
}


#' Fit a BRMS contrast per protein or across proteins
#'
#' Fit Bayesian Student-t mixed models with \pkg{brms}, compute posterior
#' contrasts and p-values, add log2FC, adjust p-values, and annotate significance.
#'
#' @param df Data frame in long format.
#' @param formula A brms formula.
#' @param priors Priors list (optional).
#' @param by_protein Logical; if TRUE fit per protein, else hierarchical.
#' @param alpha Significance threshold.
#' @param ... Extra args to \code{brm()}.
#' @return A list with elements:
#' \itemize{
#'   \item \code{res}: Tibble of per-protein results.
#'   \item \code{fit}: The fitted brms model(s).
#' }
#' @export
brms_contrast <- function(df,
                          formula,
                          priors     = NULL,
                          by_protein = FALSE,
                          alpha      = 0.10,
                          ...) {

  # Suppress rstan browseURL
  old_browser <- getOption("browser")
  options(browser = function(...) invisible(TRUE))
  on.exit(options(browser = old_browser), add = TRUE)

  # Ensure protein key
  if (!"protein" %in% names(df)) {
    df <- dplyr::mutate(
      df,
      protein = interaction(Assay, Panel, OlinkID, UniProt, sep = "|", drop = TRUE)
    )
  }

  outcome   <- all.vars(formula)[1]
  variable  <- all.vars(formula)[2]
  form_chr  <- paste(deparse(formula), collapse = " ")

  raw_fc <- function(dat) {
    lv <- levels(dat[[variable]])
    mean(dat[[outcome]][dat[[variable]] == lv[1]], na.rm = TRUE) -
      mean(dat[[outcome]][dat[[variable]] == lv[2]], na.rm = TRUE)
  }

  # Per-protein mode
  if (by_protein) {
    res <- df |>
      dplyr::group_by(Assay, Panel, OlinkID, UniProt, protein) |>
      dplyr::group_map(~{
        dat <- .x
        fit <- brms::brm(
          formula = formula,
          data    = dat,
          family  = brms::student(),
          prior   = priors,
          silent  = TRUE, refresh = 0, ...
        )

        lvl     <- levels(dat[[variable]])
        coef_n  <- paste0("b_", variable, lvl[2])
        posts   <- posterior::as_draws_df(fit, variables = coef_n)[[coef_n]]
        pval    <- 2 * min(mean(posts > 0), mean(posts < 0))

        est <- emmeans::emmeans(fit, specs = variable) |>
               emmeans::contrast(method = "pairwise", levels = lvl) |>
               summary(infer = FALSE) |>
               dplyr::pull(estimate)

        tibble::tibble(
          Assay    = .y$Assay, Panel = .y$Panel,
          OlinkID  = .y$OlinkID, UniProt = .y$UniProt,
          p.value  = pval,
          estimate = est,
          log2FC   = raw_fc(dat)
        )
      }, .keep = TRUE) |>
      dplyr::bind_rows()

  # Hierarchical mode
  } else {
    fit  <- brms::brm(
      formula = formula,
      data    = df,
      family  = brms::student(),
      prior   = priors,
      silent  = TRUE, refresh = 0, ...
    )

  # Set-up
  lvl      <- levels(df[[variable]])                  # c("control","disease")   etc.
  draws    <- posterior::as_draws_df(fit)
  fix_n    <- paste0("b_", variable, lvl[2])          # e.g. "b_Groupdisease"

  rand_cols <- grep(
    paste0("r_protein\\[.*", variable, lvl[2], "\\]"),
    colnames(draws), value = TRUE
  )
  prots <- sub("r_protein\\[([^,]+),.*", "\\1", rand_cols)

  # Summarise each protein in one shot
  stat_vec <- function(col) {
    tot <- draws[[fix_n]] + draws[[col]]          # posterior draws of θ_p
    pd  <- max(mean(tot > 0), mean(tot < 0))      # probability of direction

    c(                                                     # named numeric vector
      estimate = mean(tot),                               # == mu
      p.value  = 2 * (1 - pd),                            # 2 × min tail area
      cil      = quantile(tot, 0.025),
      cih      = quantile(tot, 0.975),
      pd       = pd
    )
  }

  stat_mat <- t(vapply(rand_cols, stat_vec,
                       FUN.VALUE = setNames(numeric(5),
                                            c("estimate","p.value",
                                              "cil","cih","pd"))))

  res <- tibble::as_tibble(stat_mat, .name_repair = "minimal") |>
    dplyr::mutate(protein = prots, .before = 1)

  # Add raw log2 fold-change and annotation columns
  res <- res |>
    dplyr::left_join(
      df |>
        dplyr::group_by(protein) |>
        dplyr::summarise(log2FC = raw_fc(cur_data()), .groups = "drop"),
      by = "protein"
    ) |>
    dplyr::left_join(
      dplyr::distinct(df, Assay, Panel, OlinkID, UniProt, protein),
      by = "protein"
    )
  }

  res <- dplyr::mutate(res,
    Adjusted_pval = stats::p.adjust(p.value, method = "BH"),
    Threshold = ifelse(Adjusted_pval < alpha, "Significant", "Non-significant"),
    formula = form_chr
  )

  # Optionally add raw NPX values to results
  if (all(c("SampleID", "NPX") %in% names(df))) {
    res <- add_raw_npx(res, df) %>% sort_results()
  } else {
    res <- sort_results(res)
  }

  list(res = res, fit = fit)
}


#' ANCOVA via afex::aov_car
#'
#' @inheritParams one_mixed
#' @return List with anova table and emm_int table.
#' @export
one_car <- function(
  df_sub, key, formula, form_obj, observed,
  variable, alt_lfc_col, report_lfc, emm,
  alpha, outcome, treatment, control
) {
  # 1) Fit model & extract base ANOVA table (aov_car)
  mdl <- suppressMessages(
    afex::aov_car(
      formula   = form_obj,
      data      = df_sub,
      factorize = FALSE,
      observed  = observed
    )
  )
  fit <- mdl$lm
  anova_src <- afex::nice(
    mdl, es = "pes", observed = observed,
    sig_symbols = rep("", 4), round_ps = identity
  ) |>
    tibble::as_tibble()

  # Add estimate for all df = 1 effects
  ## 1. get coefficients and the ‘assign’ attribute
  beta        <- stats::coef(fit)
  mm          <- stats::model.matrix(fit)
  term_id     <- attr(mm, "assign")           # one integer per model column
  coef_map    <- split(seq_along(beta), term_id)

  ## 2. pick terms that have exactly *one* column -> 1-df in the ANOVA table
  # keep only real terms (assign > 0), no intercept
  one_df_idx  <- names(coef_map)[vapply(coef_map, length, 1L) == 1L &
                              names(coef_map) != "0"]

  coef_tbl <- tibble::tibble(
    Effect   = attr(stats::terms(fit), "term.labels")[as.integer(one_df_idx)],
    estimate = beta[unlist(coef_map[one_df_idx])]
  )

  ## 3. join – every 1-df (two-level factor) Effect gets its coefficient, everything else stays NA
  anova_src <- anova_src %>%
                 dplyr::left_join(coef_tbl, by = "Effect") %>%
                 dplyr::mutate(
                   estimate = dplyr::coalesce(estimate, NA_real_)  # keep column consistent
                 )

  # 2) Add partial eta²
  anova_src <- add_pes(anova_src, model = fit)

   # 3) Diagnostics & ICs
  # Use only the rows actually used in the model (after NA dropping etc.)
  mf       <- stats::model.frame(fit)            # one row per observation used
  res_vec  <- stats::residuals(fit)
  idx      <- as.integer(rownames(mf))           # original row indices in df_sub
  res_ids  <- df_sub$SampleID[idx]

  n_obs    <- length(res_vec)
  p_parm   <- length(stats::coef(fit))
  aic_val  <- stats::AIC(fit)
  bic_val  <- stats::BIC(fit)
  aicc_val <- aic_val + (2 * p_parm * (p_parm + 1)) / (n_obs - p_parm - 1)

  anova_tbl <- anova_src |>
    dplyr::mutate(
      resids = list(stats::setNames(res_vec, res_ids)),
      AIC    = aic_val,
      BIC    = bic_val,
      AICc   = aicc_val
    )

  # 4) Normality & homogeneity tests
  # Take the grouping variable from the *model frame* to stay aligned with res_vec
  grp_full <- mf[[variable]]

  tmp <- data.frame(resid = res_vec, grp = grp_full) |>
    dplyr::filter(grp %in% c(treatment, control)) |>
    droplevels()

  sw_p <- tryCatch(stats::shapiro.test(tmp$resid)$p.value, error = \(e) NA_real_)
  ad_p <- tryCatch(nortest::ad.test(tmp$resid)$p.value, error = \(e) NA_real_)
  lt_p <- tryCatch(
    car::leveneTest(resid ~ grp, data = tmp, center = median)[1, "Pr(>F)"] |> as.numeric(),
    error = \(e) NA_real_
  )
  ft_p <- tryCatch(
    stats::var.test(resid ~ grp, data = tmp)$p.value |> as.numeric(),
    error = \(e) NA_real_
  )

  anova_tbl <- anova_tbl |>
    dplyr::mutate(
      SW_pval = sw_p,
      AD_pval = ad_p,
      LT_pval = lt_p,
      FT_pval = ft_p
    ) |>
    dplyr::bind_cols(key)

  # 5) Optional raw log2 fold-change
  if (isTRUE(report_lfc)) {
    lfc_col <- if (!is.null(alt_lfc_col) && alt_lfc_col %in% names(df_sub))
                 alt_lfc_col else variable
    log2fc <- if (
      lfc_col %in% names(df_sub) &&
      is.factor(df_sub[[lfc_col]]) &&
      nlevels(df_sub[[lfc_col]]) == 2
    ) {
      lv <- levels(df_sub[[lfc_col]])
      mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[1]], na.rm = TRUE) -
        mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[2]], na.rm = TRUE)
    } else {
      NA_real_
    }
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = log2fc)
  } else {
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = NA_real_)
  }

  # EMM post-hoc handling
  emm_int  <- tibble::tibble()
  repl_tbl <- tibble::tibble()
  if (is.list(emm)) {
    emm_res   <- posthoc_emm_handler(anova_tbl, fit, emm, variable, alpha, key)
    anova_tbl <- emm_res$anova_tbl
    emm_int   <- emm_res$emm_int
    repl_tbl  <- emm_res$repl_tbl
  }

  # 7) splice-in the emmeans main-effect row (lmrob only)
  if (nrow(repl_tbl)) {
    anova_tbl <- anova_tbl |>
      dplyr::left_join(repl_tbl, by = "Effect", suffix = c("", ".emm")) |>
      dplyr::mutate(
        estimate = dplyr::coalesce(estimate.emm, estimate),
        df       = dplyr::coalesce(df.emm,       as.numeric(df)),
        p.value  = dplyr::coalesce(p.value.emm,  p.value),
        F        = dplyr::case_when(
                     !is.na(p.value.emm) ~ NA_real_,
                     TRUE                ~ as.numeric(as.character(.data$F))
                   ),
        t.ratio  = dplyr::case_when(
                     !is.na(p.value.emm) ~ t.ratio,
                     TRUE                ~ NA_real_
                   )
      ) |>
      dplyr::select(-dplyr::ends_with(".emm"))
  }

  list(anova = anova_tbl, emm_int = emm_int)
}


#' Mixed-effects ANOVA via afex::mixed
#'
#' @inheritParams one_mixed
#' @return List with anova table and emm_int table.
#' @export
one_mixed <- function(
  df_sub, key, formula, form_obj, observed,
  variable, alt_lfc_col, report_lfc, emm,
  alpha, outcome, treatment, control
) {
  # 1) Fit model & extract base ANOVA table (mixed)
  mdl <- suppressMessages(
    afex::mixed(
      formula  = form_obj,
      data     = df_sub,
      progress = FALSE,
      verbose  = FALSE
    )
  )
  fit <- mdl$full_model
  anova_src <- afex::nice(
    mdl, es = NULL, observed = observed,
    sig_symbols = rep("", 4), round_ps = identity
  ) |> tibble::as_tibble()

  # Add estimate for all df = 1 effects
  # 1. Get all fixed effect coefficients and clean up their names
  fixed_effects <- lme4::fixef(fit)
  coef_names_raw <- names(fixed_effects)
  # Clean: remove trailing "1" or "1:" in factor names
  coef_names_clean <- gsub("1(?=:|$)", "", coef_names_raw, perl = TRUE)

  coefficient_table <- tibble::tibble(
    Effect   = coef_names_clean,
    estimate = as.numeric(fixed_effects)
  )
  # Join to afex::nice() anova table output
  anova_src <- anova_src |>
    dplyr::left_join(coefficient_table, by = "Effect") |>
    dplyr::mutate(
      estimate = dplyr::coalesce(estimate, NA_real_)
    )

  # 2) Add partial eta²
  anova_src <- add_pes(anova_src, model = fit)

  # 3) Diagnostics & ICs
  res_vec  <- stats::residuals(fit)
  n_obs    <- length(res_vec)
  p_parm   <- length(stats::coef(fit))
  aic_val  <- stats::AIC(fit)
  bic_val  <- stats::BIC(fit)
  aicc_val <- aic_val + (2 * p_parm * (p_parm + 1)) / (n_obs - p_parm - 1)

  anova_tbl <- anova_src |>
    dplyr::mutate(
      resids = list(setNames(res_vec, df_sub$SampleID)),
      AIC    = aic_val,
      BIC    = bic_val,
      AICc   = aicc_val
    )

  # 4) Normality & homogeneity tests
  tmp <- data.frame(resid = res_vec, grp = df_sub[[variable]]) |>
    dplyr::filter(grp %in% c(treatment, control)) |>
    droplevels()

  sw_p <- tryCatch(stats::shapiro.test(tmp$resid)$p.value, error = \(e) NA_real_)
  ad_p <- tryCatch(nortest::ad.test(tmp$resid)$p.value, error = \(e) NA_real_)
  lt_p <- tryCatch(
    car::leveneTest(resid ~ grp, data = tmp, center = median)[1, "Pr(>F)"] |> as.numeric(),
    error = \(e) NA_real_
  )
  ft_p <- tryCatch(
    stats::var.test(resid ~ grp, data = tmp)$p.value |> as.numeric(),
    error = \(e) NA_real_
  )

  anova_tbl <- anova_tbl |>
    dplyr::mutate(
      SW_pval = sw_p,
      AD_pval = ad_p,
      LT_pval = lt_p,
      FT_pval = ft_p
    ) |>
    dplyr::bind_cols(key)

  # 5) Optional raw log2 fold-change
  if (isTRUE(report_lfc)) {
    lfc_col <- if (!is.null(alt_lfc_col) && alt_lfc_col %in% names(df_sub))
                 alt_lfc_col else variable
    log2fc <- if (
      lfc_col %in% names(df_sub) &&
      is.factor(df_sub[[lfc_col]]) &&
      nlevels(df_sub[[lfc_col]]) == 2
    ) {
      lv <- levels(df_sub[[lfc_col]])
      mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[1]], na.rm = TRUE) -
        mean(df_sub[[outcome]][df_sub[[lfc_col]] == lv[2]], na.rm = TRUE)
    } else {
      NA_real_
    }
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = log2fc)
  } else {
    anova_tbl <- dplyr::mutate(anova_tbl, log2FC = NA_real_)
  }

  # 6) EMM post-hoc handling
  emm_int  <- tibble::tibble()
  repl_tbl <- tibble::tibble()
  if (is.list(emm)) {
    emm_res   <- posthoc_emm_handler(anova_tbl, fit, emm, variable, alpha, key)
    anova_tbl <- emm_res$anova_tbl
    emm_int   <- emm_res$emm_int
    repl_tbl  <- emm_res$repl_tbl
  }

  # 7) splice-in the emmeans main-effect row
  if (nrow(repl_tbl)) {
    anova_tbl <- anova_tbl |>
      dplyr::left_join(repl_tbl, by = "Effect", suffix = c("", ".emm")) |>
      dplyr::mutate(
        estimate = dplyr::coalesce(estimate.emm, estimate),
        df       = dplyr::coalesce(df.emm, as.numeric(df)),
        p.value  = dplyr::coalesce(p.value.emm, p.value),
        F        = dplyr::case_when(
                     !is.na(p.value.emm) ~ NA_real_,
                     TRUE                ~ as.numeric(as.character(.data$F))
        ),
        t.ratio  = dplyr::case_when(
                     !is.na(p.value.emm) ~ t.ratio,
                     TRUE                ~ NA_real_
                   )
      ) |>
      dplyr::select(-dplyr::ends_with(".emm"))
  }

  list(anova = anova_tbl, emm_int = emm_int)
}


#' Wrapper for afex ANOVA/LMM
#'
#' Parallelized wrapper for \code{one_car()} or \code{one_mixed()} across all proteins.
#'
#' @inheritParams robust_test
#' @return List with combined anova and emm_int.
#' @export
afex_test <- function(
  long_df, formula, treatment, control,
  variable      = NULL,
  observed      = NULL,
  alpha         = 0.10,
  report_lfc    = FALSE,
  emm           = NULL,
  alt_lfc_col   = NULL,
  lmer          = FALSE,
  rm_incomplete = TRUE
) {
  if (missing(long_df) || missing(formula)) {
    stop("Both 'long_df' and 'formula' must be supplied.")
  }

  ## 0) set up
  form_obj <- stats::as.formula(formula)
  outcome  <- all.vars(form_obj)[1]
  vars_nom <- all.vars(form_obj)
  extras   <- c("SampleID", "Assay", "OlinkID", "UniProt", "Panel")
  data_sel <- dplyr::select(
    long_df,
    dplyr::all_of(c(vars_nom, extras, alt_lfc_col))
  )

  ## 1) optional NA pruning
  if (isTRUE(rm_incomplete)) {
    preds    <- setdiff(vars_nom, outcome)
    cmpl_idx <- stats::complete.cases(data_sel[preds])
    if (any(!cmpl_idx)) {
      warning("Removed ", sum(!cmpl_idx), " incomplete rows.")
      data_sel <- data_sel[cmpl_idx, ]
      data_sel <- data_sel %>% droplevels()
    }
  }

  ## 2) per‐protein analysis in parallel batches of 100
  grp  <- data_sel %>%
    dplyr::group_by(Assay, OlinkID, UniProt, Panel)
  keys <- grp %>% dplyr::group_keys()
  dfs  <- grp %>% dplyr::group_split()

  params <- purrr::map2(
    dfs,
    split(keys, seq(nrow(keys))),
    \(df, key_row) list(df_sub = df, key = key_row)
  )

  # Use 1 core to process 400 proteins
  chunks <- split(params, ceiling(seq_along(params) / 400))
  # parallel::detectCores() will grab every core on your system
  # so we rescrict that number here.
  # Hard code to 1 for test data so there is only one level of
  # paralellization (the contrast level, not the protein level
  # within each contrast)
  ncores <- 1

  results_list <- parallel::mclapply(
    X        = chunks,
    FUN      = function(chunk) {
      lapply(chunk, function(p) {
        if (isTRUE(lmer)) {
          one_mixed(
            df_sub      = p$df_sub,
            key         = p$key,
            formula     = formula,
            form_obj    = form_obj,
            observed    = observed,
            variable    = variable,
            alt_lfc_col = alt_lfc_col,
            report_lfc  = report_lfc,
            emm         = emm,
            alpha       = alpha,
            outcome     = outcome,
            treatment   = treatment,
            control     = control
          )
        } else {
          one_car(
            df_sub      = p$df_sub,
            key         = p$key,
            formula     = formula,
            form_obj    = form_obj,
            observed    = observed,
            variable    = variable,
            alt_lfc_col = alt_lfc_col,
            report_lfc  = report_lfc,
            emm         = emm,
            alpha       = alpha,
            outcome     = outcome,
            treatment   = treatment,
            control     = control
          )
        }
      })
    },
    mc.cores = ncores
  )
  # If errors occured, they will be burried in results_list
  results_list <- unlist(results_list, recursive = FALSE)

  ## 3) collect & post‐process
  anova_raw   <- purrr::map(results_list, "anova")   %>% dplyr::bind_rows()
  emm_int_tbl <- purrr::map(results_list, "emm_int") %>% purrr::compact() %>% dplyr::bind_rows()

  anova_final <- anova_raw %>%
    dplyr::group_by(Effect) %>%
    dplyr::mutate(
      Adjusted_pval = stats::p.adjust(p.value, method = "BH"),
      SW_padj       = stats::p.adjust(SW_pval, method = "BH"),
      AD_padj       = stats::p.adjust(AD_pval, method = "BH"),
      LT_padj       = stats::p.adjust(LT_pval, method = "BH"),
      FT_padj       = stats::p.adjust(FT_pval, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Threshold = ifelse(Adjusted_pval < alpha, "Significant", "Non-significant"),
      formula   = formula
    ) %>%
    add_raw_npx(long_df) %>%
    sort_results()

  list(anova   = anova_final,
       emm_int = emm_int_tbl)
}


#' Compute interaction contrasts via emmeans::emtrends
#'
#' @param mdl Fitted afex model.
#' @param interaction_spec Character like "Factor:Covariate".
#' @return List of \code{trends}, \code{contrasts}.
#' @export
interaction_emm <- function(mdl, interaction_spec) {
  parts     <- strsplit(interaction_spec, ":")[[1]]
  factor    <- parts[1]
  covariate <- parts[2]

  grid      <- emmeans::emtrends(mdl, specs = factor, var = covariate)
  trends    <- emmeans:::summary.emmGrid(grid, infer = TRUE)
  ctr       <- emmeans::contrast(grid, method = "pairwise")
  contrasts <- emmeans:::summary.emmGrid(ctr, infer = TRUE)
  list(trends = trends, contrasts = contrasts)
}

#' Add partial η² (pes) to ANOVA tables for robust and standard models
#'
#' Compute partial eta‐squared (pes) effect sizes for ANOVA tables derived from
#' robust linear models (\pkg{robustbase}::\code{lmrob}), robust mixed‐effects
#' models (\pkg{robustlmm}::\code{rlmer}), or standard \code{lm}/\code{lmer}
#' models.
#'
#' @param anova_tbl A \code{data.frame} or \code{tibble} containing an ANOVA
#'   table. Must include a \code{df} column and either:
#'   \itemize{
#'     \item \code{F} or \code{Chisq} for \code{lmrob} or \code{rlmerMod}
#'       models, respectively.
#'     \item \code{Sum Sq} (and \code{Residuals}) for standard \code{lm}/\code{lmer}
#'       models.
#'   }
#' @param model A fitted model object. Behavior depends on class:
#'   \itemize{
#'     \item If \code{inherits(model, "lmrob")}, computes
#'       \eqn{pes = (statistic * df) / (statistic * df + df.residual(model))},
#'       where \code{statistic} is \code{F} or \code{Chisq}.
#'     \item If \code{inherits(model, "rlmerMod")}, computes
#'       \eqn{pes = Chisq / (Chisq + df)}.
#'     \item Otherwise, uses \code{stats::anova(model)} to extract sum of
#'       squares and computes
#'       \eqn{pes = Sum Sq / (Sum Sq + Residuals)}.
#'   }
#'
#' @return The input \code{anova_tbl} augmented with a new numeric column
#'   \code{pes} giving the partial η² for each effect.
#' @export
add_pes_rob <- function(anova_tbl, model) {
  # 1) lmrob branch
  if (inherits(model, "lmrob")) {

    ## Which column holds the test statistic?
    stat_col <- dplyr::case_when(
      "F"     %in% names(anova_tbl) ~ "F",
      "Chisq" %in% names(anova_tbl) ~ "Chisq",
      TRUE                          ~ NA_character_
    )
    if (is.na(stat_col))
      stop("Cannot locate test statistic in the supplied ANOVA table.")

    out <- anova_tbl %>%
      dplyr::mutate(
        pes = .data[[stat_col]] * df /
              (.data[[stat_col]] * df + df.residual(model))
      )

    return(out)
  }

  # 2) rlmerMod branch: Wald‐χ²: pes = χ²/(χ²+df1)
  if (inherits(model, "rlmerMod")) {

    ## sanity check: required columns
    stopifnot(all(c("df", "Chisq") %in% names(anova_tbl)))

    out <- anova_tbl %>%
      dplyr::mutate(
        pes = Chisq / (Chisq + df)   # df == df1 in the formula
      )

    return(out)
  }

  # 3) vanilla lm / lme4 branch: Sum Sq pes
  aov_tab <- stats::anova(model)
  aov_tab$Effect <- rownames(aov_tab)

  anova_tbl |>
    dplyr::left_join(
      aov_tab |>
        dplyr::transmute(
          Effect,
          pes = `Sum Sq` / (`Sum Sq` + `Residuals`[["Sum Sq"]][1])
        ),
      by = "Effect"
    )
}


#' Add partial η² (pes) to standard ANOVA tables
#'
#' Uses effectsize::eta_squared() to compute partial η² and merge into results.
#'
#' @param results_tbl Results table with Effect column.
#' @param model Fitted model.
#' @return Results table with pes column.
#' @export
add_pes <- function(results_tbl, model) {
  if ("pes" %in% colnames(results_tbl)) return(results_tbl)

  # Compute partial η² per term (Type-3, full-model)
  pes_tbl <- effectsize::eta_squared(
    model,
    partial     = TRUE,  # partial eta-squeared since generaized isnt supported for lme
    ci          = NULL   # Dont report confidence interval
  ) %>% tibble::as_tibble()

  # Prepare for join: Parameter -> Effect, Eta2_g -> pes
  pes_for_join <- pes_tbl %>%
    dplyr::select(Parameter, Eta2_partial) %>%
    dplyr::rename(
      Effect = Parameter,
      pes    = Eta2_partial
    )

  # Left-join onto results_tbl
  out <- dplyr::left_join(results_tbl, pes_for_join, by = "Effect")

  # Warn if any Effects had no matching pes
  missing_eff <- setdiff(unique(results_tbl$Effect), pes_for_join$Effect)
  if (length(missing_eff) > 0) {
    warning(
      "add_pes(): no pes computed for these Effects: ",
      paste(missing_eff, collapse = ", ")
    )
  }

  out
}

#' Two-sample or paired t-test (Welch)
#'
#' @param long_df Long-format NPX data.
#' @param variable Grouping factor (2 levels).
#' @param pair_id Optional pairing column.
#' @param alpha FDR threshold.
#' @return Tibble of t-test results.
#' @export
t_test <- function(
  long_df,
  variable,
  pair_id = NULL,
  alpha   = 0.05,
  ...
) {
  # Basic checks
  if (missing(long_df) || missing(variable)) {
    stop("Both 'long_df' and 'variable' must be specified.")
  }

  var_levels <- levels(long_df[[variable]])
  if (length(var_levels) != 2) {
    stop("Grouping variable must have exactly 2 levels for Wilcoxon test.")
  }

  # Pre-flight validation
  if (!is.null(pair_id)) {
    # Paired: pivot so each (Assay,Panel,PatientID) is one row
    wide_df <- long_df %>%
      tidyr::pivot_wider(
        id_cols     = c(Assay, OlinkID, UniProt, Panel, PatientID),
        names_from  = !!rlang::ensym(variable),
        values_from = "NPX"
      )

    # count complete pairs by assay
    too_few <- wide_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::summarise(
        n_pairs = sum(
          !is.na(.data[[ var_levels[1] ]]) &
          !is.na(.data[[ var_levels[2] ]]),
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_pairs < 2)

    if (nrow(too_few) > 0) {
      msgs <- sprintf(
        "%s/%s: only %d complete pairs",
        too_few$Assay, too_few$OlinkID, too_few$n_pairs
      )
      stop(
        "Insufficient paired observations for assays:\n",
        paste0("  - ", msgs, collapse = "\n")
      )
    }

    # Run the *paired* t-test
    message(sprintf("Paired t-test on %s vs %s.", var_levels[1], var_levels[2]))
    results <- wide_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::do(broom::tidy(
        t.test(
          x      = .[[ var_levels[1] ]],
          y      = .[[ var_levels[2] ]],
          paired = TRUE,
          ...
        )
      )) %>%
      dplyr::ungroup()

  } else {
    # Unpaired: count directly in the long data
    counts <- long_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel, !!rlang::ensym(variable)) %>%
      dplyr::summarise(
        n_obs = sum(!is.na(NPX)),
        .groups = "drop"
      )

    counts_wide <- counts %>%
      tidyr::pivot_wider(
        names_from  = !!rlang::ensym(variable),
        values_from = "n_obs",
        values_fill = 0
      )

    too_few <- counts_wide %>%
      dplyr::filter(
        .data[[ var_levels[1] ]] < 2 |
        .data[[ var_levels[2] ]] < 2
      )

    if (nrow(too_few) > 0) {
      msgs <- sprintf(
        "%s/%s: %s(%d), %s(%d)",
        too_few$Assay, too_few$OlinkID,
        var_levels[1], too_few[[ var_levels[1] ]],
        var_levels[2], too_few[[ var_levels[2] ]]
      )
      stop(
        "Insufficient observations for assays:\n",
        paste0("  - ", msgs, collapse = "\n")
      )
    }

    # Run the *unpaired* t-test
    message(sprintf("t-test on %s vs %s.", var_levels[1], var_levels[2]))
    results <- long_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::do(broom::tidy(
        t.test(
          NPX ~ !!rlang::ensym(variable),
          data = .,
          ...
        )
      )) %>%
      dplyr::ungroup() %>%
      dplyr::select(-dplyr::any_of(c("estimate1","estimate2")))
  }

  # Annotate & return
  results <- results %>%
    dplyr::mutate(
      Adjusted_pval   = p.adjust(p.value, method = "fdr"),
      Threshold       = ifelse(Adjusted_pval < alpha, "Significant", "Non-significant"),
      formula         = paste(var_levels[1], "-", var_levels[2]),
      estimate_method = "mean difference",
      log2FC          = estimate
    ) %>%
    sort_results()
  }


#' Wilcoxon rank-sum or signed-rank test
#'
#' @param long_df Long-format NPX data.
#' @param variable Grouping factor (2 levels).
#' @param pair_id Optional pairing column.
#' @param alpha FDR threshold.
#' @return Tibble of Wilcoxon results.
#' @export
wilcox_test <- function(
  long_df,
  variable,
  pair_id = NULL,
  alpha   = 0.05,
  ...
) {
  # basic input checks
  if (missing(long_df) || missing(variable)) {
    stop("Both 'long_df' and 'variable' must be specified.")
  }
  var_levels <- levels(long_df[[variable]])

  if (length(var_levels) != 2) {
    stop("Grouping variable must have exactly 2 levels for Wilcoxon test.")
  }

  # decide paired vs unpaired
  if (!is.null(pair_id)) {
    # pivot so each (Assay,Panel,PatientID) is one row
    wide_df <- long_df %>%
      tidyr::pivot_wider(
        id_cols     = c(Assay, OlinkID, UniProt, Panel, PatientID),
        names_from  = !!rlang::ensym(variable),
        values_from = "NPX"
      )

    # count complete pairs
    too_few <- wide_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::summarise(
        n_pairs = sum(
          !is.na(.data[[ var_levels[1] ]]) &
          !is.na(.data[[ var_levels[2] ]]),
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_pairs < 2)

    if (nrow(too_few) > 0) {
      msgs <- sprintf(
        "%s/%s: only %d complete pairs",
        too_few$Assay, too_few$OlinkID, too_few$n_pairs
      )
      stop(
        "Insufficient paired observations for assays:\n",
        paste0("  - ", msgs, collapse = "\n")
      )
    }

    message(sprintf("Paired Wilcoxon test on %s vs %s.", var_levels[1], var_levels[2]))
    results <- wide_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::do(broom::tidy(
        stats::wilcox.test(
          x        = .[[ var_levels[1] ]],
          y        = .[[ var_levels[2] ]],
          paired   = TRUE,
          conf.int = TRUE,
          ...
        )
      )) %>%
      dplyr::ungroup()

  } else {
    # unpaired: require ≥2 obs per group
    counts <- long_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel, !!rlang::ensym(variable)) %>%
      dplyr::summarise(n_obs = sum(!is.na(NPX)), .groups = "drop")

    counts_wide <- counts %>%
      tidyr::pivot_wider(
        names_from  = !!rlang::ensym(variable),
        values_from = "n_obs",
        values_fill = 0
      )

    too_few <- counts_wide %>%
      dplyr::filter(
        .data[[ var_levels[1] ]] < 2 |
        .data[[ var_levels[2] ]] < 2
      )

    if (nrow(too_few) > 0) {
      msgs <- sprintf(
        "%s/%s: %s(%d), %s(%d)",
        too_few$Assay, too_few$OlinkID,
        var_levels[1], too_few[[ var_levels[1] ]],
        var_levels[2], too_few[[ var_levels[2] ]]
      )
      stop(
        "Insufficient observations for assays:\n",
        paste0("  - ", msgs, collapse = "\n")
      )
    }

    message(sprintf("Wilcoxon test on %s vs %s.", var_levels[1], var_levels[2]))
    results <- long_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::do(broom::tidy(
        stats::wilcox.test(
          NPX      ~ !!rlang::ensym(variable),
          data     = .,
          conf.int = TRUE,
          ...
        )
      )) %>%
      dplyr::ungroup()
  }

  # annotate, compute true log2FC, and return
  results <- results %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Adjusted_pval    = p.adjust(p.value, method = "fdr"),
      Threshold        = ifelse(Adjusted_pval < alpha, "Significant", "Non-significant"),
      formula          = paste(var_levels[1], "-", var_levels[2]),
      estimate_method  = "Hodges–Lehmann location shift",
      log2FC = {
        subdf <- long_df[long_df$Assay == Assay, ]
        mean(subdf$NPX[subdf[[variable]] == var_levels[1]], na.rm = TRUE) -
        mean(subdf$NPX[subdf[[variable]] == var_levels[2]], na.rm = TRUE)
      }
    ) %>%
    dplyr::ungroup() %>%
    sort_results()
}


#' Kruskal-Wallis or Friedman test
#'
#' @param long_df Long-format NPX data.
#' @param variable Grouping factor.
#' @param alpha FDR threshold.
#' @param dependence If TRUE run Friedman test with subject ID.
#' @param subject Subject column for Friedman.
#' @return Tibble of results with epsilon² or Kendall's W.
#' @export
kruskal_test <- function(
  long_df,
  variable,
  alpha      = 0.05,
  dependence = FALSE,
  subject    = NULL,
  verbose    = TRUE
) {
  # Basic checks
  if (missing(long_df) || missing(variable)) {
    stop("Both 'long_df' and 'variable' must be specified.")
  }

  var_levels <- levels(long_df[[variable]])

  k <- length(var_levels)

  if (dependence) {
    # Friedman: need subject blocks
    if (is.null(subject))
      stop("'subject' must be provided for Friedman test.")
    long_df[[subject]] <- factor(long_df[[subject]])

    if (verbose) {
      message("Performing Friedman test on ", variable)
    }
    # Pivot so each (Assay,Panel,subject) is one row with one column per level
    wide_df <- long_df %>%
      tidyr::pivot_wider(
        id_cols     = c(Assay, OlinkID, UniProt, Panel, !!rlang::ensym(subject)),
        names_from  = !!rlang::ensym(variable),
        values_from = "NPX"
      )

    # Pre-flight: count complete blocks
    too_few <- wide_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::summarise(
        n_blocks = sum(
          !is.na(.data[[var_levels[1]]]) &
          !is.na(.data[[var_levels[2]]]),
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_blocks < 2)

    if (nrow(too_few) > 0) {
      msgs <- sprintf(
        "%s/%s: only %d complete blocks",
        too_few$Assay, too_few$OlinkID, too_few$n_blocks
      )
      stop(
        "Insufficient complete blocks for Friedman test:\n",
        paste0("  - ", msgs, collapse = "\n")
      )
    }

    # Run Friedman per assay
    res <- wide_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      do({
        d <- .
        ft <- friedman.test(as.matrix(d[, var_levels]))
        tibble(
          term      = variable,
          statistic = unname(ft$statistic),
          df        = unname(ft$parameter),
          p.value   = ft$p.value,
          N         = nrow(d)
        )
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        estimate        = statistic / (N * (k - 1)),
        estimate_method = "kendall_W",
        formula         = paste0("NPX ~ ", variable, " | ", subject)
      ) %>%
      dplyr::select(-N)

  } else {
    # Kruskal-Wallis
    if (verbose) {
      message("Performing Kruskal-Wallis test on ", variable)
    }
    # Run the test per assay
    res <- long_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::do(broom::tidy(
        stats::kruskal.test(NPX ~ .data[[variable]], data = .)
      )) %>%
      dplyr::ungroup() %>%
      dplyr::rename(df = parameter)

    # Compute per‐assay n_obs
    counts <- long_df %>%
      dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
      dplyr::summarise(
        n_assay = sum(!is.na(NPX)),
        .groups = "drop"
      )

    # Join in and compute epsilon² as estimate
    res <- res %>%
      dplyr::left_join(counts, by = c("Assay","OlinkID","UniProt","Panel")) %>%
      dplyr::mutate(
        estimate        = (statistic - k + 1) / (n_assay - k),
        estimate_method = "epsilon_squared",
        formula         = paste0("NPX ~ ", variable)
      ) %>%
      dplyr::select(-n_assay)
  }

  # Adjust p-values and flag significance
  res <- res %>%
    dplyr::mutate(
      Adjusted_pval = p.adjust(p.value, method = "fdr"),
      Threshold     = ifelse(Adjusted_pval < alpha, "Significant", "Non-significant")
    )

  # Add log2FC for two-level factors
  if (k == 2) {
    res <- res %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        log2FC = {
          subdf <- long_df[long_df$Assay == Assay, ]
          lv1   <- var_levels[1]
          lv2   <- var_levels[2]
          mean(subdf$NPX[subdf[[variable]] == lv1], na.rm = TRUE) -
          mean(subdf$NPX[subdf[[variable]] == lv2], na.rm = TRUE)
        }
      ) %>%
      dplyr::ungroup()
  } else {
    res <- res %>% dplyr::mutate(log2FC = NA_real_)
  }

  res <- res %>%
    sort_results()
}


#' Ordinal regression (ordinal ANOVA)
#'
#' @param long_df Long-format NPX data.
#' @param variable Factor(s).
#' @param covariates Covariates.
#' @param alpha Threshold.
#' @param outcome Ordinal outcome column.
#' @return Tibble with ordinal ANOVA results.
#' @export
ordinal_regression <- function(
  long_df,
  variable,
  covariates       = NULL,
  alpha            = 0.05,
  outcome          = "NPX",
  return.covariates= FALSE,
  verbose          = TRUE
) {
  # Basic checks
  if (missing(long_df) || missing(variable)) {
    stop("Both 'long_df' and 'variable' must be supplied.")
  }

  var_levels <- levels(long_df[[variable]])

  # Rank outcome
  long_df[[outcome]] <- as.numeric(
    factor(long_df[[outcome]], ordered = TRUE)
  )

  # Build formula
  if (is.null(covariates) || length(covariates)==0) {
    formula_string <- paste0(outcome, " ~ ", paste(variable, collapse="*"))
  } else {
    formula_string <- paste0(
      outcome, " ~ ", paste(variable, collapse="*"), " + ", paste(covariates, collapse="+")
    )
  }
  if (verbose) message("Ordinal ANOVA model: ", formula_string)

  # PER-ASSAY FIT & TYPE III ANOVA
  results <- long_df %>%
    dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
    dplyr::do({
      mod <- stats::lm(as.formula(formula_string), data = .)
      tbl <- car::Anova(mod, type = 3)
      data.frame(
        term      = rownames(tbl),
        df        = tbl[,"Df"],
        statistic = tbl[,"F value"],
        p.value   = tbl[,"Pr(>F)"],
        resid_df  = mod$df.residual,
        stringsAsFactors = FALSE
      )
    }) %>%
    dplyr::ungroup() %>%
    # drop intercept/residual
    dplyr::filter(!term %in% c("(Intercept)","Residuals"))

  # FDR + Partial ETA² + formating
  results <- results %>%
    dplyr::mutate(
      Adjusted_pval   = p.adjust(p.value, method="fdr"),
      Threshold       = ifelse(Adjusted_pval < alpha, "Significant","Non-significant"),
      estimate        = (statistic * df)/(statistic * df + resid_df),
      estimate_method = "partial η²",
      formula         = formula_string
    ) %>%
    dplyr::select(-resid_df)

  # Optionally drop covariates
  if (!return.covariates && !is.null(covariates)) {
    results <- results %>% dplyr::filter(!term %in% covariates)
  }

  # add log2FC for two-level factors
  # compute levels once
  lvls_list <- lapply(variable, function(v) levels(long_df[[v]]))
  names(lvls_list) <- variable

  results <- results %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      log2FC = {
        if (term %in% variable && length(lvls_list[[term]])==2) {
          subdf <- long_df[long_df$Assay==Assay, ]
          lvl1 <- lvls_list[[term]][1]
          lvl2 <- lvls_list[[term]][2]
          mean(subdf[[outcome]][subdf[[term]]==lvl1], na.rm=TRUE) -
          mean(subdf[[outcome]][subdf[[term]]==lvl2], na.rm=TRUE)
        } else {
          NA_real_
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    sort_results()
}


#' Linear mixed model via lmerTest
#'
#' @param long_df Long-format NPX data.
#' @param variable Fixed factors.
#' @param random Random effect column.
#' @param covariates Optional covariates.
#' @param outcome Outcome column.
#' @param alpha Threshold.
#' @return Tibble of LMM results.
#' @export
linear_mixed_model <- function(
  long_df,
  variable,
  random,
  covariates       = NULL,
  outcome          = "NPX",
  alpha            = 0.05,
  return.covariates= FALSE,
  verbose          = TRUE,
  ...
) {
  # Input checks
  if (missing(long_df) || missing(variable) || missing(random)) {
    stop("`long_df`, `variable`, and `random` must all be specified.")
  }

  var_levels <- levels(long_df[[variable]])

  # Build formula
  if (is.null(covariates) || length(covariates)==0) {
    formula_string <- paste0(
      outcome, " ~ ", paste(variable, collapse="*"),
      " + (1|", random, ")"
    )
  } else {
    formula_string <- paste0(
      outcome, " ~ ", paste(variable, collapse="*"),
      " + ", paste(covariates, collapse = "+"),
      " + (1|", random, ")"
    )
  }
  if (verbose) message("LMM formula: ", formula_string)

  # Per Assay model fit & type III ANOVA
  raw <- long_df %>%
    dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
    dplyr::do({
      # Check if this test is basically a paired t-test
      if (lme4::isSingular(lmm, tol = 1e-4)) {
        warning(
          sprintf(
            "Assay %s: singular LMM fit (random intercept variance ≈ 0).
            This is equivalent to a paired t-test. Perhaps you only have one observation per %s?",
            unique(.$Assay),
            random
          )
        )
      }

      aov <- stats::anova(lmm, type = 3)
      data.frame(
        term      = rownames(aov),
        df        = aov$NumDF,
        denDF     = aov$DenDF,
        statistic = aov$`F value`,
        p.value   = aov$`Pr(>F)`,
        stringsAsFactors = FALSE
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term != "(Intercept)")

  # Add FDR, Threshold and Effect size
  results <- raw %>%
    dplyr::mutate(
      Adjusted_pval   = p.adjust(p.value, method = "fdr"),
      Threshold       = ifelse(Adjusted_pval < alpha, "Significant", "Non-significant"),
      estimate        = (statistic * df) / (statistic * df + denDF),
      estimate_method = "partial η²",
      formula         = formula_string
    )

  # Optionally drop covariates
  if (!return.covariates && !is.null(covariates)) {
    results <- results %>% dplyr::filter(!term %in% covariates)
  }

  # Prepare levels for log2FC calculation
  # only compute for fixed‐effect terms
  lvl_map <- lapply(variable, function(v) levels(long_df[[v]]))
  names(lvl_map) <- variable

  # Compute log2FC for two-level terms
  results <- results %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      log2FC = {
        if (term %in% variable && length(lvl_map[[term]])==2) {
          subdf <- long_df[long_df$Assay==Assay, ]
          lvl1  <- lvl_map[[term]][1]
          lvl2  <- lvl_map[[term]][2]
          mean(subdf[[outcome]][subdf[[term]]==lvl1], na.rm=TRUE) -
          mean(subdf[[outcome]][subdf[[term]]==lvl2], na.rm=TRUE)
        } else {
          NA_real_
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    sort_results()
}
