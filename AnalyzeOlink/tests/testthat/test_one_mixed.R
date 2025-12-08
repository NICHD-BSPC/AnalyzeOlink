testthat::context("Testing afex_test()'s call to one_mixed() with paired pre/post treatment data")

# ----------------------------------------------------------------------
# 1) Create a hardcoded toy paired dataset that matches the
#    structure expected by afex_test()/one_mixed()
# ----------------------------------------------------------------------

base::set.seed(123)
n_patients  <- 11
patient_ids <- base::paste0("dis", base::seq_len(n_patients))   # dis1, dis2, ...

# suffix "a" for pre, "b" for post
sample_suffix <- base::rep(c("a", "b"), times = n_patients)

toy_df <- tibble::tibble(
  PatientID = base::rep(patient_ids, each = 2),
  PrePost   = base::factor(
    base::rep(c("pre", "post"), times = n_patients),
    levels = c("pre", "post")
  ),
  NPX      = stats::runif(2 * n_patients, min = -10, max = 10),
  Age_c    = stats::rnorm(2 * n_patients),
  NSS_17_c = stats::rnorm(2 * n_patients),
  ANO_c    = stats::rnorm(2 * n_patients),
  Sex      = base::factor(base::sample(c("M", "F"), 2 * n_patients, replace = TRUE)),
  Assay    = "SESTD1",
  OlinkID  = "OID20790",
  UniProt  = "Q86VW0",
  Panel    = "Neurology"
)

toy_df <- dplyr::mutate(
  toy_df,
  SampleID = base::paste0(PatientID, sample_suffix)
)

long_df <- toy_df

formula_str <- "NPX ~ PrePost * Age_c + NSS_17_c + Sex + ANO_c + (1 | PatientID)"
observed    <- c("Sex", "Age_c", "NSS_17_c", "ANO_c")
variable    <- "PrePost"
alt_lfc_col <- NULL
report_lfc  <- TRUE
emm         <- NULL
alpha       <- 0.10
outcome     <- "NPX"
treatment   <- "post"
control     <- "pre"

# ----------------------------------------------------------------------
# 2) Run afex_test with lmer = TRUE which calls one_mixed() internally
# ----------------------------------------------------------------------

res_all <- AnalyzeOlink::afex_test(
  long_df       = long_df,
  formula       = formula_str,
  treatment     = treatment,
  control       = control,
  variable      = variable,
  observed      = observed,
  alpha         = alpha,
  report_lfc    = report_lfc,
  emm           = emm,
  alt_lfc_col   = alt_lfc_col,
  lmer          = TRUE,
  rm_incomplete = TRUE
)

testthat::test_that("afex_test(lmer=TRUE) returns list with 'anova' and 'emm_int'", {
  testthat::expect_type(res_all, "list")
  testthat::expect_named(res_all, c("anova", "emm_int"))
})

res <- res_all$anova

testthat::test_that("anova tibble is not empty", {
  testthat::expect_gt(base::nrow(res), 0)
})

testthat::test_that("ANOVA table has expected columns", {
  expected_cols <- c(
    "Effect", "df", "F", "p.value", "pes", "estimate",
    "AIC", "BIC", "AICc", "resids",
    "SW_pval", "AD_pval", "LT_pval", "FT_pval",
    "Assay", "OlinkID", "UniProt", "Panel", "log2FC"
  )
  testthat::expect_true(base::all(expected_cols %in% base::colnames(res)))
})

testthat::test_that("log2FC matches difference of means (pre - post)", {
  pre_mean  <- base::mean(long_df$NPX[long_df$PrePost == "pre"],  na.rm = TRUE)
  post_mean <- base::mean(long_df$NPX[long_df$PrePost == "post"], na.rm = TRUE)
  expected  <- pre_mean - post_mean

  # log2FC is computed from alt_lfc_col = Pre, so it should match pre - post
  log2fc_vals <- base::unique(res$log2FC[!base::is.na(res$log2FC)])
  testthat::expect_gt(base::length(log2fc_vals), 0)
  testthat::expect_true(base::all(base::abs(log2fc_vals - expected) < 1e-5))
})

testthat::test_that("diagnostic p-values are between 0 and 1", {
  for (col in c("SW_pval", "AD_pval", "LT_pval", "FT_pval")) {
    pvs <- base::unique(res[[col]])
    testthat::expect_true(base::all(!base::is.na(pvs)))
    testthat::expect_true(base::all(pvs >= 0 & pvs <= 1))
  }
})
