context("Testing afex_test()/one_car() on PEA test data")

# ----------------------------------------------------------------------
# 1) Read in PEA test data & metadata and join
# ----------------------------------------------------------------------

pea_path  <- testthat::test_path("..", "..", "..", "data", "PEA", "pea_test_data.tsv")
meta_path <- testthat::test_path("..", "..", "..", "data", "PEA", "pea_test_metadata.tsv")

pea_df  <- readr::read_tsv(pea_path,  show_col_types = FALSE)
meta_df <- readr::read_tsv(meta_path, show_col_types = FALSE)

joined <- pea_df %>%
  dplyr::left_join(meta_df, by = "SampleID") %>%
  # keep samples with complete predictors
  dplyr::filter(
    !is.na(Group),
    !is.na(Cov1),
    !is.na(Cov2)
  ) %>%
  # restrict to the two levels used for log2FC / EMM tests
  dplyr::filter(Group %in% c("control", "disease"))

if (nrow(joined) == 0) {
  stop("Joined test data has no rows with Group in {control, disease}.")
}

# ----------------------------------------------------------------------
# 2) Restrict to the first 10 assays with both control and disease
# ----------------------------------------------------------------------

group_info <- joined %>%
  dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>%
  dplyr::summarise(
    n       = dplyr::n(),
    n_treat = dplyr::n_distinct(Group),
    .groups = "drop"
  )

valid_group <- group_info %>%
  dplyr::filter(n_treat == 2) %>%
  head(10)

if (nrow(valid_group) == 0) {
  stop("No assay in pea_test_data.tsv has both control and disease samples.")
}

long_df <- joined %>%
  dplyr::semi_join(
    valid_group,
    by = c("Assay", "OlinkID", "UniProt", "Panel")
  ) %>%
  dplyr::mutate(
    SampleID  = factor(SampleID),
    Group = factor(Group, levels = c("control", "disease")),
    Cov2       = factor(Cov2),
    Cov1_c     = Cov1 - mean(Cov1, na.rm = TRUE)
  )

if (nrow(long_df) < 3) {
  stop("Not enough observations in the selected assays to run afex_test().")
}

# ----------------------------------------------------------------------
# 3) Run afex_test() which calls one_car() internally
# ----------------------------------------------------------------------

afex_res <- afex_test(
  long_df       = long_df,
  formula       = "NPX ~ Group * Cov1_c + Cov2 + Error(SampleID)",
  treatment     = "disease",
  control       = "control",
  variable      = "Group",
  observed      = c("Cov2", "Cov1_c"),
  alpha         = 0.10,
  report_lfc    = TRUE,
  emm           = list(Group = c("disease", "control")),
  alt_lfc_col   = NULL,
  lmer          = FALSE,
  rm_incomplete = TRUE
)

res <- afex_res$anova

# ----------------------------------------------------------------------
# 4) Perform tests
# ----------------------------------------------------------------------

test_that("afex_test() returns list with 'anova' and 'emm_int'", {
  expect_type(afex_res, "list")
  expect_named(afex_res, c("anova", "emm_int"))
})

test_that("anova tibble is not empty", {
  expect_gt(nrow(res), 0)
})

test_that("ANOVA table has expected columns", {
  expected_cols <- c(
    "Effect", "df", "F", "p.value", "pes", "estimate",
    "AIC", "BIC", "AICc", "resids",
    "SW_pval", "AD_pval", "LT_pval", "FT_pval",
    "Assay", "OlinkID", "UniProt", "Panel", "log2FC"
  )
  expect_true(all(expected_cols %in% colnames(res)))
})

test_that("log2FC matches difference of means (control - disease) per assay", {
  # use only rows where log2FC is defined
  res_non_na <- res %>% dplyr::filter(!is.na(log2FC))

  assays <- res_non_na %>%
    dplyr::distinct(Assay, OlinkID, UniProt, Panel)

  expect_gt(nrow(assays), 0)

  for (i in seq_len(nrow(assays))) {
    a <- assays[i, ]

    df_sub <- long_df %>%
      dplyr::semi_join(a, by = c("Assay", "OlinkID", "UniProt", "Panel"))

    control_mean <- mean(df_sub$NPX[df_sub$Group == "control"], na.rm = TRUE)
    disease_mean <- mean(df_sub$NPX[df_sub$Group == "disease"], na.rm = TRUE)
    expected     <- control_mean - disease_mean

    res_sub <- res_non_na %>%
      dplyr::semi_join(a, by = c("Assay", "OlinkID", "UniProt", "Panel"))

    log2fc_vals <- unique(res_sub$log2FC[!is.na(res_sub$log2FC)])

    expect_gt(length(log2fc_vals), 0)
    expect_true(all(abs(log2fc_vals - expected) < 1e-8))
  }
})

test_that("diagnostic p-values are between 0 and 1", {
  for (col in c("SW_pval", "AD_pval", "LT_pval", "FT_pval")) {
    pvs <- unique(res[[col]])
    expect_true(all(!is.na(pvs)))
    expect_true(all(pvs >= 0 & pvs <= 1))
  }
})
