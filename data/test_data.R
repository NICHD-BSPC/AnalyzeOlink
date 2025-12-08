# Make synthetic, predictable PEA & ELISA test data and matching metadata
# Data contains:
# - Two "Group" levels: "disease" and "control"
# - All downstream tests focus on disease_vs_control
# - Metadata columns are: Group, Cov1, Cov2
#     * Group: "disease" / "control"
#     * Cov1: continuous covariate (analogous to Age)
#     * Cov2: binary covariate (analogous to Sex)
# - We add QC_Warning and Assay_Warning columns to PEA as mostly "PASS" with a few deterministic "WARN"
# - We add deterministic outliers:
#     PEA:   SampleID "A1" NPX set to 7 for all assays
#     ELISA: SampleID "A1" ELISA set to 12.00 for all assays
# - PEA/ELISA file pairs share the same structure

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(OlinkAnalyze)
})

set.seed(123)

# ------------------------------------------------------------------------------
# Load base PEA test data from the OlinkAnalyze package to act as our base
# ------------------------------------------------------------------------------

pea_test_data_raw <- get("npx_data1")

# ------------------------------------------------------------
# SAFELY REDUCE DATASET TO 10 assays per panel (we have two panels)
# ------------------------------------------------------------

set.seed(123)

# Number of assays to keep per panel
assays_per_panel <- 10L

# List assays by panel
panel_assays <- pea_test_data_raw %>%
  distinct(Panel, Assay) %>%
  group_by(Panel) %>%
  summarize(assays = list(unique(Assay)), .groups = "drop")

# Sample 50 from each panel
panel_assays$assays_keep <- lapply(panel_assays$assays, function(x) {
  sample(x, size = assays_per_panel)
})

# Combine into a single vector
assays_keep <- unlist(panel_assays$assays_keep, use.names = FALSE)

# Filter raw PEA data to just these assays
pea_test_data_raw <- pea_test_data_raw %>%
  filter(Assay %in% assays_keep)

# ------------------------------------------------------------------------------
# Build metadata
# ------------------------------------------------------------------------------

pea_test_metadata <- pea_test_data_raw %>%
  distinct(SampleID) %>%
  arrange(SampleID) %>%
  mutate(
    row   = dplyr::row_number(),
    Group = if_else(row <= floor(n()/2), "disease", "control"),
    Cov1  = runif(n(), min = 0.9, max = 68.0),
    Cov2  = sample(c("M","F"), size = n(), replace = TRUE)
  ) %>%
  select(SampleID, Group, Cov1, Cov2)

# ------------------------------------------------------------------------------
# Per-assay effects by Group: continuous effect sizes for disease_vs_control
# ------------------------------------------------------------------------------

all_assays <- sort(unique(pea_test_data_raw$Assay))
n_assay    <- length(all_assays)

# Small deterministic per-assay offset to prevent ties
a_off <- if (n_assay > 1) {
  (seq_along(all_assays) - 1) / (n_assay - 1) * 0.4 - 0.2
} else {
  0
}
assay_offsets <- tibble::tibble(Assay = all_assays, a_off = a_off)

# Assign assay-level patterns to create a realistic mix of null/up/down
pattern_levels <- c("null", "disease_up", "disease_down")
pattern_probs  <- c(0.68, 0.16, 0.16)

assay_patterns <- tibble::tibble(
  Assay   = all_assays,
  pattern = sample(pattern_levels, n_assay, replace = TRUE, prob = pattern_probs)
)

# Continuous baseline means per assay
mu_control <- rnorm(n_assay, mean = 0, sd = 0.30)
mu_disease <- mu_control + rnorm(n_assay, mean = 0, sd = 0.10)

idx_d_up <- assay_patterns$pattern == "disease_up"
idx_d_dn <- assay_patterns$pattern == "disease_down"

# disease-specific patterns with continuous magnitudes
if (any(idx_d_up)) {
  mu_disease[idx_d_up] <- mu_control[idx_d_up] +
    abs(rnorm(sum(idx_d_up), mean = 0.9, sd = 0.45))
}
if (any(idx_d_dn)) {
  mu_disease[idx_d_dn] <- mu_control[idx_d_dn] -
    abs(rnorm(sum(idx_d_dn), mean = 0.9, sd = 0.45))
}

assay_patterns <- assay_patterns %>%
  mutate(
    mu_control = mu_control,
    mu_disease = mu_disease
  )

# ------------------------------------------------------------------------------
# Build PEA test data with the continuous effects
# ------------------------------------------------------------------------------

pea_test_data <- pea_test_data_raw %>%
  select(-Treatment) %>% # safety if base dataset has Treatment
  left_join(pea_test_metadata %>% select(SampleID, Group), by = "SampleID") %>%
  left_join(assay_offsets, by = "Assay") %>%
  left_join(
    assay_patterns %>% select(Assay, mu_control, mu_disease),
    by = "Assay"
  ) %>%
  mutate(
    base = case_when(
      Group == "control" ~ mu_control,
      Group == "disease" ~ mu_disease,
      TRUE               ~ 0
    ),
    # Noise chosen to yield classic volcano spread
    NPX = base + a_off + rnorm(n(), sd = 0.5)
  ) %>%
  select(-base, -mu_control, -mu_disease)

# Deterministic PEA outlier for A1
pea_test_data <- pea_test_data %>%
  mutate(NPX = if_else(SampleID == "A1", 7, NPX))

# ------------------------------------------------------------------------------
# Add QC_Warning and Assay_Warning to PEA with mostly PASS, and few deterministic WARN
# ------------------------------------------------------------------------------

unique_samples_sorted <- sort(unique(pea_test_data$SampleID))
unique_assays_sorted  <- sort(unique(pea_test_data$Assay))

qc_warn_ids_pea <- head(unique_samples_sorted, 2)
assay_warns_pea <- head(unique_assays_sorted,  2)

pea_test_data <- pea_test_data %>%
  mutate(
    QC_Warning    = if_else(SampleID %in% qc_warn_ids_pea, "WARN", "PASS"),
    Assay_Warning = if_else(Assay %in% assay_warns_pea, "WARN", "PASS")
  ) %>%
  select(
    SampleID, Index, OlinkID, UniProt, Assay, MissingFreq, Panel, Panel_Version,
    PlateID, QC_Warning, LOD, NPX, Subject, Site, Time, Project, Assay_Warning
  )

# ------------------------------------------------------------------------------
# Build mirrored ELISA test data with no QC columns from Olink
#   (automatically respects the reduced set of assays because it
#    starts from pea_test_data_raw, which we already filtered)
# ------------------------------------------------------------------------------

elisa_test_metadata <- pea_test_metadata

elisa_skeleton <- pea_test_data_raw %>%
  select(Assay, Panel, OlinkID, UniProt, SampleID) %>%
  distinct() %>%
  left_join(elisa_test_metadata %>% select(SampleID, Group), by = "SampleID") %>%
  left_join(assay_offsets, by = "Assay") %>%
  left_join(
    assay_patterns %>% select(Assay, mu_control, mu_disease),
    by = "Assay"
  )

elisa_test_data <- elisa_skeleton %>%
  mutate(
    base_pea = case_when(
      Group == "control" ~ mu_control,
      Group == "disease" ~ mu_disease,
      TRUE               ~ 0
    ),
    # Map PEA-scale effects to ELISA-scale with noise, then clip to [1,6]
    ELISA = 3 + 1.5 * (base_pea + a_off) + rnorm(n(), sd = 0.15),
    ELISA = pmin(pmax(ELISA, 1), 6)
  ) %>%
  select(Assay, Panel, OlinkID, UniProt, SampleID, ELISA)

# Deterministic ELISA outlier for SampleID "A1"
elisa_test_data <- elisa_test_data %>%
  mutate(ELISA = if_else(SampleID == "A1", 12.00, ELISA))

# ------------------------------------------------------------------------------
# Save outputs assuming this script is run from repo root (./)
# ------------------------------------------------------------------------------

dir.create("./data/PEA",   showWarnings = FALSE, recursive = TRUE)
dir.create("./data/ELISA", showWarnings = FALSE, recursive = TRUE)

write_tsv(pea_test_data,       "./data/PEA/pea_test_data.tsv")
write_tsv(pea_test_metadata,   "./data/PEA/pea_test_metadata.tsv")
write_tsv(elisa_test_data,     "./data/ELISA/elisa_test_data.tsv")
write_tsv(elisa_test_metadata, "./data/ELISA/elisa_test_metadata.tsv")

# ------------------------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------------------------
cat("\nPEA test data:")
head(pea_test_data)
cat("\nPEA test metadata:")
head(pea_test_metadata)
cat("\nELISA test data:")
head(elisa_test_data)
cat("\nELISA test metadata:")
head(elisa_test_metadata)

cat("\nPEA means by Group (joined to metadata):\n")
pea_test_data %>%
  left_join(pea_test_metadata %>% select(SampleID, Group), by = "SampleID") %>%
  group_by(Group) %>%
  summarise(mean_NPX = mean(NPX), sd_NPX = sd(NPX), n = dplyr::n(), .groups = "drop") %>%
  print(n = Inf)

cat("\nELISA means by Group (joined to metadata):\n")
elisa_test_data %>%
  left_join(elisa_test_metadata %>% select(SampleID, Group), by = "SampleID") %>%
  group_by(Group) %>%
  summarise(mean_ELISA = mean(ELISA), sd_ELISA = sd(ELISA), n = dplyr::n(), .groups = "drop") %>%
  print(n = Inf)

cat("\nQC/Assay warnings (PEA):\n")
print(table(pea_test_data$QC_Warning))
print(table(pea_test_data$Assay_Warning))
