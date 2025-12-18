# Olink® PEA Analysis Workflow

**AnalyzeOlink** provides a configurable report-based framework for analyzing
Olink® Proximity Extension Assay (PEA) data, focusing on differential abundance
testing and downstream validation. Supported statistical methods include Type
III ANCOVA and linear mixed-effects models (via the `afex` package), t-tests
and Wilcoxon tests (paired and unpaired), as well as Kruskal–Wallis and
Friedman tests.

The package includes nine R Markdown analysis reports that cover key workflow
components: assumption checks for ANOVA and LMER models, leave-one-out
sensitivity analyses for model terms and sample-level robustness, ELISA
validation, non-parametric covariate effects, and NPX correlations with
metadata variables.

An annotations module integrates Human Protein Atlas (HPA) metadata, including
tissue- and cell-type–specific expression profiles and Tau specificity scores.
These annotations can be used to color interactive downstream visualizations
such as volcano plots and estimate-vs-estimate plots.

Each phase of the analysis is implemented as an `.Rmd` report and can be
coordinated by a [Snakemake](https://snakemake.readthedocs.io/en/stable/)
pipeline. The central `.Rds` outputs are created by:

- `differential_expression.Rmd`: Outputs the PEA `.Rds` required by all `_pea.Rmd` scripts.
- `elisa.Rmd`: Outputs the ELISA `.Rds` required by all `_elisa.Rmd` scripts.

## Analysis Components

### PEA-Focused Reports

| Report | Description |
| :--- | :--- |
| **`differential_expression.Rmd`** | QC, contrast configurations, and differential abundance testing. |
| **`assumptions_pea.Rmd`** | Parametric assumption checks for the PEA models (normal errors, homoscedasticity, outliers). |
| **`sensitivity_pea.Rmd`** | Sensitivity of model results to specific samples and model terms (leave-one-out analysis). |
| **`covariates_pea.Rmd`** | Tests whether covariates (e.g., Age, Sex) differ between contrast arms (for use with non-parametric frameworks). |
| **`correlations_pea.Rmd`** | Spearman correlations between NPX values and selected metadata variables. |
| **`annotations_pea.Rmd`** | Human Protein Atlas (HPA) IHC tissue-specific annotations of PEA results. |
| **`functional_enrichment_pea.Rmd`** | Functional enrichment analysis (ORA / GSEA). *Note: Currently under development.* |

### ELISA-Focused Reports

| Report | Description |
| :--- | :--- |
| **`elisa.Rmd`** | QC and differential analysis of ELISA measurements. |
| **`assumptions_elisa.Rmd`** | Parametric assumption checks for the ELISA models. |

---

## Quickstart using Test Data on Linux

### 1. Clone the GitHub repository

```bash
git clone git@github.com:NICHD-BSPC/AnalyzeOlink.git
```

### 2. Create the Conda environment

```bash
cd environment/
conda env create -p ./env --file env.yaml -y
```

### 3. Activate the environment

```bash
conda activate ./env
```

### 4. Generate test data

```bash
# From repo root
cd ../
Rscript data/test_data.R
```

### 5. Run the pipeline

```bash
cd snakemake
snakemake -j 8
```

### 6. View results

Open the `.html` reports generated in the `rmds/` directory in a web browser. Each report links to its exported tables and figures.

---

## Resource Requirements and Runtime

- **Installation time:** < 3 minutes to clone the repo, build the Conda environment, and create the test data.
- **Workflow runtime (test data):** < 3 minutes to run the Snakefile.
- **Hardware requirements:**
  - ≥ 2 GB RAM
  - ≥ 2 CPU cores

## Tested Platforms

- **Operating system:** Linux (Ubuntu LTS; tested on GitHub Actions `ubuntu-latest` runners, as of December 2025: Ubuntu 24.04)
- **Conda distribution:** Miniforge3 (GitHub Actions latest Miniforge3 release, as of December 2025)
- **Conda version:** 25.11.1

---

## User Configuration (for your projects data)

To run the workflow on your own PEA/ELISA data, you must configure the pipeline in two locations.

### 1. Global Settings (`config/config.yaml`)

This is the primary control center for the pipeline. **Open `config/config.yaml` and read the detailed comments inside.**

The config file handles:
* **Input Paths:** Locations for your Olink and optionally ELISA data files.
* **Metadata Requirements:** Specifics on required columns (e.g., `SampleID`, Grouping columns).
* **Pipeline Control:** Which analysis reports to run or skip via Snakemake.
* **Global Parameters:** Significance thresholds (`alpha`), parallelization (`cores`), and enrichment databases.

### 2. Analysis-Specific Settings (`.Rmd` Scripts)

While file paths are handled globally, specific analysis logic is configured directly within the individual R Markdown reports.

Open the `.Rmd` files in `rmds/` and look for the code chunks containing the text, "**USER SETTINGS**".

Key configurations found in `.Rmd` files include:
* **Outlier Removal:** Defining specific sample IDs to exclude based on QC.
* **Contrast Definitions:** Defining the specific statistical comparisons (e.g., `disease_vs_control`) and subsetting the data for each contrast.
* **Covariate Selection:** Specifying which metadata columns to use for correlation or covariate testing.
* **Sensitivity Checks:** Defining which model terms to test in leave-one-out analyses.

---

## Output Summary

Paths defined under `output_paths` in `config/config.yaml` are populated with:
- **XLSX tables**
- **PDF figures**

The `.html` reports in `rmds/` provide clickable links to these outputs.

---

## Dependencies and Package Versions

All package dependencies are managed via Conda in `environment/env.yaml`. Use
the same environment to run both Snakemake and the `.Rmd` scripts (if running
them outside of Snakemake). See the Conda specification file at
`environment/env.yaml` for exact package versions and dependency resolution.

---

## Support

For questions or issues, please reach out to the [NICHD Bioinformatics and Scientific Programming Core](https://www.nichd.nih.gov/about/org/dir/other-facilities/cores/bioinformatics).

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.
See the [LICENSE](LICENSE) file for full terms.
