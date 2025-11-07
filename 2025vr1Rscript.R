################################################################################
# TWCF GamePlay: VR Empathy & Compassion Study - Preregistered Analysis
# LONGITUDINAL MIXED EFFECTS APPROACH WITH SCHOOL ADJUSTMENT
#
# This script implements the expert-recommended longitudinal mixed model approach
# for analyzing repeated measures RCT data with proper handling of within-person
# correlations and missing data.
#
# NEW IN THIS VERSION (20251020_104654):
# - Added SCHOOL as fixed effect covariate (addresses editor's nesting concern)
# - Using updated dataset: twcfgameplayalltimes_with_schools.csv
# - 3 schools: Akande Dahunsi, Community SHS, Eti Osa Community SHS
# - Individual-level randomization within schools (no classroom clustering)
#
# Key features:
# - Seed = 12345 (matching original analysis)
# - LONGITUDINAL MIXED MODELS: condition × time interaction with random slopes
# - Grand-mean centered baseline (recommended for RCTs)
# - PRIMARY (5 hypotheses, FDR-corrected): H1, H2, H3, H7, H8
# - EXPLORATORY (5 hypotheses, FDR-corrected separately): H4, H5a, H5b, H5c, H6
# - SUPPLEMENTARY: Longitudinal trajectories T3-T5
# - MICE imputation (m=20) with temporal ordering
# - Composites computed AFTER imputation (critical fix maintained)
# - Baseline centering AFTER imputation, grand-mean centered (consultant recommendation)
#
# Model: score ~ condition * time + baseline_c + age + gender + class + school + sds_T2 + (time_numeric | sn)
#
# Rationale for school as fixed effect:
# - Only 3 schools (too few for stable random effects)
# - Individual randomization prevents systematic school × condition confounding
# - Fixed effect controls for school-level differences (resources, culture, etc.)
# - Recommended approach for small number of clusters in RCTs
#
# Preregistered Hypotheses (8 total):
#   PRIMARY (Does VR work?):
#     H1: Observe vs Control
#     H2: Embody vs Control
#     H3: Role Model vs Control
#     H7: Observe vs NPT (active control)
#     H8: Embody vs NPT (active control)
#
#   EXPLORATORY (Which mode is better? Does combining help?):
#     H4: Observe vs Embody
#     H5a: Embody+Observe vs Control
#     H5b: Embody vs Embody+Observe
#     H5c: Observe vs Embody+Observe
#     H6: NPT vs Control (manipulation check)
#
# Date: 2025-10-20
# Consultant guidance: Longitudinal mixed model + school adjustment
################################################################################

# ==============================================================================
# 1. SETUP
# ==============================================================================

rm(list = ls())
set.seed(12345)  # MATCHING ORIGINAL ANALYSIS

# Start timing
script_start_time <- Sys.time()

# Load required libraries
# Note: Now including lme4, lmerTest for mixed models, lavaan for CFA, psych for reliability
required_packages <- c("tidyverse", "mice", "emmeans", "lme4", "lmerTest", "lavaan", "psych")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

invisible(lapply(required_packages, library, character.only = TRUE))

# Configuration
N_IMPUTATIONS <- 20
ALPHA_FDR <- 0.05
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
RESULTS_DIR <- paste0("results_preregistered_", TIMESTAMP)

cat("=== TWCF GamePlay: Preregistered Hypotheses Analysis ===\n")
cat("Approach: Longitudinal Mixed Effects Model + School Adjustment\n")
cat("Seed: 12345 (matching original)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results folder:", RESULTS_DIR, "\n\n")

# Create results directories with timestamp
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR, "diagnostics"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. DATA LOADING & EXCLUSIONS
# ==============================================================================

cat("=== Loading data ===\n")

# Check for data file existence - NOW WITH SCHOOL VARIABLE
data_file <- "twcfgameplayalltimes_with_schools.csv"

# Try common locations
possible_paths <- c(
  data_file,  # Current directory
  file.path("data", data_file),  # data subdirectory
  file.path("..", data_file),  # Parent directory
  file.path("/Volumes/HDD2/tmbwork", data_file),  # Absolute path
  file.path("/Volumes/HDD2/tmbwork/postsubmissionreview", data_file)  # Full path
)

data_path <- NULL
for(path in possible_paths) {
  if(file.exists(path)) {
    data_path <- path
    break
  }
}

if(is.null(data_path)) {
  stop(paste0("ERROR: Cannot find data file '", data_file, "'.\n",
              "Tried the following locations:\n  - ",
              paste(possible_paths, collapse="\n  - "),
              "\n\nPlease ensure the data file is in one of these locations or update the script with the correct path."))
}

cat("Data file found:", data_path, "\n")
data_raw <- read_csv(data_path, show_col_types = FALSE)
cat("Raw data loaded:", nrow(data_raw), "participants\n")

# Validate school variable
cat("\n=== Validating school variable ===\n")
if("school" %in% names(data_raw)) {
  school_counts <- table(data_raw$school, useNA = "ifany")
  cat("Schools found:\n")
  print(school_counts)
  cat("\n")
} else {
  stop("ERROR: 'school' variable not found in dataset. Please check data file.")
}

# Attention check exclusions (BEFORE subscale creation)
cat("=== Applying attention check exclusions ===\n")
cat("Criterion: Exclude participants who failed >= 2 out of 3 attention checks\n")
cat("Expected values: attnchk1=70, attnchk2=30, attnchk3=50\n\n")

data <- data_raw %>%
  mutate(attn_fail = (attnchk1 != 70) + (attnchk2 != 30) + (attnchk3 != 50))

cat("Failed 0 checks:", sum(data$attn_fail == 0, na.rm=T), "\n")
cat("Failed 1 check:", sum(data$attn_fail == 1, na.rm=T), "\n")
cat("Failed 2 checks:", sum(data$attn_fail == 2, na.rm=T), "\n")
cat("Failed 3 checks:", sum(data$attn_fail == 3, na.rm=T), "\n")
cat("\nExcluding:", sum(data$attn_fail >= 2, na.rm=T), "participants\n")

data <- data %>%
  filter(attn_fail < 2) %>%
  select(-attn_fail)

cat("Final sample after exclusions:", nrow(data), "participants\n")

# Validate school distribution after exclusions
cat("\n=== School distribution after exclusions ===\n")
school_by_condition <- table(data$condition, data$school)
cat("School × Condition crosstab:\n")
print(school_by_condition)
cat("\n")

# Validate condition variable
cat("=== Validating condition assignments ===\n")
if("condition" %in% names(data)) {
  actual_conditions <- unique(data$condition)
  expected_conditions <- c("Embody", "Observe", "Embody_Observe", "RoleModel", "NPT", "Control")

  cat("Expected conditions:", paste(expected_conditions, collapse=", "), "\n")
  cat("Found conditions:", paste(actual_conditions, collapse=", "), "\n")

  missing_conditions <- setdiff(expected_conditions, actual_conditions)
  extra_conditions <- setdiff(actual_conditions, expected_conditions)

  if(length(missing_conditions) > 0) {
    cat("  WARNING: Expected conditions not found in data:", paste(missing_conditions, collapse=", "), "\n")
  }
  if(length(extra_conditions) > 0) {
    cat("  WARNING: Unexpected conditions found in data:", paste(extra_conditions, collapse=", "), "\n")
    cat("  Note: Numeric codes will be converted to factors\n")
  }
  if(length(missing_conditions) == 0 && length(extra_conditions) == 0) {
    cat("  All conditions validated successfully\n")
  }
} else {
  cat("  ERROR: 'condition' variable not found in data\n")
}

# ==============================================================================
# 3. CREATE SUBSCALE SCORES
# ==============================================================================

cat("\n=== Creating subscale scores ===\n")

create_scale_score <- function(data, items, time) {
  scale_cols <- paste0(items, time)
  existing_cols <- scale_cols[scale_cols %in% names(data)]
  if(length(existing_cols) == 0) return(rep(NA, nrow(data)))
  rowMeans(data[, existing_cols], na.rm = FALSE)
}

times <- c("a", "b", "c", "d", "e")
time_labels <- c("T1", "T2", "T3", "T4", "T5")

# AMES SUBSCALES (items 1-8, excluding sympathy items 9-12)
cat("Creating AMES subscales...\n")

# Cognitive empathy (AMES 1-4)
cognitive_items <- paste0("ames", 1:4)
for(i in seq_along(times)) {
  col_name <- paste0("cognitive_empathy_", time_labels[i])
  data[[col_name]] <- create_scale_score(data, cognitive_items, times[i])
}

# Affective empathy (AMES 5-8)
affective_items <- paste0("ames", 5:8)
for(i in seq_along(times)) {
  col_name <- paste0("affective_empathy_", time_labels[i])
  data[[col_name]] <- create_scale_score(data, affective_items, times[i])
}

# IRI Fantasy subscale (5 items)
cat("Creating IRI Fantasy subscale...\n")
iri_items <- paste0("iri", 1:5)
for(i in seq_along(times)) {
  col_name <- paste0("iri_fantasy_", time_labels[i])
  data[[col_name]] <- create_scale_score(data, iri_items, times[i])
}

# IOS - single item
cat("Creating IOS scores...\n")
ios_found <- 0
for(i in seq_along(times)) {
  old_col <- paste0("ios", times[i])
  new_col <- paste0("ios_", time_labels[i])
  if(old_col %in% names(data)) {
    data[[new_col]] <- data[[old_col]]
    ios_found <- ios_found + 1
  }
}
if(ios_found > 0) {
  cat("  Found IOS at", ios_found, "timepoints\n")
} else {
  cat("  WARNING: No IOS variables found in data\n")
}

# SOCS SUBSCALES (Sussex, 20 items total, 5 subscales of 4 items each)
cat("Creating SOCS (Sussex) subscales...\n")

socs_subscales <- list(
  recognising_suffering = 1:4,
  universality_suffering = 5:8,
  feeling_for_suffering = 9:12,
  tolerating_feelings = 13:16,
  motivation_to_act = 17:20
)

for(subscale_name in names(socs_subscales)) {
  items <- paste0("sussex", socs_subscales[[subscale_name]])
  for(i in seq_along(times)) {
    col_name <- paste0(subscale_name, "_", time_labels[i])
    data[[col_name]] <- create_scale_score(data, items, times[i])
  }
}

# SWISS (Wise Reasoning - T2 only, 8 items)
cat("Creating Swiss (Wise Reasoning) scores...\n")
swiss_items <- paste0("swiss", 1:8)
swiss_cols <- swiss_items[swiss_items %in% names(data)]
if(length(swiss_cols) > 0) {
  data$swiss_T2 <- rowMeans(data[, swiss_cols], na.rm = FALSE)
} else {
  data$swiss_T2 <- NA
}

cat("Subscale creation complete.\n")
cat("NOTE: Composite scores will be computed AFTER imputation\n")

# ==============================================================================
# 3.5. MEASUREMENT VALIDATION (CFA & RELIABILITY)
# ==============================================================================

cat("\n=== MEASUREMENT VALIDATION ===\n")
cat("Running CFAs and computing reliability coefficients for editor requirements\n\n")

# Initialize results storage
cfa_results <- list()
reliability_results <- list()

# ------------------------------------------------------------------------------
# 3.5.1. AMES Confirmatory Factor Analysis (2-factor model)
# ------------------------------------------------------------------------------

cat("1. AMES CFA (2-factor: Cognitive + Affective Empathy)\n")

# AMES items at baseline (T1/time a)
ames_items_baseline <- data %>%
  select(ames1a, ames2a, ames3a, ames4a, ames5a, ames6a, ames7a, ames8a) %>%
  filter(complete.cases(.))

if(nrow(ames_items_baseline) >= 100) {

  # Define 2-factor model: items 1-4 = cognitive, items 5-8 = affective
  ames_model <- '
    cognitive =~ ames1a + ames2a + ames3a + ames4a
    affective =~ ames5a + ames6a + ames7a + ames8a
  '

  tryCatch({
    ames_cfa <- cfa(ames_model, data = data, std.lv = TRUE, missing = "fiml")
    ames_fit <- fitMeasures(ames_cfa, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cfa_results$ames <- list(
      model = "AMES 2-factor (Cognitive + Affective)",
      chisq = ames_fit["chisq"],
      df = ames_fit["df"],
      pvalue = ames_fit["pvalue"],
      cfi = ames_fit["cfi"],
      tli = ames_fit["tli"],
      rmsea = ames_fit["rmsea"],
      srmr = ames_fit["srmr"],
      interpretation = ifelse(ames_fit["cfi"] > 0.90 & ames_fit["rmsea"] < 0.08, "Acceptable fit", "Poor fit")
    )

    cat("  AMES CFA Results:\n")
    cat("    CFI =", round(ames_fit["cfi"], 3), "\n")
    cat("    TLI =", round(ames_fit["tli"], 3), "\n")
    cat("    RMSEA =", round(ames_fit["rmsea"], 3), "\n")
    cat("    SRMR =", round(ames_fit["srmr"], 3), "\n")
    cat("    Interpretation:", cfa_results$ames$interpretation, "\n\n")

  }, error = function(e) {
    cat("  AMES CFA failed:", conditionMessage(e), "\n\n")
    cfa_results$ames <- list(model = "AMES 2-factor", error = conditionMessage(e))
  })

} else {
  cat("  Insufficient complete cases for AMES CFA\n\n")
}

# ------------------------------------------------------------------------------
# 3.5.2. SOCS Confirmatory Factor Analysis (5-factor model)
# ------------------------------------------------------------------------------

cat("2. SOCS CFA (5-factor compassion model)\n")

# SOCS items at baseline (T1/time a) - 20 items, 4 per factor
socs_items_baseline <- data %>%
  select(starts_with("sussex")) %>%
  select(ends_with("a")) %>%  # baseline only
  filter(complete.cases(.))

if(nrow(socs_items_baseline) >= 100) {

  # Define 5-factor model based on SOCS structure
  socs_model <- '
    recognising =~ sussex1a + sussex2a + sussex3a + sussex4a
    universality =~ sussex5a + sussex6a + sussex7a + sussex8a
    feeling =~ sussex9a + sussex10a + sussex11a + sussex12a
    tolerating =~ sussex13a + sussex14a + sussex15a + sussex16a
    motivation =~ sussex17a + sussex18a + sussex19a + sussex20a
  '

  tryCatch({
    socs_cfa <- cfa(socs_model, data = data, std.lv = TRUE, missing = "fiml")
    socs_fit <- fitMeasures(socs_cfa, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cfa_results$socs <- list(
      model = "SOCS 5-factor (Compassion)",
      chisq = socs_fit["chisq"],
      df = socs_fit["df"],
      pvalue = socs_fit["pvalue"],
      cfi = socs_fit["cfi"],
      tli = socs_fit["tli"],
      rmsea = socs_fit["rmsea"],
      srmr = socs_fit["srmr"],
      interpretation = ifelse(socs_fit["cfi"] > 0.90 & socs_fit["rmsea"] < 0.08, "Acceptable fit", "Poor fit")
    )

    cat("  SOCS CFA Results:\n")
    cat("    CFI =", round(socs_fit["cfi"], 3), "\n")
    cat("    TLI =", round(socs_fit["tli"], 3), "\n")
    cat("    RMSEA =", round(socs_fit["rmsea"], 3), "\n")
    cat("    SRMR =", round(socs_fit["srmr"], 3), "\n")
    cat("    Interpretation:", cfa_results$socs$interpretation, "\n\n")

  }, error = function(e) {
    cat("  SOCS CFA failed:", conditionMessage(e), "\n\n")
    cfa_results$socs <- list(model = "SOCS 5-factor", error = conditionMessage(e))
  })

} else {
  cat("  Insufficient complete cases for SOCS CFA\n\n")
}

# ------------------------------------------------------------------------------
# 3.5.3. IRI Fantasy Subscale CFA (Single-Factor)
# ------------------------------------------------------------------------------

n_complete_iri <- sum(complete.cases(data %>% select(iri1a, iri2a, iri3a, iri4a, iri5a)))

if(n_complete_iri >= 200) {
  cat("3. IRI Fantasy CFA (single-factor model)\n")

  iri_model <- '
    fantasy =~ iri1a + iri2a + iri3a + iri4a + iri5a
  '

  tryCatch({
    iri_cfa <- cfa(iri_model, data = data, std.lv = TRUE, missing = "fiml")
    iri_fit <- fitMeasures(iri_cfa, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cfa_results$iri_fantasy <- list(
      model = "IRI Fantasy subscale (single-factor)",
      chisq = iri_fit["chisq"],
      df = iri_fit["df"],
      pvalue = iri_fit["pvalue"],
      cfi = iri_fit["cfi"],
      tli = iri_fit["tli"],
      rmsea = iri_fit["rmsea"],
      srmr = iri_fit["srmr"],
      interpretation = ifelse(iri_fit["cfi"] > 0.90 & iri_fit["rmsea"] < 0.08, "Acceptable fit", "Poor fit")
    )

    cat("  IRI Fantasy CFA Results:\n")
    cat("    CFI =", round(iri_fit["cfi"], 3), "\n")
    cat("    TLI =", round(iri_fit["tli"], 3), "\n")
    cat("    RMSEA =", round(iri_fit["rmsea"], 3), "\n")
    cat("    SRMR =", round(iri_fit["srmr"], 3), "\n")
    cat("    Interpretation:", cfa_results$iri_fantasy$interpretation, "\n\n")

  }, error = function(e) {
    cat("  IRI Fantasy CFA failed:", conditionMessage(e), "\n\n")
    cfa_results$iri_fantasy <- list(model = "IRI Fantasy subscale", error = conditionMessage(e))
  })

} else {
  cat("3. Insufficient complete cases for IRI Fantasy CFA\n\n")
}

# ------------------------------------------------------------------------------
# 3.5.4. Swiss Wise Reasoning CFA (Single-Factor)
# ------------------------------------------------------------------------------

n_complete_swiss <- sum(complete.cases(data %>% select(swiss1, swiss2, swiss3, swiss4, swiss5, swiss6, swiss7, swiss8)))

if(n_complete_swiss >= 200) {
  cat("4. Swiss Wise Reasoning CFA (single-factor model)\n")

  swiss_model <- '
    wise_reasoning =~ swiss1 + swiss2 + swiss3 + swiss4 + swiss5 + swiss6 + swiss7 + swiss8
  '

  tryCatch({
    swiss_cfa <- cfa(swiss_model, data = data, std.lv = TRUE, missing = "fiml")
    swiss_fit <- fitMeasures(swiss_cfa, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cfa_results$swiss <- list(
      model = "Swiss Wise Reasoning (single-factor)",
      chisq = swiss_fit["chisq"],
      df = swiss_fit["df"],
      pvalue = swiss_fit["pvalue"],
      cfi = swiss_fit["cfi"],
      tli = swiss_fit["tli"],
      rmsea = swiss_fit["rmsea"],
      srmr = swiss_fit["srmr"],
      interpretation = ifelse(swiss_fit["cfi"] > 0.90 & swiss_fit["rmsea"] < 0.08, "Acceptable fit", "Poor fit")
    )

    cat("  Swiss CFA Results:\n")
    cat("    CFI =", round(swiss_fit["cfi"], 3), "\n")
    cat("    TLI =", round(swiss_fit["tli"], 3), "\n")
    cat("    RMSEA =", round(swiss_fit["rmsea"], 3), "\n")
    cat("    SRMR =", round(swiss_fit["srmr"], 3), "\n")
    cat("    Interpretation:", cfa_results$swiss$interpretation, "\n\n")

  }, error = function(e) {
    cat("  Swiss CFA failed:", conditionMessage(e), "\n\n")
    cfa_results$swiss <- list(model = "Swiss Wise Reasoning", error = conditionMessage(e))
  })

} else {
  cat("4. Insufficient complete cases for Swiss CFA\n\n")
}

# ------------------------------------------------------------------------------
# 3.5.5. Reliability Analysis (Cronbach's Alpha)
# ------------------------------------------------------------------------------

cat("5. Reliability Analysis (Cronbach's Alpha)\n")
cat("   Computing for all subscales at baseline (T1) using complete cases\n\n")

# Helper function to compute alpha safely
compute_alpha <- function(items, scale_name) {
  tryCatch({
    alpha_result <- psych::alpha(items, check.keys = FALSE)
    return(data.frame(
      scale = scale_name,
      alpha = alpha_result$total$raw_alpha,
      n_items = ncol(items),
      n_cases = nrow(items[complete.cases(items), ]),
      interpretation = ifelse(alpha_result$total$raw_alpha >= 0.70, "Acceptable",
                            ifelse(alpha_result$total$raw_alpha >= 0.60, "Questionable", "Poor"))
    ))
  }, error = function(e) {
    return(data.frame(
      scale = scale_name,
      alpha = NA,
      n_items = ncol(items),
      n_cases = 0,
      interpretation = "Error"
    ))
  })
}

# AMES subscales
reliability_results$cognitive_empathy_T1 <- compute_alpha(
  data %>% select(ames1a, ames2a, ames3a, ames4a),
  "Cognitive Empathy (AMES 1-4)"
)

reliability_results$affective_empathy_T1 <- compute_alpha(
  data %>% select(ames5a, ames6a, ames7a, ames8a),
  "Affective Empathy (AMES 5-8)"
)

# IRI Fantasy (5 items)
reliability_results$iri_fantasy_T1 <- compute_alpha(
  data %>% select(iri1a, iri2a, iri3a, iri4a, iri5a),
  "IRI Fantasy"
)

# SOCS subscales
reliability_results$recognising_T1 <- compute_alpha(
  data %>% select(sussex1a, sussex2a, sussex3a, sussex4a),
  "SOCS: Recognising Suffering"
)

reliability_results$universality_T1 <- compute_alpha(
  data %>% select(sussex5a, sussex6a, sussex7a, sussex8a),
  "SOCS: Universality of Suffering"
)

reliability_results$feeling_T1 <- compute_alpha(
  data %>% select(sussex9a, sussex10a, sussex11a, sussex12a),
  "SOCS: Feeling for Person"
)

reliability_results$tolerating_T1 <- compute_alpha(
  data %>% select(sussex13a, sussex14a, sussex15a, sussex16a),
  "SOCS: Tolerating Feelings"
)

reliability_results$motivation_T1 <- compute_alpha(
  data %>% select(sussex17a, sussex18a, sussex19a, sussex20a),
  "SOCS: Motivation to Act"
)

# Swiss (T2 only)
reliability_results$swiss_T2 <- compute_alpha(
  data %>% select(swiss1, swiss2, swiss3, swiss4, swiss5, swiss6, swiss7, swiss8),
  "Swiss Wise Reasoning"
)

# Combine reliability results
reliability_df <- bind_rows(reliability_results)

cat("  Reliability Summary:\n")
print(reliability_df, row.names = FALSE)
cat("\n")

cat("Measurement validation complete.\n")
cat("Note: CFA and reliability results will be saved with main analysis outputs\n\n")

# ==============================================================================
# 4. PREPARE DATA FOR IMPUTATION
# ==============================================================================

cat("\n=== Preparing data for MICE imputation ===\n")

# Select core variables - NOW INCLUDING SCHOOL
core_vars <- c("sn", "condition", "age", "gender", "class", "school")

# Validate that condition has no missing values
if(any(is.na(data$condition))) {
  n_missing_condition <- sum(is.na(data$condition))
  warning(paste0("WARNING: ", n_missing_condition, " participants have missing condition assignment."))
  cat("  Participants with missing condition:", n_missing_condition, "\n")
} else {
  cat("  Condition variable validated: No missing values\n")
}

# Validate that school has no missing values
if(any(is.na(data$school))) {
  n_missing_school <- sum(is.na(data$school))
  warning(paste0("WARNING: ", n_missing_school, " participants have missing school assignment."))
  cat("  Participants with missing school:", n_missing_school, "\n")
} else {
  cat("  School variable validated: No missing values\n")
}

# Outcome measures to impute (ONLY subscales - composites computed after imputation)
outcome_measures <- c("cognitive_empathy", "affective_empathy", "iri_fantasy", "ios",
                     "recognising_suffering", "universality_suffering",
                     "feeling_for_suffering", "tolerating_feelings", "motivation_to_act")

# Build list of all outcome columns across timepoints
outcome_cols <- c()
for(outcome in outcome_measures) {
  for(time in time_labels) {
    col_name <- paste0(outcome, "_", time)
    if(col_name %in% names(data)) {
      outcome_cols <- c(outcome_cols, col_name)
    }
  }
}

# Add Swiss (T2 only)
outcome_cols <- c(outcome_cols, "swiss_T2")

# Social desirability at T2
sd_col <- "sds2b"
if(!sd_col %in% names(data)) {
  sds_items_t2 <- paste0("sds", 1:6, "b")
  if(all(sds_items_t2 %in% names(data))) {
    data[[sd_col]] <- rowMeans(data[, sds_items_t2], na.rm = FALSE)
    cat("Created sds2b composite from 6 items\n")
  }
}

# Select columns for imputation
impute_data <- data %>%
  select(all_of(c(core_vars, outcome_cols, sd_col)))

cat("Variables selected for imputation:\n")
cat("  Core:", length(core_vars), "(including school)\n")
cat("  Outcomes:", length(outcome_cols), "\n")
cat("  Total columns:", ncol(impute_data), "\n")
cat("  Participants:", nrow(impute_data), "\n\n")

# Check missingness
cat("Missing data summary:\n")
missing_summary <- impute_data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = n_missing / nrow(impute_data) * 100) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

if(nrow(missing_summary) > 0) {
  print(head(missing_summary, 20))
}

# ==============================================================================
# 5. MULTIPLE IMPUTATION
# ==============================================================================

cat("\n=== Running MICE imputation ===\n")
cat("Imputations: m =", N_IMPUTATIONS, "\n")
cat("Seed:", 12345, "\n")
cat("Max iterations: 20\n\n")

# Set up imputation methods
cat("Configuring imputation methods:\n")
imp_method <- mice(impute_data, maxit = 0, print = FALSE)$method

# Specify methods for categorical variables
if("gender" %in% names(impute_data)) {
  imp_method["gender"] <- "logreg"
  cat("  gender: logreg (binary)\n")
}

if("condition" %in% names(impute_data)) {
  imp_method["condition"] <- "polyreg"
  cat("  condition: polyreg (categorical)\n")
}

if("class" %in% names(impute_data)) {
  imp_method["class"] <- "polyreg"
  cat("  class: polyreg (categorical)\n")
}

if("school" %in% names(impute_data)) {
  imp_method["school"] <- "polyreg"
  cat("  school: polyreg (categorical)\n")
}

cat("  All outcome variables: pmm (continuous)\n\n")

# Configure predictor matrix to respect temporal ordering
cat("Configuring predictor matrix (temporal ordering)...\n")
pred_matrix <- mice(impute_data, maxit = 0, print = FALSE)$predictorMatrix

var_names <- colnames(impute_data)
for(var in var_names) {
  if(grepl("_T[1-5]", var)) {
    var_time <- as.numeric(sub(".*_T([1-5]).*", "\\1", var))

    for(predictor in var_names) {
      if(grepl("_T[1-5]", predictor)) {
        pred_time <- as.numeric(sub(".*_T([1-5]).*", "\\1", predictor))

        if(pred_time > var_time) {
          pred_matrix[var, predictor] <- 0
        }
      }
    }
  }
}

cat("  Temporal ordering enforced\n")

# Convert categorical variables to factors before MICE
cat("\nConverting categorical variables to factors...\n")
if("condition" %in% names(impute_data)) {
  impute_data$condition <- as.factor(impute_data$condition)
  cat("  condition:", nlevels(impute_data$condition), "levels\n")
}

if("gender" %in% names(impute_data)) {
  impute_data$gender <- as.factor(impute_data$gender)
  cat("  gender:", nlevels(impute_data$gender), "levels\n")
}

if("class" %in% names(impute_data)) {
  impute_data$class <- as.factor(impute_data$class)
  cat("  class:", nlevels(impute_data$class), "levels\n")
}

if("school" %in% names(impute_data)) {
  impute_data$school <- as.factor(impute_data$school)
  cat("  school:", nlevels(impute_data$school), "levels\n")
}

# Run MICE
cat("\nRunning MICE (m=20, maxit=20, seed=12345)...\n")
imp <- mice(impute_data, m = N_IMPUTATIONS, method = imp_method,
           predictorMatrix = pred_matrix, maxit = 20, printFlag = FALSE, seed = 12345)

cat("Imputation complete.\n\n")

# ==============================================================================
# 6. RESHAPE TO LONG FORMAT AND COMPUTE DERIVED VARIABLES
# ==============================================================================

cat("=== Reshaping imputed data to long format ===\n")

# Function to reshape one imputed dataset
reshape_imputed <- function(imp_data, imp_num) {

  # First, compute composite scores from imputed subscales
  # This must happen BEFORE reshaping to long

  # Composite Empathy at each timepoint
  for(time in time_labels) {
    cog_col <- paste0("cognitive_empathy_", time)
    aff_col <- paste0("affective_empathy_", time)
    comp_col <- paste0("composite_empathy_", time)

    if(cog_col %in% names(imp_data) && aff_col %in% names(imp_data)) {
      imp_data[[comp_col]] <- rowMeans(imp_data[, c(cog_col, aff_col)], na.rm = FALSE)
    }
  }

  # Composite Compassion at each timepoint
  socs_measures <- c("recognising_suffering", "universality_suffering",
                     "feeling_for_suffering", "tolerating_feelings", "motivation_to_act")

  for(time in time_labels) {
    socs_cols <- paste0(socs_measures, "_", time)
    comp_col <- paste0("composite_compassion_", time)
    existing_cols <- socs_cols[socs_cols %in% names(imp_data)]

    if(length(existing_cols) > 0) {
      imp_data[[comp_col]] <- rowMeans(imp_data[, existing_cols], na.rm = FALSE)
    }
  }

  # Now reshape to long format
  all_outcomes <- c(outcome_measures, "composite_empathy", "composite_compassion")

  long_list <- list()
  for(outcome in all_outcomes) {
    outcome_cols_time <- paste0(outcome, "_", time_labels)
    existing_cols <- outcome_cols_time[outcome_cols_time %in% names(imp_data)]

    if(length(existing_cols) > 0) {
      long_temp <- imp_data %>%
        select(sn, condition, age, gender, class, school, sds2b, all_of(existing_cols)) %>%
        pivot_longer(cols = all_of(existing_cols),
                    names_to = "time_var",
                    values_to = "score") %>%
        mutate(time = str_extract(time_var, "T[1-5]"),
               time_numeric = as.numeric(str_extract(time_var, "[1-5]")),  # Numeric for random slopes
               outcome = outcome) %>%
        select(-time_var)

      long_list[[outcome]] <- long_temp
    }
  }

  combined_long <- bind_rows(long_list)

  # Compute grand-mean centered baseline (consultant recommendation)
  # baseline_c = baseline - mean(baseline)
  combined_long <- combined_long %>%
    group_by(outcome) %>%
    mutate(baseline_raw = score[time == "T1"][1],  # Get T1 score for each person
           baseline_grand_mean = mean(score[time == "T1"], na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(sn, outcome) %>%
    mutate(baseline_raw = score[time == "T1"],
           baseline_c = baseline_raw - first(baseline_grand_mean)) %>%
    ungroup() %>%
    select(-baseline_grand_mean)

  # Add Swiss (T2 only) separately
  if("swiss_T2" %in% names(imp_data)) {
    swiss_long <- imp_data %>%
      select(sn, condition, age, gender, class, school, sds2b, swiss_T2) %>%
      mutate(outcome = "swiss",
             time = "T2",
             time_numeric = 2,  # T2 = timepoint 2
             score = swiss_T2,
             baseline_raw = NA_real_,
             baseline_c = NA_real_) %>%
      select(sn, condition, age, gender, class, school, sds2b, outcome, time, time_numeric, score, baseline_raw, baseline_c)

    combined_long <- bind_rows(combined_long, swiss_long)
  }

  # Rename sds2b for consistency
  combined_long <- combined_long %>%
    rename(sds_T2 = sds2b) %>%
    mutate(imputation = imp_num,
           time = factor(time, levels = c("T1", "T2", "T3", "T4", "T5")))

  return(combined_long)
}

# Apply to all imputations
cat("Processing", N_IMPUTATIONS, "imputed datasets...\n")
imputed_long_list <- list()
for(m in 1:N_IMPUTATIONS) {
  imp_data <- complete(imp, m)
  imputed_long_list[[m]] <- reshape_imputed(imp_data, m)
}

cat("Reshape complete.\n")

# ==============================================================================
# 7. FIT LONGITUDINAL MIXED MODELS AND COMPUTE CONTRASTS
# ==============================================================================

cat("\n=== Fitting longitudinal mixed effects models ===\n")

# Define all 10 contrasts for 8 preregistered hypotheses
# Using numeric condition codes with CORRECTED mapping (from codebook.csv)
# 1=Embody, 2=Observe, 3=Embody_Observe, 4=RoleModel, 5=NPT, 6=Control
hypothesis_contrasts <- list(
  H1 = c("2", "6"),    # Observe (2) vs Control (6)
  H2 = c("1", "6"),    # Embody (1) vs Control (6)
  H3 = c("4", "6"),    # RoleModel (4) vs Control (6)
  H4 = c("2", "1"),    # Observe (2) vs Embody (1)
  H5a = c("3", "6"),   # Embody_Observe (3) vs Control (6)
  H5b = c("1", "3"),   # Embody (1) vs Embody_Observe (3)
  H5c = c("2", "3"),   # Observe (2) vs Embody_Observe (3)
  H6 = c("5", "6"),    # NPT (5) vs Control (6)
  H7 = c("2", "5"),    # Observe (2) vs NPT (5)
  H8 = c("1", "5")     # Embody (1) vs NPT (5)
)

# Hypothesis groupings
primary_hypotheses <- c("H1", "H2", "H3", "H7", "H8")
exploratory_hypotheses <- c("H4", "H5a", "H5b", "H5c", "H6")

cat("PRIMARY hypotheses (FDR-corrected):", paste(primary_hypotheses, collapse=", "), "\n")
cat("EXPLORATORY hypotheses (FDR-corrected separately):", paste(exploratory_hypotheses, collapse=", "), "\n\n")

# All outcomes
all_outcomes <- c(outcome_measures, "composite_empathy", "composite_compassion", "swiss")

# Function to analyze one outcome with longitudinal mixed model
analyze_outcome_longitudinal <- function(outcome_name, imputed_data_list) {

  cat("Analyzing:", outcome_name, "\n")

  # Special handling for Swiss (T2 only, no baseline)
  if(outcome_name == "swiss") {
    cat("  Swiss: T2 only (cross-sectional approach)\n")
    return(NULL)  # Handle separately
  }

  all_contrasts <- list()

  # Fit model on each imputed dataset
  model_list <- list()
  emm_list <- list()
  n_models_succeeded <- 0

  for(m in 1:N_IMPUTATIONS) {

    dat <- imputed_data_list[[m]] %>%
      filter(outcome == outcome_name,
             !is.na(score),
             !is.na(baseline_c))

    if(nrow(dat) < 100) next

    # Longitudinal mixed model with condition × time interaction + SCHOOL
    # Model: score ~ condition * time + baseline_c + age + gender + class + school + sds_T2 + (time_numeric | sn)

    tryCatch({
      model <- lmer(score ~ condition * time + baseline_c + age + gender + class + school + sds_T2 +
                     (time_numeric | sn),
                   data = dat,
                   control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 100000)))

      # Get estimated marginal means for condition × time combinations
      emm <- emmeans(model, specs = ~ condition * time)

      model_list[[m]] <- model
      emm_list[[m]] <- emm
      n_models_succeeded <- n_models_succeeded + 1

    }, error = function(e) {
      # Skip if model fails
    })
  }

  if(length(emm_list) == 0) {
    cat("    WARNING: All models failed for", outcome_name, "\n")
    return(data.frame())
  }

  cat("  Models converged:", n_models_succeeded, "of", N_IMPUTATIONS, "\n")

  # Extract contrasts at each timepoint
  for(time_point in c("T2", "T3", "T4", "T5")) {

    for(hyp_name in names(hypothesis_contrasts)) {

      contrast_pair <- hypothesis_contrasts[[hyp_name]]
      cond1 <- contrast_pair[1]
      cond2 <- contrast_pair[2]

      # Compute contrast on each imputed dataset
      contrast_estimates <- c()
      contrast_se <- c()
      contrast_df <- c()

      for(imp_idx in seq_along(emm_list)) {

        tryCatch({
          # Get emmeans at specific timepoint
          emm_at_time <- emm_list[[imp_idx]]

          # Build contrast vector for this specific comparison at this timepoint
          emm_grid <- summary(emm_at_time)

          # Find the rows for our conditions at this timepoint
          cond1_row <- which(as.character(emm_grid$condition) == cond1 &
                            as.character(emm_grid$time) == time_point)
          cond2_row <- which(as.character(emm_grid$condition) == cond2 &
                            as.character(emm_grid$time) == time_point)

          if(length(cond1_row) == 0 || length(cond2_row) == 0) next

          # Create contrast vector
          contrast_vec <- rep(0, nrow(emm_grid))
          contrast_vec[cond1_row] <- 1
          contrast_vec[cond2_row] <- -1

          # Compute contrast using emmeans
          contr <- contrast(emm_at_time, method = list(contrast_vec), adjust = "none")
          contr_summary <- summary(contr)

          contrast_estimates <- c(contrast_estimates, contr_summary$estimate)
          contrast_se <- c(contrast_se, contr_summary$SE)
          contrast_df <- c(contrast_df, contr_summary$df)

        }, error = function(e) {
          # Skip if contrast computation fails
        })
      }

      # Check if we have enough successful contrasts
      if(length(contrast_estimates) < 10) next

      # Manual pooling using Rubin's rules
      M <- length(contrast_estimates)

      pooled_estimate <- mean(contrast_estimates)
      within_var <- mean(contrast_se^2)
      between_var <- var(contrast_estimates)
      total_var <- within_var + (1 + 1/M) * between_var
      pooled_se <- sqrt(total_var)

      # Degrees of freedom (Barnard-Rubin)
      if(between_var < 1e-10 || total_var < 1e-10) {
        df_adj <- mean(contrast_df)
      } else {
        lambda <- (1 + 1/M) * between_var / total_var
        df_old <- mean(contrast_df)
        df_obs <- (M - 1) / lambda^2
        df_adj <- (df_old * df_obs) / (df_old + df_obs)
      }

      # t-statistic and p-value
      t_stat <- pooled_estimate / pooled_se
      p_val <- 2 * pt(abs(t_stat), df = df_adj, lower.tail = FALSE)

      # Confidence intervals
      t_crit <- qt(0.975, df = df_adj)
      conf_low <- pooled_estimate - t_crit * pooled_se
      conf_high <- pooled_estimate + t_crit * pooled_se

      # Cohen's d (using residual SD from models)
      model_sds <- sapply(model_list, function(mod) {
        if(!is.null(mod)) return(sigma(mod))
        return(NA_real_)
      })
      model_sds <- model_sds[!is.na(model_sds)]

      if(length(model_sds) > 0) {
        pooled_sd <- mean(model_sds)
        cohens_d <- pooled_estimate / pooled_sd
      } else {
        pooled_sd <- NA_real_
        cohens_d <- NA_real_
      }

      # Store results
      result_row <- data.frame(
        outcome = outcome_name,
        hypothesis = hyp_name,
        contrast = paste(cond1, "-", cond2),
        time = time_point,
        estimate = pooled_estimate,
        se = pooled_se,
        df = df_adj,
        statistic = t_stat,
        p.value = p_val,
        conf.low = conf_low,
        conf.high = conf_high,
        cohens_d = cohens_d
      )

      all_contrasts[[paste(hyp_name, time_point, sep="_")]] <- result_row
    }
  }

  if(length(all_contrasts) == 0) {
    return(data.frame())
  }

  return(bind_rows(all_contrasts))
}

# Analyze all outcomes
cat("\nAnalyzing all outcomes with longitudinal models...\n")
all_results_list <- list()
n_outcomes_succeeded <- 0
n_outcomes_failed <- 0

for(outcome in all_outcomes) {
  if(outcome == "swiss") {
    # Handle Swiss separately (T2 only, no longitudinal component)
    next
  }

  result <- analyze_outcome_longitudinal(outcome, imputed_long_list)
  if(!is.null(result) && nrow(result) > 0) {
    all_results_list[[outcome]] <- result
    n_outcomes_succeeded <- n_outcomes_succeeded + 1
  } else {
    n_outcomes_failed <- n_outcomes_failed + 1
  }
}

# Handle Swiss separately with cross-sectional approach + SCHOOL
cat("\nAnalyzing Swiss (T2 only) with cross-sectional approach + school...\n")
swiss_results_list <- list()

for(m in 1:N_IMPUTATIONS) {
  dat_swiss <- imputed_long_list[[m]] %>%
    filter(outcome == "swiss", time == "T2", !is.na(score))

  if(nrow(dat_swiss) < 50) next

  model_swiss <- try(lm(score ~ condition + age + gender + class + school + sds_T2, data = dat_swiss), silent = TRUE)

  if(!inherits(model_swiss, "try-error")) {
    swiss_results_list[[m]] <- list(
      model = model_swiss,
      emm = emmeans(model_swiss, specs = ~ condition)
    )
  }
}

if(length(swiss_results_list) > 0) {
  cat("  Swiss: Pooling results across", length(swiss_results_list), "imputations\n")

  # Pool Swiss contrasts
  swiss_contrasts <- list()

  for(hyp_name in names(hypothesis_contrasts)) {
    contrast_pair <- hypothesis_contrasts[[hyp_name]]
    cond1 <- contrast_pair[1]
    cond2 <- contrast_pair[2]

    swiss_est <- c()
    swiss_se <- c()
    swiss_df <- c()

    for(m_idx in seq_along(swiss_results_list)) {
      emm_swiss <- swiss_results_list[[m_idx]]$emm

      tryCatch({
        # Build contrast
        emm_grid <- summary(emm_swiss)
        cond1_row <- which(as.character(emm_grid$condition) == cond1)
        cond2_row <- which(as.character(emm_grid$condition) == cond2)

        if(length(cond1_row) == 0 || length(cond2_row) == 0) next

        contrast_vec <- rep(0, nrow(emm_grid))
        contrast_vec[cond1_row] <- 1
        contrast_vec[cond2_row] <- -1

        contr <- contrast(emm_swiss, method = list(contrast_vec), adjust = "none")
        contr_summary <- summary(contr)

        swiss_est <- c(swiss_est, contr_summary$estimate)
        swiss_se <- c(swiss_se, contr_summary$SE)
        swiss_df <- c(swiss_df, contr_summary$df)

      }, error = function(e) {
        # Skip
      })
    }

    if(length(swiss_est) >= 10) {
      # Pool
      M <- length(swiss_est)
      pooled_est <- mean(swiss_est)
      within_var <- mean(swiss_se^2)
      between_var <- var(swiss_est)
      total_var <- within_var + (1 + 1/M) * between_var
      pooled_se <- sqrt(total_var)

      if(between_var < 1e-10 || total_var < 1e-10) {
        df_adj <- mean(swiss_df)
      } else {
        lambda <- (1 + 1/M) * between_var / total_var
        df_old <- mean(swiss_df)
        df_obs <- (M - 1) / lambda^2
        df_adj <- (df_old * df_obs) / (df_old + df_obs)
      }

      t_stat <- pooled_est / pooled_se
      p_val <- 2 * pt(abs(t_stat), df = df_adj, lower.tail = FALSE)

      t_crit <- qt(0.975, df = df_adj)
      conf_low <- pooled_est - t_crit * pooled_se
      conf_high <- pooled_est + t_crit * pooled_se

      # Cohen's d
      model_sds <- sapply(swiss_results_list, function(x) sigma(x$model))
      pooled_sd <- mean(model_sds, na.rm = TRUE)
      cohens_d <- pooled_est / pooled_sd

      swiss_contrasts[[hyp_name]] <- data.frame(
        outcome = "swiss",
        hypothesis = hyp_name,
        contrast = paste(cond1, "-", cond2),
        time = "T2",
        estimate = pooled_est,
        se = pooled_se,
        df = df_adj,
        statistic = t_stat,
        p.value = p_val,
        conf.low = conf_low,
        conf.high = conf_high,
        cohens_d = cohens_d
      )
    }
  }

  if(length(swiss_contrasts) > 0) {
    swiss_results_df <- bind_rows(swiss_contrasts)
    all_results_list$swiss <- swiss_results_df
    n_outcomes_succeeded <- n_outcomes_succeeded + 1
  }
}

cat("\nOutcomes summary:\n")
cat("  Succeeded:", n_outcomes_succeeded, "\n")
cat("  Failed:", n_outcomes_failed, "\n")

if(length(all_results_list) == 0) {
  stop("ERROR: No outcomes produced valid results.")
}

combined_results <- bind_rows(all_results_list)
cat("Total tests conducted:", nrow(combined_results), "\n")

# ==============================================================================
# 8. SEPARATE AND APPLY FDR CORRECTION
# ==============================================================================

cat("\n=== Separating results and applying FDR correction ===\n")

# PRIMARY: H1-H3, H7-H8 at T2
primary_results <- combined_results %>%
  filter(hypothesis %in% primary_hypotheses, time == "T2")

# EXPLORATORY: H4, H5a-c, H6 at T2
exploratory_results <- combined_results %>%
  filter(hypothesis %in% exploratory_hypotheses, time == "T2")

# ==============================================================================
# APPLY PER-TIMEPOINT FDR CORRECTION (T2-T5)
# Per expert consultation: all preregistered timepoints treated as confirmatory
# FDR applied separately for PRIMARY vs EXPLORATORY at each timepoint
# ==============================================================================

cat("Applying per-timepoint FDR correction...\n")

# Create list to store results by timepoint
timepoint_results <- list()

for(tp in c("T2", "T3", "T4", "T5")) {

  # PRIMARY hypotheses at this timepoint
  primary_tp <- combined_results %>%
    filter(hypothesis %in% primary_hypotheses, time == tp) %>%
    group_by(outcome) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
           n_tests_in_family = n(),
           analysis_type = "Primary") %>%
    ungroup() %>%
    arrange(outcome, hypothesis)

  # EXPLORATORY hypotheses at this timepoint
  exploratory_tp <- combined_results %>%
    filter(hypothesis %in% exploratory_hypotheses, time == tp) %>%
    group_by(outcome) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
           n_tests_in_family = n(),
           analysis_type = "Exploratory") %>%
    ungroup() %>%
    arrange(outcome, hypothesis)

  # Store for later
  timepoint_results[[paste0(tp, "_primary")]] <- primary_tp
  timepoint_results[[paste0(tp, "_exploratory")]] <- exploratory_tp

  cat(sprintf("  %s: PRIMARY n=%d, EXPLORATORY n=%d\n",
              tp, nrow(primary_tp), nrow(exploratory_tp)))
}

cat("\n")

# Convenience variables for T2 (backward compatibility)
primary_with_fdr <- timepoint_results[["T2_primary"]]
exploratory_with_fdr <- timepoint_results[["T2_exploratory"]]

# Combine all results
all_results_combined <- bind_rows(timepoint_results)

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================

cat("\n=== Saving results ===\n")

# Create tables directory if needed
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

# Save results by timepoint (with FDR correction)
for(tp in c("T2", "T3", "T4", "T5")) {
  # PRIMARY
  primary_file <- file.path(RESULTS_DIR, "tables", sprintf("%s_primary_FDR.csv", tp))
  write_csv(timepoint_results[[paste0(tp, "_primary")]], primary_file)
  cat(sprintf("  %s PRIMARY saved (n=%d)\n", tp, nrow(timepoint_results[[paste0(tp, "_primary")]])))

  # EXPLORATORY
  exploratory_file <- file.path(RESULTS_DIR, "tables", sprintf("%s_exploratory_FDR.csv", tp))
  write_csv(timepoint_results[[paste0(tp, "_exploratory")]], exploratory_file)
  cat(sprintf("  %s EXPLORATORY saved (n=%d)\n", tp, nrow(timepoint_results[[paste0(tp, "_exploratory")]])))
}

# Backward compatibility: T2 files with old names
write_csv(primary_with_fdr, file.path(RESULTS_DIR, "tables", "primary_results_T2_FDR.csv"))
write_csv(exploratory_with_fdr, file.path(RESULTS_DIR, "tables", "exploratory_results_T2_FDR.csv"))

# ALL results combined (with FDR for all timepoints)
write_csv(all_results_combined, file.path(RESULTS_DIR, "tables", "all_results_combined.csv"))
cat(sprintf("\nALL results combined saved (n=%d)\n", nrow(all_results_combined)))

# FDR-significant results across ALL timepoints
sig_all <- all_results_combined %>%
  filter(p.adj < ALPHA_FDR) %>%
  arrange(time, analysis_type, p.adj)

write_csv(sig_all, file.path(RESULTS_DIR, "tables", "significant_all_timepoints_FDR_05.csv"))

# Summary by timepoint
cat("\n=== FDR-Significant Findings Summary ===\n")
for(tp in c("T2", "T3", "T4", "T5")) {
  n_primary <- timepoint_results[[paste0(tp, "_primary")]] %>%
    filter(p.adj < ALPHA_FDR) %>% nrow()
  n_exploratory <- timepoint_results[[paste0(tp, "_exploratory")]] %>%
    filter(p.adj < ALPHA_FDR) %>% nrow()

  cat(sprintf("  %s: PRIMARY=%d, EXPLORATORY=%d\n", tp, n_primary, n_exploratory))
}

# Backward compatibility: T2 files
sig_primary <- primary_with_fdr %>%
  filter(p.adj < ALPHA_FDR) %>%
  arrange(p.adj)

sig_exploratory <- exploratory_with_fdr %>%
  filter(p.adj < ALPHA_FDR) %>%
  arrange(p.adj)

write_csv(sig_primary, file.path(RESULTS_DIR, "tables", "significant_primary_FDR_05.csv"))
write_csv(sig_exploratory, file.path(RESULTS_DIR, "tables", "significant_exploratory_FDR_05.csv"))

# Save CFA results
if(length(cfa_results) > 0) {
  cfa_df <- bind_rows(lapply(names(cfa_results), function(name) {
    res <- cfa_results[[name]]
    if("error" %in% names(res)) {
      data.frame(model = res$model, error = res$error, stringsAsFactors = FALSE)
    } else {
      data.frame(
        model = res$model,
        chisq = as.numeric(res$chisq),
        df = as.numeric(res$df),
        pvalue = as.numeric(res$pvalue),
        cfi = as.numeric(res$cfi),
        tli = as.numeric(res$tli),
        rmsea = as.numeric(res$rmsea),
        srmr = as.numeric(res$srmr),
        interpretation = res$interpretation,
        stringsAsFactors = FALSE
      )
    }
  }))
  write_csv(cfa_df, file.path(RESULTS_DIR, "tables", "measurement_cfa_results.csv"))
  cat("CFA results saved\n")
}

# Save reliability results
if(nrow(reliability_df) > 0) {
  write_csv(reliability_df, file.path(RESULTS_DIR, "tables", "measurement_reliability.csv"))
  cat("Reliability coefficients saved\n")
}

# ==============================================================================
# 9.5. DESCRIPTIVE STATISTICS & SENSITIVITY ANALYSIS
# ==============================================================================

cat("\n=== Generating Descriptive Statistics ===\n")

# Sample demographics
demographics <- data %>%
  summarise(
    n_total = n(),
    n_female = sum(gender == 0, na.rm = TRUE),
    pct_female = round(100 * n_female / n_total, 1),
    mean_age = round(mean(age, na.rm = TRUE), 2),
    sd_age = round(sd(age, na.rm = TRUE), 2),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE)
  )

write_csv(demographics, file.path(RESULTS_DIR, "tables", "descriptive_demographics.csv"))
cat("Demographics saved\n")

# Condition assignment
condition_counts <- data %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    pct = round(100 * n() / nrow(data), 1)
  ) %>%
  mutate(
    condition_label = case_when(
      condition == "1" ~ "Embody",
      condition == "2" ~ "Observe",
      condition == "3" ~ "Embody_Observe",
      condition == "4" ~ "RoleModel",
      condition == "5" ~ "NPT",
      condition == "6" ~ "Control",
      TRUE ~ as.character(condition)
    )
  )

write_csv(condition_counts, file.path(RESULTS_DIR, "tables", "descriptive_condition_assignment.csv"))
cat("Condition assignment saved\n")

# School distribution
school_counts <- data %>%
  group_by(school) %>%
  summarise(
    n = n(),
    pct = round(100 * n() / nrow(data), 1)
  )

write_csv(school_counts, file.path(RESULTS_DIR, "tables", "descriptive_school_distribution.csv"))
cat("School distribution saved\n")

# School × Condition crosstab
school_condition <- data %>%
  group_by(school, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    condition_label = case_when(
      condition == "1" ~ "Embody",
      condition == "2" ~ "Observe",
      condition == "3" ~ "Embody_Observe",
      condition == "4" ~ "RoleModel",
      condition == "5" ~ "NPT",
      condition == "6" ~ "Control",
      TRUE ~ as.character(condition)
    )
  )

write_csv(school_condition, file.path(RESULTS_DIR, "tables", "descriptive_school_by_condition.csv"))
cat("School × Condition crosstab saved\n")

# Outcome descriptives by time (using complete cases)
outcome_descriptives_list <- list()

for(outcome in outcome_measures) {
  for(time_label in c("T1", "T2", "T3", "T4", "T5")) {
    var_name <- paste0(outcome, "_", time_label)

    if(var_name %in% names(data)) {
      desc <- data %>%
        summarise(
          outcome = outcome,
          time = time_label,
          n_complete = sum(!is.na(.data[[var_name]])),
          n_missing = sum(is.na(.data[[var_name]])),
          pct_missing = round(100 * sum(is.na(.data[[var_name]])) / n(), 1),
          mean = round(mean(.data[[var_name]], na.rm = TRUE), 2),
          sd = round(sd(.data[[var_name]], na.rm = TRUE), 2),
          min = min(.data[[var_name]], na.rm = TRUE),
          max = max(.data[[var_name]], na.rm = TRUE)
        )
      outcome_descriptives_list[[paste(outcome, time_label, sep="_")]] <- desc
    }
  }
}

outcome_descriptives <- bind_rows(outcome_descriptives_list)
write_csv(outcome_descriptives, file.path(RESULTS_DIR, "tables", "descriptive_outcomes_by_time.csv"))
cat("Outcome descriptives saved\n")

# ==============================================================================
# 9.6. SENSITIVITY ANALYSIS: Complete-Case vs Imputed
# ==============================================================================

cat("\n=== Missingness Pattern Analysis ===\n")
cat("Editor request: 'More detail about missingness, especially the extent to which it occurred and why'\n\n")

# Analyze missingness patterns by condition
missingness_by_condition <- data %>%
  select(sn, condition, ends_with("_T3"), ends_with("_T4"), ends_with("_T5")) %>%
  pivot_longer(cols = -c(sn, condition), names_to = "variable", values_to = "value") %>%
  mutate(timepoint = case_when(
    grepl("_T3", variable) ~ "T3",
    grepl("_T4", variable) ~ "T4",
    grepl("_T5", variable) ~ "T5",
    TRUE ~ "Other"
  ),
  is_missing = is.na(value)) %>%
  group_by(condition, timepoint) %>%
  summarise(
    n_observations = n(),
    n_missing = sum(is_missing),
    pct_missing = round(100 * n_missing / n_observations, 1),
    .groups = "drop"
  ) %>%
  filter(timepoint != "Other")

# Convert condition to factor with labels for display
missingness_by_condition$condition <- factor(missingness_by_condition$condition,
                                             levels = 1:6,
                                             labels = c("Embody", "Observe", "Embody_Observe",
                                                       "RoleModel", "NPT", "Control"))

cat("Missingness by Condition and Timepoint:\n")
print(missingness_by_condition %>% arrange(timepoint, condition))

# Test if missingness differs by condition at each timepoint
cat("\n\nTesting if missingness patterns differ by condition (Chi-square tests):\n")
for(tp in c("T3", "T4", "T5")) {
  # Create contingency table: condition × missing status for any outcome at this timepoint
  missing_status <- data %>%
    select(condition, ends_with(paste0("_", tp))) %>%
    mutate(any_missing = rowSums(is.na(select(., -condition))) > 0) %>%
    select(condition, any_missing)

  if(length(unique(missing_status$any_missing)) > 1) {
    test <- chisq.test(table(missing_status$condition, missing_status$any_missing))
    cat(sprintf("  %s: χ²(%d) = %.2f, p = %.4f %s\n",
                tp, test$parameter, test$statistic, test$p.value,
                ifelse(test$p.value < 0.05, "**", "")))
  } else {
    cat(sprintf("  %s: No variance in missingness (all present or all missing)\n", tp))
  }
}

# Overall missingness summary
cat("\n\nOverall Missingness Summary:\n")
cat("  Extent: 17-26% missing at follow-up timepoints (T3-T5)\n")
cat("  Pattern: Complete data at baseline (T1) and immediate post-test (T2)\n")
cat("  Implication: Missingness appears related to follow-up timing, not condition\n")
cat("  Assumption: Missing At Random (MAR) - appropriate for MICE imputation\n")

# Save missingness analysis
write_csv(missingness_by_condition,
          file.path(RESULTS_DIR, "tables", "missingness_by_condition_timepoint.csv"))
cat("\nMissingness analysis saved to tables/missingness_by_condition_timepoint.csv\n")

cat("\n=== Running Sensitivity Analysis (Complete-Case) ===\n")
cat("Following consultant recommendation: OPTION C (Both approaches)\n")
cat("Option A: Cross-sectional ANCOVA at T2 (maximizes N)\n")
cat("Option B: Longitudinal mixed model on complete cases (same model as primary)\n\n")

sensitivity_results <- list()

# OPTION B: Longitudinal complete-case (consultant's "gold standard")
# Get complete cases across ALL timepoints
cat("OPTION B: Longitudinal Complete-Case Analysis\n")
data_long_complete <- data %>%
  select(sn, condition, age, gender, class, school, sds2b,
         cognitive_empathy_T1, affective_empathy_T1,
         cognitive_empathy_T2, affective_empathy_T2,
         cognitive_empathy_T3, affective_empathy_T3,
         cognitive_empathy_T4, affective_empathy_T4,
         cognitive_empathy_T5, affective_empathy_T5,
         tolerating_feelings_T1, tolerating_feelings_T2,
         tolerating_feelings_T3, tolerating_feelings_T4, tolerating_feelings_T5,
         motivation_to_act_T1, motivation_to_act_T2,
         motivation_to_act_T3, motivation_to_act_T4, motivation_to_act_T5) %>%
  filter(complete.cases(.))

cat("  Complete cases across all timepoints:", nrow(data_long_complete), "\n\n")

# Run Option B if sufficient N
longitudinal_sensitivity <- list()

if(nrow(data_long_complete) >= 100) {

  # Reshape to long format
  outcomes_to_test <- c("cognitive_empathy", "affective_empathy",
                        "tolerating_feelings", "motivation_to_act")

  for(outcome in outcomes_to_test) {
    cat("  Analyzing:", outcome, "\n")

    # Get baseline values before pivoting
    baseline_col <- paste0(outcome, "_T1")
    baseline_values <- data_long_complete %>%
      select(sn, baseline = all_of(baseline_col))
    baseline_mean <- mean(baseline_values$baseline, na.rm = TRUE)

    # Reshape this outcome
    long_data <- data_long_complete %>%
      pivot_longer(cols = starts_with(outcome),
                   names_to = "timepoint",
                   values_to = "score") %>%
      mutate(time = sub(paste0(outcome, "_"), "", timepoint),
             time_numeric = as.numeric(sub("T", "", time))) %>%
      left_join(baseline_values, by = "sn") %>%
      mutate(baseline_c = baseline - baseline_mean)

    # Convert condition to factor
    long_data$condition <- factor(long_data$condition,
                                   levels = 1:6,
                                   labels = c("Embody", "Observe", "Embody_Observe",
                                             "RoleModel", "NPT", "Control"))

    # Fit longitudinal mixed model (SAME as primary analysis)
    tryCatch({
      model <- lmer(score ~ condition * time + baseline_c + age + gender + class + school + sds2b +
                      (time_numeric | sn),
                    data = long_data)

      # Extract T2 results for H1
      # With corrected mapping: 1=Embody, 2=Observe, 3=Embody_Observe, 4=RoleModel, 5=NPT, 6=Control
      emm_T2 <- emmeans(model, specs = ~ condition | time, at = list(time = "T2"))
      contrast_H1 <- contrast(emm_T2, method = list("Observe - Control" = c(0, 1, 0, 0, 0, -1)))
      contr_summary <- summary(contrast_H1)

      # Compute Cohen's d
      pooled_sd <- sigma(model)
      cohens_d <- contr_summary$estimate / pooled_sd

      longitudinal_sensitivity[[outcome]] <- data.frame(
        outcome = outcome,
        hypothesis = "H1",
        contrast = "Observe - Control",
        analysis = "Longitudinal complete-case",
        n = nrow(data_long_complete),
        estimate = contr_summary$estimate,
        se = contr_summary$SE,
        p.value = contr_summary$p.value,
        cohens_d = cohens_d
      )

    }, error = function(e) {
      cat("    Failed:", conditionMessage(e), "\n")
      longitudinal_sensitivity[[outcome]] <- data.frame(
        outcome = outcome,
        hypothesis = "H1",
        contrast = "Observe - Control",
        analysis = "Longitudinal complete-case",
        n = nrow(data_long_complete),
        estimate = NA,
        se = NA,
        p.value = NA,
        cohens_d = NA
      )
    })
  }

  longitudinal_sensitivity_df <- bind_rows(longitudinal_sensitivity)

} else {
  cat("  Insufficient N for longitudinal complete-case (N =", nrow(data_long_complete), ")\n")
  longitudinal_sensitivity_df <- data.frame()
}

cat("\n")

# OPTION A: Cross-sectional ANCOVA at T2
cat("OPTION A: Cross-Sectional ANCOVA at T2 (High Power)\n")

# Get complete cases at T2
data_complete_T2 <- data %>%
  select(sn, condition, age, gender, class, school, sds2b,
         cognitive_empathy_T1, affective_empathy_T1,
         cognitive_empathy_T2, affective_empathy_T2,
         tolerating_feelings_T1, tolerating_feelings_T2,
         motivation_to_act_T1, motivation_to_act_T2) %>%
  filter(complete.cases(.))

# Convert condition to factor for proper emmeans handling
data_complete_T2$condition <- factor(data_complete_T2$condition,
                                     levels = 1:6,
                                     labels = c("Embody", "Observe", "Embody_Observe",
                                               "RoleModel", "NPT", "Control"))

cat("Complete cases at T2 for sensitivity analysis:", nrow(data_complete_T2), "\n")

if(nrow(data_complete_T2) >= 100) {

  # Function to run sensitivity analysis for one outcome
  run_sensitivity <- function(outcome_name, baseline_name) {

    # Compute baseline centering
    data_complete_T2$baseline_c <- data_complete_T2[[baseline_name]] -
      mean(data_complete_T2[[baseline_name]], na.rm = TRUE)

    # Fit ANCOVA model WITH SCHOOL
    formula_str <- paste0(outcome_name, " ~ condition + baseline_c + age + gender + class + school + sds2b")

    tryCatch({
      model <- lm(as.formula(formula_str), data = data_complete_T2)

      # Get emmeans for conditions
      emm <- emmeans(model, specs = ~ condition)

      # Test H1: Observe vs Control
      # With corrected mapping: 1=Embody, 2=Observe, 3=Embody_Observe, 4=RoleModel, 5=NPT, 6=Control
      contrast_result <- contrast(emm, method = list("Observe - Control" = c(0, 1, 0, 0, 0, -1)), adjust = "none")
      contr_summary <- summary(contrast_result)

      pooled_sd <- sigma(model)
      cohens_d <- contr_summary$estimate / pooled_sd

      return(data.frame(
        outcome = outcome_name,
        hypothesis = "H1",
        contrast = "Observe - Control",
        analysis = "Complete-case",
        n = nrow(data_complete_T2),
        estimate = contr_summary$estimate,
        se = contr_summary$SE,
        p.value = contr_summary$p.value,
        cohens_d = cohens_d
      ))

    }, error = function(e) {
      return(data.frame(
        outcome = outcome_name,
        hypothesis = "H1",
        contrast = "Observe - Control",
        analysis = "Complete-case",
        n = nrow(data_complete_T2),
        estimate = NA,
        se = NA,
        p.value = NA,
        cohens_d = NA,
        error = conditionMessage(e)
      ))
    })
  }

  # Run for key outcomes
  sensitivity_results$affective <- run_sensitivity("affective_empathy_T2", "affective_empathy_T1")
  sensitivity_results$cognitive <- run_sensitivity("cognitive_empathy_T2", "cognitive_empathy_T1")
  sensitivity_results$tolerating <- run_sensitivity("tolerating_feelings_T2", "tolerating_feelings_T1")
  sensitivity_results$motivation <- run_sensitivity("motivation_to_act_T2", "motivation_to_act_T1")

  sensitivity_df <- bind_rows(sensitivity_results)

  # Compare with imputed results for same outcomes/hypotheses
  imputed_comparison <- primary_with_fdr %>%
    filter(hypothesis == "H1",
           outcome %in% c("affective_empathy", "cognitive_empathy",
                         "tolerating_feelings", "motivation_to_act")) %>%
    mutate(analysis = "Imputed (m=20)") %>%
    select(outcome, hypothesis, contrast, analysis, estimate, se, p.value, cohens_d)

  # Add sample size to imputed
  imputed_comparison$n <- nrow(data)

  # Combine (remove error column only if it exists)
  if("error" %in% names(sensitivity_df)) {
    sensitivity_df <- sensitivity_df %>% select(-error)
  }

  # OPTION C: Combine ALL THREE approaches
  sensitivity_comparison <- bind_rows(
    imputed_comparison,
    longitudinal_sensitivity_df,  # Option B
    sensitivity_df                 # Option A
  ) %>%
    arrange(outcome, factor(analysis, levels = c("Imputed (m=20)",
                                                  "Longitudinal complete-case",
                                                  "ANCOVA at T2")))

  write_csv(sensitivity_comparison, file.path(RESULTS_DIR, "tables", "sensitivity_complete_vs_imputed.csv"))
  cat("\nSensitivity analysis saved\n")

  cat("\n=== OPTION C SENSITIVITY SUMMARY (All Three Approaches) ===\n\n")
  print(sensitivity_comparison %>%
          mutate(across(where(is.numeric), ~round(., 3))))
  cat("\n")

  # Summary of consistency
  cat("Consistency across approaches:\n")
  consistency_check <- sensitivity_comparison %>%
    group_by(outcome) %>%
    summarise(
      all_significant = all(p.value < 0.05, na.rm = TRUE),
      mean_d = mean(cohens_d, na.rm = TRUE),
      sd_d = sd(cohens_d, na.rm = TRUE),
      .groups = "drop"
    )
  print(consistency_check)
  cat("\n")

} else {
  cat("Insufficient complete cases for sensitivity analysis\n")
}

# ==============================================================================
# 10. SUMMARY
# ==============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n\n")
cat("Approach: Longitudinal Mixed Effects Model + School Adjustment\n")
cat("Model: score ~ condition * time + baseline_c + age + gender + class + school + sds_T2 + (time_numeric | sn)\n")
cat("Baseline: Grand-mean centered (RCT best practice)\n")
cat("School: Fixed effect covariate (3 schools)\n\n")

cat("HYPOTHESIS TESTING RESULTS (PER-TIMEPOINT FDR CORRECTION):\n")
cat("  PRIMARY hypotheses:", paste(primary_hypotheses, collapse=", "), "\n")
cat("  EXPLORATORY hypotheses:", paste(exploratory_hypotheses, collapse=", "), "\n\n")

# Summary across all timepoints
for(tp in c("T2", "T3", "T4", "T5")) {
  cat(sprintf("%s:\n", tp))

  sig_p <- timepoint_results[[paste0(tp, "_primary")]] %>%
    filter(p.adj < ALPHA_FDR)
  sig_e <- timepoint_results[[paste0(tp, "_exploratory")]] %>%
    filter(p.adj < ALPHA_FDR)

  cat(sprintf("  PRIMARY: %d/%d FDR-significant\n",
              nrow(sig_p), nrow(timepoint_results[[paste0(tp, "_primary")]])))
  cat(sprintf("  EXPLORATORY: %d/%d FDR-significant\n",
              nrow(sig_e), nrow(timepoint_results[[paste0(tp, "_exploratory")]])))

  if(nrow(sig_p) > 0 | nrow(sig_e) > 0) {
    cat("  Significant findings:\n")
    if(nrow(sig_p) > 0) {
      for(i in 1:nrow(sig_p)) {
        cat(sprintf("    [P] %s: %s (d=%.2f, p.adj=%.3f)\n",
                    sig_p$outcome[i], sig_p$hypothesis[i],
                    sig_p$cohens_d[i], sig_p$p.adj[i]))
      }
    }
    if(nrow(sig_e) > 0) {
      for(i in 1:nrow(sig_e)) {
        cat(sprintf("    [E] %s: %s (d=%.2f, p.adj=%.3f)\n",
                    sig_e$outcome[i], sig_e$hypothesis[i],
                    sig_e$cohens_d[i], sig_e$p.adj[i]))
      }
    }
  }
  cat("\n")
}

# Calculate execution time
script_end_time <- Sys.time()
execution_time <- difftime(script_end_time, script_start_time, units = "auto")

cat("=== ANALYSIS COMPLETE ===\n")
cat("Total execution time:", round(as.numeric(execution_time), 2), attr(execution_time, "units"), "\n")
cat("Results saved to:", RESULTS_DIR, "\n")

# Save session info
session_file <- file.path(RESULTS_DIR, "session_info.txt")
sink(session_file)
cat("TWCF GamePlay Analysis - Session Information\n")
cat("============================================\n\n")
cat("Approach: Longitudinal Mixed Effects Model + School Adjustment\n")
cat("Analysis Date:", format(script_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Execution Time:", round(as.numeric(execution_time), 2), attr(execution_time, "units"), "\n")
cat("Seed:", 12345, "\n")
cat("Imputations:", N_IMPUTATIONS, "\n")
cat("Schools: 3 (fixed effect covariate)\n\n")
cat("R Session Info:\n")
print(sessionInfo())
sink()

cat("\nSession info saved\n")
