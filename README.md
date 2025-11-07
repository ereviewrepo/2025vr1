# VR Empathy and Compassion Training: Data and Analysis Code

This repository contains the data and analysis code for the manuscript:

**"Virtual Reality-Based Empathy and Compassion Training for Nigerian Adolescents: A Randomized Controlled Trial"**

Submitted to: *Technology, Mind, and Behavior*

## Repository Contents

### Data File
- **`2025vr1data.csv`** - Anonymized dataset (N=377 participants after attention check exclusions)
  - Participant names removed
  - Schools anonymized (sch1, sch2, sch3)
  - 251 columns including demographics, outcome measures across 5 timepoints

### Analysis Script
- **`2025vr1Rscript.R`** - Complete statistical analysis script
  - Longitudinal mixed effects models
  - Multiple imputation for missing data (MICE, m=20 imputations)
  - FDR correction for multiple comparisons
  - Measurement validation (CFA and reliability)
  - School included as fixed effect covariate

## Study Design

- **Design**: 6-condition randomized controlled trial × 5 timepoints
- **Sample**: N=377 Nigerian adolescents (ages 14-18)
- **Conditions**:
  1. Control (standard curriculum)
  2. Observe (VR perspective-taking)
  3. Embody (VR embodiment)
  4. RoleModel (VR role modeling)
  5. Embody+Observe (combined conditions)
  6. NPT (neutral procedural training, active control)
- **Timepoints**: T1 (baseline), T2 (immediate post-test), T3 (3 weeks), T4 (6 weeks), T5 (3 months)
- **Schools**: 3 Nigerian secondary schools (individual-level randomization within schools)

## Primary Outcomes

1. Affective empathy (AMES)
2. Cognitive empathy (IRI Perspective-Taking)
3. Composite empathy (standardized combination)
4. Tolerating distressing feelings (Swiss Compassion Scale)
5. Motivation to act compassionately (Swiss Compassion Scale)

## Statistical Approach

### Primary Analysis
- **Model**: Longitudinal mixed effects with condition × time interaction
- **Formula**: `outcome ~ condition * time + baseline_centered + age + gender + class + school + (time | participant)`
- **Missing Data**: Multiple imputation (MICE, m=20, seed=12345)
- **Multiple Comparisons**: Benjamini-Hochberg FDR correction (α=.05)
  - Applied separately within each timepoint and outcome
  - T2 = primary endpoint (N=377, complete data)
  - T3-T5 = secondary endpoints (17-26% attrition)
  - PRIMARY hypotheses (5): H1, H2, H3, H7, H8
  - EXPLORATORY hypotheses (5): H4, H5a, H5b, H5c, H6

### Sensitivity Analysis
Three approaches implemented to assess robustness:
1. Multiple imputation (primary analysis, N=377)
2. Longitudinal complete-case analysis (same model, N=201)
3. Cross-sectional ANCOVA at T2 (N=371)

### School as Covariate
- School included as fixed effect covariate (k=3 schools)
- Controls for school-level differences in resources and culture
- Individual-level randomization within schools prevents systematic confounding

## Running the Analysis

### Requirements
R packages required:
- tidyverse
- mice
- emmeans
- lme4
- lmerTest
- lavaan
- psych

### Instructions

1. Ensure `2025vr1data.csv` is in the same directory as the R script
2. Run: `Rscript 2025vr1Rscript.R`
3. Expected runtime: ~19 minutes
4. Output: `results_preregistered_[timestamp]/` folder with:
   - 15 CSV tables with results
   - Diagnostic plots
   - Session info for reproducibility

### Reproducibility
- Fixed seed: 12345 throughout
- Temporal ordering enforced in imputation
- Complete session info saved with results

## Key Findings

### T2 (Immediate Post-Intervention)
**6 PRIMARY hypotheses reached FDR-corrected significance (p.adj < .05):**
- Affective empathy: Observe vs Control (d=0.51), Observe vs NPT (d=0.62)
- Composite empathy: Observe vs Control (d=0.55), Observe vs NPT (d=0.62)
- Tolerating feelings: Observe vs Control (d=0.53), Observe vs NPT (d=0.51)

**2 EXPLORATORY hypotheses significant:**
- Affective empathy: Observe vs Embody (d=0.51), Observe vs Embody+Observe (d=0.52)

All significant findings involve the Observe (third-person VR observation) condition.

### T3 (2-Month Follow-Up)
**No findings reached FDR-corrected significance.** Two comparisons showed raw p<.05 with medium-to-large effect sizes but did not survive FDR correction:
- H2: Embody vs Control on IOS closeness (d=0.61, p=.026, p_FDR=.129)
- H5a: Observe+Embody vs Control on compassionate feeling (d=0.63, p=.024, p_FDR=.119)

These findings suggest a pattern shift from third-person observation (T2) to first-person embodiment modes (T3).

### T4-T5 (4-6 Month Follow-Ups)
No effects reached significance - complete attenuation by 4 months.

**Overall**: Third-person VR observation shows robust immediate effects that do not persist beyond 2 months. Single-session intervention appears insufficient for lasting change.

See manuscript for complete results and interpretation.

## Data Dictionary

### Participant Identifiers
- `sn`: Participant serial number (unique ID)
- `idcode`: Study ID code
- `condition`: Experimental condition (1=Embody, 2=Observe, 3=Embody_Observe, 4=RoleModel, 5=NPT, 6=Control)
- `school`: School identifier (sch1, sch2, sch3 - anonymized)

### Demographics
- `gender`: 0=Male, 1=Female
- `age`: Age in years
- `class`: Grade level (1=JSS3, 2=SS1, etc.)
- `height`: Height in cm
- `weight`: Weight in kg
- Additional demographic variables (father_origin, religion, sibling counts, parental occupation/education)

### Outcome Measures

All measures repeated across timepoints (suffix: a=T1 baseline, b=T2 immediate, c=T3 3wk, d=T4 6wk, e=T5 3mo)

**Sussex Oxford Compassion Scales (SOCS)**: `sussex1a`-`sussex20e` (20 items × 5 timepoints)
- Tolerating distressing feelings
- Recognizing suffering
- Understanding universality
- Feeling for person
- Motivation to act

**Affective and Cognitive Measure of Empathy (AMES)**: `ames1a`-`ames8e` (8 items × 5 timepoints)
- Affective empathy (vicarious emotion)
- Cognitive empathy (perspective-taking)

**Self-assessment of sympathetic responding**: `symp1a`-`symp3e` (3 items × 5 timepoints)

**Interpersonal Reactivity Index (IRI)**: `iri1a`-`iri5e` (5 items × 5 timepoints)
- Perspective-taking subscale

**Social Desirability Scale (SDS)**: `sds1a`-`sds6b` (6 items, T1 and T2 only)

**Swiss Compassion Scale**: `swiss1`-`swiss8` (8 items, T2 only)
- Tolerating distressing feelings
- Motivation to act compassionately

**Health measures**: `healtha`-`healthe` (5 items × 5 timepoints)

**Attention checks**: `attnchk1`, `attnchk2`, `attnchk3` (quality control items)

**Additional measures**: Inclusion of Other in Self (ios), Identification with All Humanity (sps), Embodiment (emb)

## Citation

If you use this data or code, please cite:

[Citation to be added upon publication]

## License

**IMPORTANT**: This data and code are provided for **transparency and replication purposes ONLY**.

See [LICENSE](LICENSE) file for complete terms. Key restrictions:
- ✅ Permitted: Verification and exact replication of published analyses
- ❌ Prohibited: Use in new research, redistribution, commercial use, or derivative works

For any use beyond exact replication, you must obtain written permission from the corresponding author.

## Contact

For questions about the data or analysis, please contact the corresponding author via the journal.

---

**Data Version**: 2025vr1
**Analysis Date**: October 2025
 
