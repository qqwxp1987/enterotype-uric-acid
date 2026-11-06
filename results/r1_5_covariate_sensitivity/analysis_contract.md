# R1.5 Extended-Covariate Sensitivity Analysis Contract

Date frozen: 2026-09-01

## Scope

This supplement addresses Reviewer #1 Major Comment 5 without replacing the
original sex-and-age-adjusted descriptive analyses.

1. Differential-abundance sensitivity for the original four contrasts and all
   four data layers (species, pathways, ECs, and serum metabolites).
2. Cross-omics sensitivity for the exact pre-existing Figure 5 feature sets:
   the main species panel, the functional panel, and the global partial db-RDA.
3. Direct feature-HUA logistic sensitivity for the same 22 pre-specified
   candidates used in the existing feature-attribution workflow.
4. Complete-case selection diagnostics for BMI missingness, because BMI is the
   only extended covariate with material missingness and can change the analysis
   population independently of covariate adjustment.

The manuscript and response letter are out of scope until the user reviews the
results.

## Locked analysis sets

- Baseline adjustment: sex + age.
- Extended adjustment: sex + age + BMI + creatinine (Cr).
- Cr and eGFR are not entered together. No additional model is selected after
  seeing the results.
- Differential-abundance contrasts: E1 vs E2, E1 vs E3, E4 vs E2, E4 vs E3.
- LinDA and limma settings, transformations, prevalence filters, multiplicity
  families, and feature-selection rules match the original scripts.
- Cross-omics feature sets are copied verbatim from the existing Figure 5
  scripts; no feature is added or removed based on the new results.
- BH correction is applied within the same analysis families as the original
  workflows.
- Matched complete-case base and extended models use identical samples so the
  covariate effect is not mixed with a sample-composition change.

## Baseline gate

Before extended models are interpreted, the script must freshly reproduce:

- the four stored differential-analysis tables;
- the stored Figure 5 species and functional partial-correlation matrices;
- the stored partial db-RDA result; and
- the stored sex+age direct feature-HUA logistic results.

If tested features, comparison keys, core sets, or numerical results do not
match within the scripted tolerance, the run stops before extended-model
interpretation.

## Pre-run hypotheses and decision rules

### H1: Enterotype differential signatures

Prediction: most original core-feature contrast estimates will retain their
direction after BMI and Cr are added, although FDR retention may fall because
of complete-case restriction and added covariates.

- Directionally supportive: at least 80% of original core features retain the
  original direction in both defining contrasts for each principal layer.
- Directionally refuting: fewer than 70% retain direction in both contrasts, or
  two or more named central features reverse direction.
- Otherwise: ambiguous.
- FDR retention is reported descriptively and is not used alone to declare
  biological stability.

### H2: Cross-omics structure

Prediction: matched-base and extended partial correlations will remain strongly
concordant, and the metabolome-constrained db-RDA signal will remain detectable.

- Supportive: Spearman correlation between matched-base and extended rho values
  is at least 0.90, direction concordance is at least 90%, and the extended
  db-RDA permutation P value is below 0.05.
- Refuting: rho concordance is below 0.70 or the extended db-RDA signal is not
  detectable (P at least 0.05).
- Otherwise: ambiguous.

### H3: Direct HUA associations

Prediction: the directions for PWY0-1296, EC 2.6.1.1, and
5-methoxyindoleacetate will remain concordant between matched-base and extended
logistic models.

- Supportive: all three retain direction and at least one of the 22 candidates
  remains BH-significant.
- Refuting: two or more of the three reverse direction.
- Otherwise: ambiguous.

## Reporting and stop rules

- All tested features, all contrasts, all model failures, and all unfavorable
  results are retained in the output directory.
- Directional stability, FDR stability, and precision are reported separately.
- No causal or mechanistic interpretation is permitted from these cross-sectional
  sensitivity analyses.
- BMI-complete versus BMI-missing participants are compared on enterotype, HUA,
  UA, sex, age, and Cr to distinguish covariate effects from selection effects.
- A complete execution log, session information, input/output checksums, and
  machine-readable summary tables are required before any result claim.
