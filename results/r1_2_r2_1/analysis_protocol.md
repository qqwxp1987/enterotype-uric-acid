# R1.2 / R2.1 enterotype robustness analysis protocol

Date frozen: 2026-08-28

## Scope

This is a focused reviewer-repair analysis. It does not redefine the primary
enterotypes or use UA/HUA to choose the number of clusters. It separates:

1. assignment confidence conditional on the prespecified `k = 4` solution;
2. recovery of each `k = 4` cluster under bootstrap perturbation;
3. uncertainty in the selected number of clusters across DMM initializations;
4. sensitivity of UA/HUA associations to a `k = 3` representation; and
5. sensitivity of the continuous Prevotella/Bacteroides analysis to zeros and
   pseudocount handling.

## Frozen hypotheses and interpretation rules

### H1. Conditional k=4 recovery

- Prediction: the four prespecified components will show non-random recovery
  under 100 fixed-seed bootstrap refits.
- Primary metrics: per-component Jaccard similarity and adjusted Rand index
  (ARI) on unique in-bag participants.
- Interpretation: median Jaccard >= 0.75 is treated as stable conditional
  recovery; 0.60-0.75 as pattern recovery with membership uncertainty; < 0.60
  as weak recovery. These metrics cannot establish that `k = 4` is the unique
  optimal cluster count.
- All 100 iterations, failed fits, and all four component results will be
  retained and reported.

### H2. Cluster-count uncertainty

- Prediction: the existing 100-start scan will reproduce the known competition
  between `k = 3` and `k = 4`.
- Metrics: winner frequencies, per-start `Laplace(k=4)-Laplace(k=3)`, and the
  best (minimum) Laplace value attained for each of `k = 3, 4, 5`.
- Interpretation: the lowest attained Laplace value is the primary multi-start
  comparison. A fixed-seed `k = 4` result will not be described as proof that
  four is the uniquely supported count.

### H3. Outcome sensitivity to k=3

- Prediction: the Type-1-like signal will remain detectable, whereas the
  Type-4-specific HUA pattern may be attenuated when the Prevotella branch is
  represented as one merged `k = 3` component.
- Models: the same baseline sex/age GLS and logistic models used for the main
  analysis, plus the prespecified renal/cardiometabolic adjustment model.
- Interpretation: HUA/UA results are sensitivity outcomes only; they will not
  be used to select `k` or to claim that one clustering solution is superior.

### H4. Continuous-gradient sensitivity

- Prediction: the modest P/B-UA association may persist, but the nonlinear
  P-value may vary with zero replacement because Prevotella has structural or
  sampling zeros in many participants.
- Specifications: pseudocounts `1e-6`, `1e-5`, `1e-4`, and `1e-3`; half of the
  minimum observed positive Prevotella abundance; and a Prevotella-detected
  subset.
- Interpretation: all specifications will be reported. A significant smooth
  term will be described as evidence for a nonlinear continuous association,
  not as proof of discrete states. Absence of a binary-HUA association will be
  reported as a null result, not as evidence of equivalence or threshold
  specificity.

## Baselines to reproduce before new inference

- mean maximum DMM posterior probability: approximately 0.988;
- 100-start winners: `k=3: 79`, `k=4: 18`, `k=5: 3`;
- legacy row-weighted bootstrap agreement: approximately 0.931 +/- 0.049;
- P/B-UA partial Spearman: approximately -0.152, P approximately 0.008;
- main P/B nonlinear LRT: approximately P=0.006;
- linear AIC comparison: enterotype approximately 3515.9 and P/B approximately
  3516.1.

If these baselines do not reproduce within rounding tolerance, the analysis
stops for debugging before manuscript changes.

## Protocol amendment after baseline reproduction

The legacy GAM comparison was reproduced before replacement. Inspection then
showed that the penalized smooth had effective degrees of freedom approximately
equal to 1, while `anova.gam(linear, smooth)` had been interpreted as a test of
curvature. That interpretation is not defensible: an edf near 1 indicates an
approximately linear fitted relationship. The primary nonlinearity test was
therefore changed, before the full bootstrap run, to a nested comparison of the
linear model against fixed-degree natural-spline models (`df = 3`, with `df =
4` as sensitivity). The legacy GAM estimate and P-value remain in the audit
table, explicitly labelled as retired diagnostics, so the unfavorable change
is fully traceable rather than deleted.

The first smoke run also revealed that a shared complete-case filter would have
reduced the continuous-UA analysis to the 225 participants with BMI data. The
canonical implementation now uses model-specific complete cases: 305 for the
UA/gradient models and 225 for the extended HUA model, matching the original
estimands and baseline outputs.

The completed 100-bootstrap run reproduced 96/100 legacy row-agreement values
exactly. Four local-optimum fits differed, yielding a current mean of 0.9273
versus the legacy 0.9310 (0.39% relative difference). This is within the frozen
<=1% reproduction tolerance. The iteration-level comparison is retained in
`legacy_bootstrap_baseline_comparison.csv`; the current complete run, rather
than the legacy cache, is the evidence source for revised text.

## Reproducibility

- DMM bootstrap resampling seeds: 1001-1100.
- DMM fit seeds: 1-100.
- k=3 sensitivity fit: seed 27, the minimum-Laplace k=3 fit in the frozen
  100-start scan.
- Raw data are read only.
- The canonical script will write complete iteration-level tables, summary
  tables, figure files, an execution log, and `sessionInfo()`.
