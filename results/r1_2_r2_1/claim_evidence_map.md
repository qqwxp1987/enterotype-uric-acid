# R1.2 / R2.1 claim-evidence map

| Intended claim | Evidence source | Statistical support | Allowed wording |
|---|---|---|---|
| Most participants have confident assignments conditional on the prespecified k=4 fit. | `posterior_assignment_details.csv`; Figure S4A | Mean maximum posterior 0.988; 7/305 below 0.70; minimum top-two margin 0.036 | “high assignment confidence for most participants, conditional on k=4” |
| k=4 components are often recovered under bootstrap perturbation, but recovery differs by component. | `bootstrap_cluster_stability_summary.csv`; Figure S4B | Unique-participant agreement mean 0.925, SD 0.054; ARI median 0.885, 2.5th–97.5th percentiles 0.681–0.963; median Jaccard Type 1/2/3/4 = 0.957/0.920/0.794/0.877 | “conditional recovery was strongest for Types 1 and 2 and less uniform for the Type 3/4 boundary” |
| The number of components is not uniquely resolved. | `k_sensitivity_summary.csv`; Figure S4C | Winners across 100 starts: k=3 79, k=4 18, k=5 3; best Laplace k=3 148829.5 vs k=4 148845.1 | “k=3 and k=4 are competing resolutions; k=3 is modestly favoured by the multi-start fit criterion” |
| The Type-1-like HUA signal is retained under k=3, whereas the merged Prevotella branch does not separate from the reference. | `k4_k3_membership_crosswalk.csv`; `k3_outcome_sensitivity.csv` | Sex/age HUA OR C1 vs C2 2.48 (95% CI 1.17–5.25), P=0.018; C3 vs C2 1.33 (0.55–3.20), P=0.524. Renal/cardiometabolic ORs 3.38, P=0.003 and 1.43, P=0.466. | “Type 1 is resolution-robust; the Type 4-specific clinical signal is resolution-sensitive and exploratory” |
| The full-sample P/B ratio has a modest negative association with UA, but no reliable nonlinear evidence. | `pb_ratio_sensitivity.csv`; Figure S4D | Partial Spearman rho=-0.152, P=0.008; natural-spline df=3 P=0.865 and df=4 P=0.348; linear P/B and enterotype AIC 3516.1 vs 3515.9 | “modest full-sample association; no evidence that curvature improves on a linear gradient” |
| The P/B result is sensitive to the zero-containing population structure. | `pb_ratio_sensitivity.csv` | Prevotella absent in 122/305; detected-only rho=-0.102, P=0.170; extended HUA OR=0.991, P=0.654 | “the full-sample association partly reflects Prevotella absence/presence and does not establish an HUA association” |
| There is no direct external replication of the four-state structure plus Type 1/Type 4 HUA pattern. | Manuscript Limitations; Chu 2021 primary paper | Chu 2021 used PAM and gout vs controls; enterotype distribution P=0.8541 | “contextual taxonomic/functional concordance only; not external validation” |

## Forbidden formulations after this repair

- “four stable enterotypes” or “the four-type structure is stable” without the
  qualifier “conditional on k=4”;
- “no participant was ambiguous”;
- “significant nonlinear P/B-UA relationship”;
- “equivalent explanatory power” from Delta AIC alone;
- “Chu 2021 externally validates our core conclusions”;
- “the negative Chu enterotype result does not contradict our findings”;
- any use of UA/HUA outcomes to justify selection of k=4 over k=3.
