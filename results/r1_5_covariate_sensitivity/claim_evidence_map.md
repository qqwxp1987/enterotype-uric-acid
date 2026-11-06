# R1.5 Claim–Evidence Map

| Intended statement | Evidence source | Formal paper location | Statistical support / boundary |
|---|---|---|---|
| Covariate sets were selected according to the estimand, with BMI and creatinine added in matched complete-case sensitivity analyses. | `analysis_contract.md`; `r15_extended_covariate_sensitivity.R` | Methods 2.10 | Primary descriptive models remain sex+age; additional adjustment is sensitivity, not a replacement estimand. |
| Original differential-feature directions were mostly retained, but FDR-level core retention varied. | `differential_core_stability_summary.csv` | Results 3.4; Table S10 | Direction retention 82.7%–100% across non-empty layer-by-enterotype sets; FDR retention ranges 0%–95.8%. |
| Loss of Type 4 taxonomic core status cannot be attributed solely to BMI/Cr adjustment. | `differential_core_stability_summary.csv`; `bmi_missingness_by_enterotype.csv` | Results 3.4; Limitations; Table S10 | The three Type 4 taxa were already absent from the matched sex+age core set; Type 4 BMI completeness was 22/39. |
| Cross-omics structure remained highly concordant after BMI/Cr adjustment. | `crossomics_correlation_summary.csv`; `crossomics_dbrda_summary.csv` | Results 3.5; Table S10 | Species/function coefficient concordance 0.999/0.998; direction concordance 98.9%/99.4%; adjusted db-RDA R2=0.159, P=0.010. |
| Direct candidate-HUA associations were directionally stable, with EC 2.6.1.1 retaining FDR significance. | `direct_feature_HUA_sensitivity.csv` | Results 3.5; Table S8 | 21/22 directions retained; EC 2.6.1.1 OR=14.82, 95% CI 3.23–67.96, FDR=0.0115. OR is scale-dependent. |
| Mediation sensitivity supports directional rather than complete significance-level robustness. | `mediation_extended_covariate_sensitivity.csv` | Results 3.6; Table S9 | 42/44 directions retained; original 10/10 FDR-positive directions retained but 2/10 remained FDR-positive in the smaller adjusted models. |
| Complete-case attenuation may reflect non-random BMI availability. | `bmi_complete_case_selection_tests.csv`; `bmi_missingness_by_enterotype.csv` | Limitations; Table S10 note | BMI missing 80/305; completeness differed by enterotype (P=0.0368); missing participants had higher mean UA (P=0.0491). |

Forbidden summary claims: “all findings remained significant,” “renal/metabolic confounding was excluded,” or “attenuation was caused by BMI/creatinine.”
