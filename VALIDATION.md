# Validation Record

## Completed checks

### Script and path checks

- All **54 R scripts** passed `parse()` syntax checking.
- Python table builders ran successfully from the portable repository layout.
- Public scripts contain no local absolute paths, `Table_SX` / `Figure_SX` placeholders, or dependencies on the manuscript-revision directories.

### Final CSVs and workbook

- The continuous-gradient, QC/missingness, nested-model/Firth, and purine-module Python builders regenerated their final CSVs successfully.
- Regenerated tables matched the frozen reference CSVs; the S8h/S8i assignments were explicitly checked against the submitted workbook: **S8h = shared ECs, S8i = Type 4-exclusive ECs**.
- The rebuilt workbook contains **29 sheets**: Contents, S1–S19, and S8a–S8i.
- Every shared data cell in the 28 data sheets matched the submitted workbook. Historical title rows in S8 sheets were excluded from the data-cell comparison. S16 intentionally retains the complete `rho_at_sign_change` column required by the final table description.
- The workbook compatibility gate reported **0 formulas, 0 missing caches, 0 errors, and 0 warnings**. The delivered workbook contains static result values rather than formulas.

### Independent R reruns of the new second-round analyses

- `scripts/figure3/nested_model_stability.R`: **23/23 internal QC assertions passed**, including cohort structure, original-model fidelity, nested-model checks, and Firth-model checks.
- `scripts/figure4/purine_module_analysis.R`: **75/75 internal QC assertions passed**, including all 474 pathway identifiers, the frozen 34-pathway set, EC-union calculations, denominator provenance, dual implementations, model contrasts, and source-table fidelity.
- Key regenerated estimates and final CSVs agreed with the frozen analysis outputs.

### Frozen long-running analyses

The 5,000-replicate mediation analyses, 100 DMM bootstrap refits, and 100 DMM initialization refits were not repeated during release assembly. Their authoritative outputs are retained, and structural/key-number checks are implemented in `scripts/verification/validate_release.py`.

## R runtime caveat observed during this release check

The manuscript analyses were originally run in **R 4.3.1**. That exact runtime was not available on the release-check machine. The available R 4.5.1 installation loaded packages built under R 4.5.2 and 4.5.3. Both new second-round analyses completed their calculations, wrote their outputs, printed their final success messages, and passed all internal assertions; however, the mismatched R runtime then raised a segmentation fault during interpreter shutdown. Simple package-load tests reproduced the shutdown fault independently of the analysis scripts.

Accordingly, this release does **not** claim an error-free end-to-end rerun under R 4.5.1. Exact reproduction should use R 4.3.1 with compatible package versions, guided by the relevant `session_info.txt` records. The shutdown fault is an environment limitation, not evidence of a statistical discrepancy.

## Reproducing the release checks

From the repository root:

```bash
python scripts/verification/validate_release.py
```

The validator checks required files, portable paths, final-table coverage, cohort and core-feature counts, key numerical anchors, mediation and sensitivity summaries, workbook structure, and all 13 final PDFs. It writes `results/verification/release_validation.json`.

`MANIFEST.sha256` records hashes for the exact file set prepared for publication.
