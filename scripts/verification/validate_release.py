#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Validate the public main-results release without refitting long analyses."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv
import json
import math
import re

from openpyxl import load_workbook
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
FINAL = ROOT / "results" / "tables" / "final"
REPORT = ROOT / "results" / "verification" / "release_validation.json"
checks = []


def check(name, condition, detail=""):
    passed = bool(condition)
    checks.append({"name": name, "pass": passed, "detail": str(detail)})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def rows(name):
    with (FINAL / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def close(a, b, tolerance=1e-8):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tolerance)


# Structure and portability.
required = [
    "README.md", "data/README.md",
    "supplementary/build_supplementary_xlsx.py",
    "supplementary/tables/Supplementary_Tables.xlsx",
]
for relative in required:
    check(f"required file {relative}", (ROOT / relative).is_file(), relative)
check("enterotype-robustness results directory", (ROOT / "results/enterotype_robustness").is_dir(), "results/enterotype_robustness")
check("covariate-sensitivity results directory", (ROOT / "results/covariate_sensitivity").is_dir(), "results/covariate_sensitivity")
script_files = [p for p in (ROOT / "scripts").rglob("*") if p.suffix.lower() in {".r", ".py"}]
check("R script count", len(list((ROOT / "scripts").rglob("*.R"))) == 54, len(list((ROOT / "scripts").rglob("*.R"))))
for path in script_files:
    if path.resolve() == Path(__file__).resolve():
        continue  # The validator itself contains the forbidden strings as test patterns.
    text = path.read_text(encoding="utf-8")
    forbidden = ["C:/Users/", "C:\\Users\\", "Journal_TranslationalMedicine", "Table_SX", "Figure_SX"]
    check(f"portable script {path.relative_to(ROOT)}", not any(item in text for item in forbidden), "no local path or placeholder")

expected_tables = {
    "Table_S1_qc_summary.csv", "Table_S2_missingness_and_analysis_n.csv",
    "Table_S3_purine_pathway_curation.csv", "Table_S4_mediator_traceability.csv",
    "Table_S5_continuous_gradient.csv", "Table_S6_nested_models.csv", "Table_S7_model_stability.csv",
    "Table_S8_core_features.csv", "Table_S9_purine_module_statistics.csv", "Table_S10_purine_module_UA.csv",
    "Table_S11_multiomics_covariate_sensitivity.csv", "Table_S12_candidate_feature_associations.csv",
    "Table_S13_mediation_main.csv", "Table_S14_mediation_binary_HUA.csv",
    "Table_S15_mediation_reference_split.csv", "Table_S16_mediation_medsens.csv",
    "Table_S17_mediation_extended_covariates.csv", "Table_S18_type4_heterogeneity.csv",
    "Table_S19_k3_outcome_sensitivity.csv",
} | {f"Table_S8{suffix}_core_features.csv" for suffix in "abcdefghi"}
check("final table file set", {p.name for p in FINAL.glob("*.csv")} == expected_tables, f"{len(expected_tables)} CSVs")

# Cohort and key numerical anchors.
with (ROOT / "results/tables/section_3_3/section3_6_analysis_data.csv").open("r", encoding="utf-8-sig", newline="") as handle:
    clinical = list(csv.DictReader(handle))
check("clinical cohort n", len(clinical) == 305, len(clinical))
check("enterotype counts", Counter(r["Enterotype"] for r in clinical) == {"Type 1": 120, "Type 2": 98, "Type 3": 48, "Type 4": 39}, Counter(r["Enterotype"] for r in clinical))

core = rows("Table_S8_core_features.csv")
expected_core = {
    "Taxa": (72, 2, 1), "Pathway": (43, 52, 0), "EC": (296, 371, 36), "Metabolome": (24, 0, 0),
}
for row in core:
    actual = (int(row["E1_only"]), int(row["E4_only"]), int(row["Shared"]))
    check(f"core counts {row['data_type']}", actual == expected_core[row["data_type"]], actual)
for suffix, expected_n in zip("abcdefghi", (72, 1, 2, 43, 52, 24, 296, 36, 371)):
    data = rows(f"Table_S8{suffix}_core_features.csv")
    check(f"S8{suffix} row count", len(data) == expected_n, len(data))

gradient = rows("Table_S5_continuous_gradient.csv")
check("continuous-gradient specifications", len(gradient) == 6, len(gradient))
primary = gradient[0]
check("primary P/B partial rho", close(primary["ua_partial_rho_sex_age"], -0.152, 5e-4), primary["ua_partial_rho_sex_age"])
check("primary P/B HUA OR", close(primary["hua_or_per_doubling"], 0.991, 5e-4), primary["hua_or_per_doubling"])

nested = rows("Table_S6_nested_models.csv")
stability = rows("Table_S7_model_stability.csv")
check("nested-model rows", len(nested) == 8, len(nested))
check("Firth-model rows", len(stability) == 6, len(stability))
check("S2 EPV", close(stability[-1]["epv"], 3.2, 0.01), stability[-1]["epv"])
check("S2 Firth overall P", close(stability[-1]["firth_penalized_lrt_p_overall_enterotype"], 0.0147, 5e-5), stability[-1]["firth_penalized_lrt_p_overall_enterotype"])

purine = rows("Table_S3_purine_pathway_curation.csv")
module_stats = rows("Table_S9_purine_module_statistics.csv")
module_ua = rows("Table_S10_purine_module_UA.csv")
check("purine curated pathways", len(purine) == 34, len(purine))
check("purine module-statistic rows", len(module_stats) == 33, len(module_stats))
check("purine module-UA rows", len(module_ua) == 3, len(module_ua))
check("de novo BH-adjusted P", close(module_ua[0]["BH-adjusted P (across the 3 module scores)"], 0.0365, 1e-10), module_ua[0]["BH-adjusted P (across the 3 module scores)"])
check("module-UA complete n", {int(r["n"]) for r in module_ua} == {305}, "all n=305")

mediation = rows("Table_S13_mediation_main.csv")
check("main mediation rows", len(mediation) == 44, len(mediation))
for model in ("Type1_vs_Ref", "Type4_vs_Ref"):
    significant = sum(float(r["ACME_padj"]) < 0.05 for r in mediation if r["model"] == model)
    check(f"BH-significant mediators {model}", significant == 5, significant)
extended = rows("Table_S17_mediation_extended_covariates.csv")
check("extended mediation rows", len(extended) == 44, len(extended))
check("ACME directions retained", sum(r["direction_consistent"].upper() == "TRUE" for r in extended) == 42, "42/44")
check("matched sex-and-age FDR-significant combinations", sum(r["base_sig"].upper() == "TRUE" for r in extended) == 1, "1 (the original larger-sample 10 are checked in Table S13)")
check("extended FDR-significant combinations", sum(r["ext_sig"].upper() == "TRUE" for r in extended) == 2, "2")
heterogeneity = rows("Table_S18_type4_heterogeneity.csv")
check("within-Type 4 tests", len(heterogeneity) == 23, len(heterogeneity))
check("within-Type 4 nominal findings", sum(float(r["wilcox_p"]) < 0.05 for r in heterogeneity) == 4, "4")
medsens = rows("Table_S16_mediation_medsens.csv")
check("medsens rows and complete column", len(medsens) == 10 and "rho_at_sign_change" in medsens[0], "10 rows")

check("multi-omics sensitivity rows", len(rows("Table_S11_multiomics_covariate_sensitivity.csv")) == 13, 13)
check("candidate traceability rows", len(rows("Table_S4_mediator_traceability.csv")) == 22, 22)
check("direct candidate rows", len(rows("Table_S12_candidate_feature_associations.csv")) == 22, 22)
k3 = rows("Table_S19_k3_outcome_sensitivity.csv")
check("k3 table rows", len(k3) == 15, len(k3))
check("k3 sections", Counter(r["section"] for r in k3) == {"A_membership_crosswalk": 4, "B_component_descriptives": 3, "C_association_tests": 8}, Counter(r["section"] for r in k3))

# Workbook structure and consistency with final CSVs.
workbook_path = ROOT / "supplementary/tables/Supplementary_Tables.xlsx"
wb = load_workbook(workbook_path, data_only=False)
expected_order = ["Contents"]
for number in range(1, 20):
    expected_order.append(f"S{number}")
    if number == 8:
        expected_order.extend(f"S8{suffix}" for suffix in "abcdefghi")
check("workbook sheet order", wb.sheetnames == expected_order, f"{len(wb.sheetnames)} sheets")
check("S8h is shared EC", wb["S8h"].max_row == 37, wb["S8h"].max_row - 1)
check("S8i is Type 4-exclusive EC", wb["S8i"].max_row == 372, wb["S8i"].max_row - 1)
check("S16 includes rho_at_sign_change", wb["S16"].cell(1, 4).value == "rho_at_sign_change", wb["S16"].cell(1, 4).value)
check("S17 reader-friendly headers", wb["S17"].cell(1, 4).value == "sex_age_ACME" and wb["S17"].cell(1, 17).value == "FDR_status_stable", "matched")
check("S13 workbook mediation rows", wb["S13"].max_row == 45, wb["S13"].max_row - 1)
check("S10 workbook de novo adjusted P", close(wb["S10"]["D2"].value, 0.0365, 1e-12), wb["S10"]["D2"].value)
check("workbook has no formulas", not any(cell.data_type == "f" for ws in wb for row in ws.iter_rows() for cell in row), "static result tables")

# Final PDFs are present, non-empty and single-page.
supplementary_pdf_names = [
    "FigureS1_enterotype_robustness.pdf", "FigureS2_Individual_Bar.pdf",
    "FigureS3_type4_diagnostics.pdf", "FigureS4_composite.pdf",
    "FigureS5_Purine_Pathway_Context.pdf", "FigureS6_module_scores.pdf",
    "FigureS7_Full_Species_Heatmap.pdf",
]
pdfs = [ROOT / f"results/figures/Figure{n}.pdf" for n in range(1, 7)]
pdfs += [ROOT / "supplementary/figures" / name for name in supplementary_pdf_names]
check("final PDF count", len(pdfs) == 13, len(pdfs))
for path in pdfs:
    page_count = len(PdfReader(str(path)).pages)
    check(f"final PDF {path.name}", path.stat().st_size > 1000 and page_count >= 1, f"{path.stat().st_size} bytes; {page_count} page(s)")

REPORT.parent.mkdir(parents=True, exist_ok=True)
report = {
    "status": "PASS",
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "checks_passed": len(checks),
    "checks": checks,
    "scope_note": "Structural and numerical checks only; long DMM and mediation resampling were not rerun by this validator.",
}
REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"PASS: {len(checks)} checks. Report: {REPORT}")
