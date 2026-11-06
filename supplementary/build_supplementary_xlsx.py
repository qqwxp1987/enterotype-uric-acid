#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build Additional file 8 (Tables S1-S19) from final, numbered CSV sources.

No statistical estimates are calculated here. The workbook contains static values
only; every value is read from results/tables/final/.
"""
from pathlib import Path
import csv
import re

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "tables" / "final"
OUTPUT = ROOT / "supplementary" / "tables" / "Supplementary_Tables.xlsx"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

FILES = {
    "S1": "Table_S1_qc_summary.csv",
    "S2": "Table_S2_missingness_and_analysis_n.csv",
    "S3": "Table_S3_purine_pathway_curation.csv",
    "S4": "Table_S4_mediator_traceability.csv",
    "S5": "Table_S5_continuous_gradient.csv",
    "S6": "Table_S6_nested_models.csv",
    "S7": "Table_S7_model_stability.csv",
    "S8": "Table_S8_core_features.csv",
    **{f"S8{suffix}": f"Table_S8{suffix}_core_features.csv" for suffix in "abcdefghi"},
    "S9": "Table_S9_purine_module_statistics.csv",
    "S10": "Table_S10_purine_module_UA.csv",
    "S11": "Table_S11_multiomics_covariate_sensitivity.csv",
    "S12": "Table_S12_candidate_feature_associations.csv",
    "S13": "Table_S13_mediation_main.csv",
    "S14": "Table_S14_mediation_binary_HUA.csv",
    "S15": "Table_S15_mediation_reference_split.csv",
    "S16": "Table_S16_mediation_medsens.csv",
    "S17": "Table_S17_mediation_extended_covariates.csv",
    "S18": "Table_S18_type4_heterogeneity.csv",
    "S19": "Table_S19_k3_outcome_sensitivity.csv",
}

DESCRIPTIONS = {
    "S1": ("Quality-control summary for the metagenomic and metabolomic datasets", "Sequencing and post-filtering statistics, metabolomics QC design, filtering, normalization, batch-effect assessment, and retained feature counts."),
    "S2": ("Missingness summary and analysis-specific sample sizes", "Variable-level missingness and sample sizes for each prespecified analysis."),
    "S3": ("Composition-based curation of purine-related pathways", "The 34 curated pathways with module assignment, curation role, and reaction counts; includes the PWY-6595 annotation note."),
    "S4": ("Concise selection traceability for the 22 upstream candidates", "Data layer, upstream source, inclusion basis, statistical support, effect-size role, bounded rationale, and reference or database support."),
    "S5": ("Continuous Prevotella/Bacteroides gradient association tests", "UA and HUA association tests under six zero-handling specifications, including nonlinearity and model-fit comparisons."),
    "S6": ("Nested-model comparison within fixed complete-case populations", "The adjustment sequence refitted in fixed complete-case populations, with UA coefficients and HUA odds ratios."),
    "S7": ("Model stability assessment for binary HUA", "Events per variable, standard likelihood-ratio tests, and Firth penalized refits for each HUA model."),
    "S8": ("Independent core multi-omics features by enterotype", "Layer-count summary; companion sheets S8a-S8i contain the per-layer core-feature lists."),
    "S9": ("Purine-module score statistics across enterotypes", "Group medians, overall enterotype effects, and six Holm-adjusted pairwise contrasts per module."),
    "S10": ("Associations between purine-module scores and serum uric acid", "Sex-and-age-adjusted partial Spearman correlations with nominal and BH-adjusted P-values across the three module scores."),
    "S11": ("Sensitivity of differential features and cross-omics associations to BMI and creatinine adjustment", "Matched complete-case comparison of differential features, cross-omics correlations, partial db-RDA, and selection diagnostics."),
    "S12": ("Candidate feature attribution and direct UA/HUA associations", "Attribution and direct UA/HUA associations for the 22 candidates under sex-and-age and sex-age-BMI-creatinine adjustment."),
    "S13": ("Mediation analysis - main results", "Bootstrap estimates for 22 candidates in Type 1 and exploratory Type 4 models, with 5,000 resamples and BH-adjusted ACME P-values."),
    "S14": ("Mediation sensitivity analysis - binary hyperuricemia outcome", "Mediation models refitted with binary HUA status instead of continuous serum UA."),
    "S15": ("Mediation sensitivity analysis - reference-group split", "Type 1 and Type 4 each compared separately against Type 2 and Type 3."),
    "S16": ("Mediation sensitivity analysis - sequential ignorability", "Sensitivity parameter rho at which the estimated ACME would change sign."),
    "S17": ("Mediation sensitivity to BMI and serum-creatinine adjustment", "Matched complete-case mediation models with 5,000 nonparametric bootstrap resamples."),
    "S18": ("Within-Type 4 heterogeneity - HUA+ versus HUA-", "Wilcoxon P-values and Cliff's delta for 22 candidates plus Prevotella copri; nominal, unadjusted exploratory tests."),
    "S19": ("Outcome sensitivity under the best-fitting k = 3 enterotype solution", "Membership crosswalk, component descriptives, and covariate-adjusted UA/HUA associations."),
}
for suffix, label in zip("abcdefghi", (
    "Type 1-exclusive core species", "Shared core species", "Type 4-exclusive core species",
    "Type 1-exclusive core pathways", "Type 4-exclusive core pathways", "Type 1 core serum metabolites",
    "Type 1-exclusive core enzymes", "Shared core enzymes", "Type 4-exclusive core enzymes",
)):
    DESCRIPTIONS[f"S8{suffix}"] = (label, f"Per-feature differential-abundance estimates and BH-adjusted P-values: {label}.")

ORDER = ["Contents"]
for number in range(1, 20):
    ORDER.append(f"S{number}")
    if number == 8:
        ORDER.extend(f"S8{suffix}" for suffix in "abcdefghi")

NUMBER = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")
def typed(value, blank_na=False):
    if value == "" or (blank_na and value == "NA"):
        return None
    if value in {"TRUE", "FALSE"}:
        return value == "TRUE"
    if NUMBER.fullmatch(value):
        try:
            return int(value) if re.fullmatch(r"[+-]?\d+", value) else float(value)
        except ValueError:
            pass
    # Prevent text annotations from being interpreted as spreadsheet formulas.
    return "'" + value if value.startswith("=") else value

wb = Workbook()
contents = wb.active
contents.title = "Contents"
contents.append(["Table", "Title", "Description"])
for sheet in ORDER[1:]:
    title, description = DESCRIPTIONS[sheet]
    contents.append([sheet, title, description])
    contents.cell(contents.max_row, 1).hyperlink = f"#'{sheet}'!A1"
    contents.cell(contents.max_row, 1).font = Font(color="0563C1", underline="single")

header_fill = PatternFill("solid", fgColor="1F4E78")
header_font = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
body_font = Font(name="Calibri", size=10)
for sheet in ORDER[1:]:
    source = SOURCE / FILES[sheet]
    if not source.is_file():
        raise FileNotFoundError(source)
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2:
        raise ValueError(f"Table has no data rows: {source}")
    ws = wb.create_sheet(sheet)
    for r_idx, row in enumerate(rows, 1):
        for c_idx, value in enumerate(row, 1):
            s17_headers = {
                "base_ACME": "sex_age_ACME", "base_CI_lo": "sex_age_CI_lo", "base_CI_hi": "sex_age_CI_hi",
                "base_p": "sex_age_p", "base_padj": "sex_age_FDR",
                "ext_ACME": "sex_age_BMI_Cr_ACME", "ext_CI_lo": "sex_age_BMI_Cr_CI_lo",
                "ext_CI_hi": "sex_age_BMI_Cr_CI_hi", "ext_p": "sex_age_BMI_Cr_p", "ext_padj": "sex_age_BMI_Cr_FDR",
                "base_sig": "sex_age_FDR_significant", "ext_sig": "sex_age_BMI_Cr_FDR_significant",
                "sig_stable": "FDR_status_stable",
            }
            header_value = s17_headers.get(value, value) if sheet == "S17" and r_idx == 1 else value
            cell = ws.cell(r_idx, c_idx, header_value if r_idx == 1 else typed(value, blank_na=(sheet == "S11")))
            cell.alignment = Alignment(vertical="top", wrap_text=True)
            cell.font = header_font if r_idx == 1 else body_font
            if r_idx == 1:
                cell.fill = header_fill
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(ws.max_column)}{len(rows)}"
    for c_idx in range(1, ws.max_column + 1):
        values = [str(ws.cell(r, c_idx).value or "") for r in range(1, ws.max_row + 1)]
        ws.column_dimensions[get_column_letter(c_idx)].width = min(48, max(10, max(map(len, values)) + 2))
    ws.sheet_properties.pageSetUpPr.fitToPage = True
    ws.page_setup.orientation = "landscape"
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0
    ws.print_title_rows = "1:1"

note = (SOURCE / "Table_S3_purine_pathway_curation_footnote.txt").read_text(encoding="utf-8").strip()
ws = wb["S3"]
note_row = ws.max_row + 2
ws.cell(note_row, 1, note)
ws.merge_cells(start_row=note_row, start_column=1, end_row=note_row, end_column=5)
ws.cell(note_row, 1).alignment = Alignment(wrap_text=True, vertical="top")
ws.cell(note_row, 1).font = Font(name="Calibri", size=10, italic=True)
ws.row_dimensions[note_row].height = 48

for ws in wb.worksheets:
    ws.freeze_panes = "A2"
    for cell in ws[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(vertical="top", wrap_text=True)
    for c_idx in range(1, ws.max_column + 1):
        letter = get_column_letter(c_idx)
        if ws.title == "Contents":
            widths = {1: 10, 2: 70, 3: 100}
            ws.column_dimensions[letter].width = widths.get(c_idx, 30)
    ws.auto_filter.ref = f"A1:{get_column_letter(ws.max_column)}{ws.max_row}"
wb._sheets = [wb[name] for name in ORDER]
wb.active = 0
wb.save(OUTPUT)
print(f"Wrote {OUTPUT}")
print(f"Sheets ({len(ORDER)}): {', '.join(ORDER)}")
