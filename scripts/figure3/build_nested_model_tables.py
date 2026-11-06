# -*- coding: utf-8 -*-
"""
Build the two supplementary-table CSVs for Major Comment 3 from the VERIFIED
pipeline outputs in results/ (no re-fitting; repackaging only).

Outputs (results/tables/final/):
  Table_S6_nested_models.csv
  Table_S7_model_stability.csv
"""
import csv
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(BASE, "..", ".."))
RES = os.path.join(ROOT, "results", "tables", "section_3_3", "nested_model_stability")
OUT = os.path.join(ROOT, "results", "tables", "final")
os.makedirs(OUT, exist_ok=True)

Z = 1.959964  # 97.5% normal quantile; matches the Wald-type CIs reported in the manuscript


def fmt_p(p):
    p = float(p)
    return f"{p:.4f}" if p >= 0.0001 else f"{p:.2e}"


def fmt_beta(x):
    return f"{float(x):.1f}"


def fmt_or(x):
    return f"{float(x):.2f}"


def read_csv(name):
    with open(os.path.join(RES, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


summary = read_csv("02_nested_models_summary.csv")      # overall P per panel x model
ua = read_csv("03_nested_models_ua_contrasts.csv")      # UA contrasts (estimate, se)
hua = read_csv("04_nested_models_hua_contrasts.csv")    # HUA contrasts (OR, CI, p)
epv = read_csv("05_epv_by_model.csv")                   # events per variable
firth = read_csv("06_firth_hua_models.csv")             # Firth refits
register = read_csv("00_model_register.csv")

MODEL_LABEL = {
    "M1_Demographic": "M1 Demographic",
    "M2_BMI": "M2 +BMI",
    "M3_Renal_Cr": "M3 +Renal (creatinine, urea)",
    "M4_Cardiomet": "M4 +Cardiometabolic (GLU, TG, HDL-C, LDL-C)",
    "S2_Full_Clinical": "S2 Full clinical (+BMI, CRP)",
}
COV_ADDED = {
    "M1_Demographic": "sex, age (baseline)",
    "M2_BMI": "+ BMI",
    "M3_Renal_Cr": "+ creatinine, urea",
    "M4_Cardiomet": "+ glucose, triglycerides, HDL-C, LDL-C",
    "S2_Full_Clinical": "+ BMI, CRP (full clinical set)",
}
PANEL_LABEL = {
    "FULL223": "Fixed complete-case population (n = 223)",
    "BMI225": "Fixed BMI complete-case population (n = 225)",
}
PANEL_ORDER = {"FULL223": 1, "BMI225": 2}
MODEL_ORDER = {"M1_Demographic": 1, "M2_BMI": 2, "M3_Renal_Cr": 3,
               "M4_Cardiomet": 4, "S2_Full_Clinical": 5}
TERMS = ["Type2_vs_Type1", "Type3_vs_Type1", "Type4_vs_Type1"]

# ---------------------------------------------------------------- Table S14 --
key = lambda r: (r["panel"], r["model"])
sum_by = {(r["panel"], r["model"]): r for r in summary}
ua_by = {(r["panel"], r["model"], r["term"]): r for r in ua}
hua_by = {(r["panel"], r["model"], r["term"]): r for r in hua}

header = (
    ["panel", "model", "covariates_added", "n", "hua_events",
     "ua_overall_enterotype_p"]
    + [f"ua_{t.lower()}_{x}" for t in TERMS for x in ("beta", "ci_low", "ci_high")]
    + ["hua_overall_enterotype_p"]
    + [f"hua_{t.lower()}_{x}" for t in TERMS for x in ("or", "ci_low", "ci_high", "p")]
)

rows = []

# Anchor row: full-cohort M1 as published (manuscript Sections 3.2; overall P
# values re-verified by the fidelity rerun in 01_fidelity_rerun_summary.csv).
anchor = {
    "panel": "Full cohort (reference; M1 as published)",
    "model": "M1 Demographic",
    "covariates_added": "sex, age (baseline)",
    "n": "305", "hua_events": "64",
    "ua_overall_enterotype_p": "0.0372",
    "ua_type2_vs_type1_beta": "-20.0", "ua_type2_vs_type1_ci_low": "-40.2", "ua_type2_vs_type1_ci_high": "0.2",
    "ua_type3_vs_type1_beta": "-33.4", "ua_type3_vs_type1_ci_low": "-57.3", "ua_type3_vs_type1_ci_high": "-9.6",
    "ua_type4_vs_type1_beta": "-3.5", "ua_type4_vs_type1_ci_low": "-36.7", "ua_type4_vs_type1_ci_high": "29.7",
    "hua_overall_enterotype_p": "0.0106",
    "hua_type2_vs_type1_or": "0.46", "hua_type2_vs_type1_ci_low": "0.22", "hua_type2_vs_type1_ci_high": "0.98", "hua_type2_vs_type1_p": "0.045",
    "hua_type3_vs_type1_or": "0.23", "hua_type3_vs_type1_ci_low": "0.07", "hua_type3_vs_type1_ci_high": "0.70", "hua_type3_vs_type1_p": "0.010",
    "hua_type4_vs_type1_or": "1.01", "hua_type4_vs_type1_ci_low": "0.44", "hua_type4_vs_type1_ci_high": "2.33", "hua_type4_vs_type1_p": "0.98",
}
rows.append(anchor)

nested_keys = sorted(sum_by, key=lambda k: (PANEL_ORDER[k[0]], MODEL_ORDER[k[1]]))
for panel, model in nested_keys:
    s = sum_by[(panel, model)]
    row = {
        "panel": PANEL_LABEL[panel],
        "model": MODEL_LABEL[model],
        "covariates_added": COV_ADDED[model],
        "n": s["n"], "hua_events": s["hua_events"],
        "ua_overall_enterotype_p": fmt_p(s["ua_enterotype_p"]),
        "hua_overall_enterotype_p": fmt_p(s["hua_enterotype_p"]),
    }
    for t in TERMS:
        u = ua_by[(panel, model, t)]
        est, se = float(u["estimate"]), float(u["se"])
        row[f"ua_{t.lower()}_beta"] = fmt_beta(est)
        row[f"ua_{t.lower()}_ci_low"] = fmt_beta(est - Z * se)
        row[f"ua_{t.lower()}_ci_high"] = fmt_beta(est + Z * se)
        h = hua_by[(panel, model, t)]
        row[f"hua_{t.lower()}_or"] = fmt_or(h["OR"])
        row[f"hua_{t.lower()}_ci_low"] = fmt_or(h["ci_low"])
        row[f"hua_{t.lower()}_ci_high"] = fmt_or(h["ci_high"])
        row[f"hua_{t.lower()}_p"] = fmt_p(h["p"])
    rows.append(row)

out1 = os.path.join(OUT, "Table_S6_nested_models.csv")
with open(out1, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=header)
    w.writeheader()
    w.writerows(rows)

# ---------------------------------------------------------------- Table S15 --
firth_model = {}
for r in firth:
    m = r["model"]
    if m not in firth_model:
        firth_model[m] = {"plrt": r["enterotype_plrt_p"], "glm": r["enterotype_glm_lrt_p"],
                          "n": r["n"], "events": r["events"]}
    firth_model[m][r["term"]] = r

epv_by = {r["model"]: r for r in epv}

header15 = (["model", "n", "hua_events", "n_predictors", "epv",
             "standard_lrt_p_overall_enterotype",
             "firth_penalized_lrt_p_overall_enterotype"]
            + [f"firth_{t.lower()}_{x}" for t in TERMS for x in ("or", "ci_low", "ci_high")])

rows15 = []
for m in ["M1_Demographic", "M2_BMI", "M3_Renal_Cr", "M4_Cardiomet", "S1_Renal_eGFR", "S2_Full_Clinical"]:
    e = epv_by[m]
    fm = firth_model[m]
    row = {
        "model": MODEL_LABEL.get(m, "S1 +eGFR (replaces creatinine)"),
        "n": e["n"], "hua_events": e["events"],
        "n_predictors": e["predictors"], "epv": f"{float(e['epv']):.1f}",
        "standard_lrt_p_overall_enterotype": fmt_p(fm["glm"]),
        "firth_penalized_lrt_p_overall_enterotype": fmt_p(fm["plrt"]),
    }
    for t in TERMS:
        fr = fm[t]
        row[f"firth_{t.lower()}_or"] = fmt_or(fr["OR"])
        row[f"firth_{t.lower()}_ci_low"] = fmt_or(fr["ci_low"])
        row[f"firth_{t.lower()}_ci_high"] = fmt_or(fr["ci_high"])
    rows15.append(row)

out2 = os.path.join(OUT, "Table_S7_model_stability.csv")
with open(out2, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=header15)
    w.writeheader()
    w.writerows(rows15)

# -------------------------------------------------------------------- QC ----
assert len(rows) == 8, f"Table S14 must have 8 rows (1 anchor + 5 + 2), got {len(rows)}"
assert len(rows15) == 6, f"Table S15 must have 6 rows, got {len(rows15)}"
# spot checks against the verified sources
s2 = sum_by[("FULL223", "S2_Full_Clinical")]
assert abs(float(s2["ua_enterotype_p"]) - 0.507062604382589) < 1e-12
assert fmt_p(s2["ua_enterotype_p"]) == "0.5071"
assert rows[0]["ua_overall_enterotype_p"] == "0.0372"
m1_223 = [r for r in rows if r["panel"].startswith("Fixed complete-case") and r["model"] == "M1 Demographic"][0]
assert m1_223["ua_overall_enterotype_p"] == "0.0759"
assert m1_223["ua_type3_vs_type1_beta"] == "-38.4"
s2_row = rows15[-1]
assert s2_row["epv"] == "3.2" and s2_row["firth_penalized_lrt_p_overall_enterotype"] == "0.0147"
# BMI225 M2 row must equal the published M2 (fidelity anchor of the original model)
m2_225 = [r for r in rows if r["panel"].startswith("Fixed BMI") and r["model"] == "M2 +BMI"][0]
assert m2_225["ua_overall_enterotype_p"] == "0.1328" and m2_225["n"] == "225"
print("QC PASS")
print(f"Wrote {out1}  ({len(rows)} rows x {len(header)} cols)")
print(f"Wrote {out2}  ({len(rows15)} rows x {len(header15)} cols)")
