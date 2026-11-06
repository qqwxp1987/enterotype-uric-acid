# -*- coding: utf-8 -*-
"""Build Tables S5 and S19 from the frozen enterotype-robustness outputs.

The continuous-gradient and k=3 sensitivity statistics are read from the
authoritative result CSVs. Hard assertions prevent silent numerical drift."""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "tables" / "supplementary"
OUT = ROOT / "results" / "tables" / "final"
OUT.mkdir(parents=True, exist_ok=True)


def approx(a, b, tol=1e-6):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(p, fieldnames, rows):
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def fmt(x, nd):
    return f"{float(x):.{nd}f}"


def fmt_p(x):
    x = float(x)
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}".rstrip("0").rstrip(".")


# --------------------------------------------------------------------------
# Table S16: continuous Prevotella/Bacteroides gradient association tests
# --------------------------------------------------------------------------
pb = read_csv(SRC / "pb_ratio_sensitivity.csv")
assert len(pb) == 6, f"expected 6 specifications, got {len(pb)}"

spec_labels = {
    "pc_1e_6": "Pseudocount 1e-6 (primary)",
    "pc_1e_5": "Pseudocount 1e-5",
    "pc_1e_4": "Pseudocount 1e-4",
    "pc_1e_3": "Pseudocount 1e-3",
    "half_min_positive": "Half-minimum positive replacement",
    "prevotella_detected": "Prevotella-detected participants only",
}

# --- QC on frozen values quoted in the response letter ---
r0 = pb[0]
assert r0["specification"] == "pc_1e_6"
assert approx(r0["partial_rho"], -0.151825421762326, 1e-9)
assert approx(r0["partial_p"], 0.00811495095607351, 1e-9)
assert approx(r0["hua_or"], 0.990568518098423, 1e-9)
assert approx(r0["hua_ci_lower"], 0.950324803993722, 1e-9)
assert approx(r0["hua_ci_upper"], 1.03251644587632, 1e-9)
assert approx(r0["hua_p"], 0.654289046443617, 1e-9)
assert approx(r0["nonlinear_ns3_p"], 0.864582663780821, 1e-9)
assert approx(r0["nonlinear_ns4_p"], 0.347593858672415, 1e-9)
assert approx(r0["aic_enterotype"], 3515.87234428002, 1e-6)
assert approx(r0["aic_pb_linear"], 3516.12582305862, 1e-6)
assert approx(r0["prevotella_zero_n"], 122, 1e-9)
det = [r for r in pb if r["specification"] == "prevotella_detected"][0]
assert approx(det["n_ua"], 183, 1e-9)
assert approx(det["partial_rho"], -0.102319266145344, 1e-9)
assert approx(det["partial_p"], 0.170490964801156, 1e-9)

s16_rows = []
for r in pb:
    s16_rows.append({
        "specification": spec_labels[r["specification"]],
        "n_ua_analysis": int(float(r["n_ua"])),
        "n_hua_analysis": int(float(r["n_hua"])),
        "prevotella_undetected_n": int(float(r["prevotella_zero_n"])),
        "ua_partial_rho_sex_age": fmt(r["partial_rho"], 3),
        "ua_partial_p": fmt_p(r["partial_p"]),
        "nonlinearity_spline_df3_p": fmt(r["nonlinear_ns3_p"], 3),
        "nonlinearity_spline_df4_p": fmt(r["nonlinear_ns4_p"], 3),
        "hua_or_per_doubling": fmt(r["hua_or"], 3),
        "hua_or_ci_low": fmt(r["hua_ci_lower"], 3),
        "hua_or_ci_high": fmt(r["hua_ci_upper"], 3),
        "hua_p": fmt(r["hua_p"], 3),
        "aic_enterotype_k4": fmt(r["aic_enterotype"], 1),
        "aic_pb_linear": fmt(r["aic_pb_linear"], 1),
        "aic_pb_spline_df3": fmt(r["aic_pb_spline3"], 1),
        "aic_pb_spline_df4": fmt(r["aic_pb_spline4"], 1),
    })

s16_fields = ["specification", "n_ua_analysis", "n_hua_analysis",
              "prevotella_undetected_n", "ua_partial_rho_sex_age",
              "ua_partial_p", "nonlinearity_spline_df3_p",
              "nonlinearity_spline_df4_p",
              "hua_or_per_doubling", "hua_or_ci_low", "hua_or_ci_high",
              "hua_p", "aic_enterotype_k4", "aic_pb_linear",
              "aic_pb_spline_df3", "aic_pb_spline_df4"]
write_csv(OUT / "Table_S5_continuous_gradient.csv", s16_fields, s16_rows)

# --------------------------------------------------------------------------
# Table S17: outcome sensitivity under the best-fitting k=3 partition
# --------------------------------------------------------------------------
cross = read_csv(SRC / "k4_k3_membership_crosswalk.csv")
desc = read_csv(SRC / "k3_outcome_descriptives.csv")
sens = read_csv(SRC / "k3_outcome_sensitivity.csv")

# --- QC on frozen values quoted in the response letter ---
cw = {r["k4_enterotype"]: r for r in cross}
assert int(cw["Type 1"]["C1"]) == 120 and int(cw["Type 2"]["C2"]) == 98
assert int(cw["Type 3"]["C3"]) == 36 and int(cw["Type 4"]["C3"]) == 39
assert int(cw["Type 3"]["C1"]) == 1 and int(cw["Type 3"]["C2"]) == 11
dd = {r["k3"]: r for r in desc}
assert int(dd["C1"]["n"]) == 121 and int(dd["C2"]["n"]) == 109 and int(dd["C3"]["n"]) == 75
assert int(dd["C3"]["HUA_n"]) == 15
sa_hua = [r for r in sens if r["adjustment"] == "sex_age"
          and r["outcome"] == "HUA" and r["contrast"] == "C1_vs_C2"][0]
assert approx(sa_hua["estimate"], 2.47959415751843, 1e-9)
assert approx(sa_hua["p_value"], 0.0176311837483939, 1e-9)
rc_hua = [r for r in sens if r["adjustment"] == "renal_cardiometabolic"
          and r["outcome"] == "HUA" and r["contrast"] == "C1_vs_C2"][0]
assert approx(rc_hua["estimate"], 3.3755360310461, 1e-9)
assert approx(rc_hua["p_value"], 0.00308933944407195, 1e-9)
c3_sa = [r for r in sens if r["adjustment"] == "sex_age"
         and r["outcome"] == "HUA" and r["contrast"] == "C3_vs_C2"][0]
assert approx(c3_sa["estimate"], 1.33052609118811, 1e-9)
assert approx(c3_sa["p_value"], 0.524284193101432, 1e-9)

adj_labels = {
    "sex_age": "Sex + age",
    "renal_cardiometabolic": ("Renal + cardiometabolic (creatinine, urea, "
                              "glucose, triglycerides, HDL-C, LDL-C)"),
}
outcome_labels = {"UA": "Serum UA (umol/L, GLS beta)",
                  "HUA": "Binary HUA (logistic OR)"}
contrast_labels = {
    "C1_vs_C2": "Type 1-like (C1) vs Type 2-like (C2)",
    "C3_vs_C2": "Merged Prevotella-oriented (C3) vs Type 2-like (C2)",
}

s17_fields = ["section", "k4_enterotype", "k3_component",
              "k3_C1_n", "k3_C2_n", "k3_C3_n",
              "n", "mean_ua", "sd_ua", "hua_n", "hua_pct",
              "adjustment", "outcome", "contrast", "estimate",
              "ci_low", "ci_high", "p_value", "global_p", "scale"]
s17_rows = []
for r in cross:
    s17_rows.append({
        "section": "A_membership_crosswalk",
        "k4_enterotype": r["k4_enterotype"],
        "k3_C1_n": r["C1"], "k3_C2_n": r["C2"], "k3_C3_n": r["C3"],
    })
for r in desc:
    s17_rows.append({
        "section": "B_component_descriptives",
        "k3_component": r["k3"],
        "n": r["n"], "mean_ua": fmt(r["mean_UA"], 1),
        "sd_ua": fmt(r["sd_UA"], 1), "hua_n": r["HUA_n"],
        "hua_pct": fmt(r["HUA_pct"], 1),
    })
for r in sens:
    is_or = r["scale"] == "odds_ratio"
    nd = 2 if is_or else 1
    s17_rows.append({
        "section": "C_association_tests",
        "adjustment": adj_labels[r["adjustment"]],
        "outcome": outcome_labels[r["outcome"]],
        "contrast": contrast_labels[r["contrast"]],
        "n": r["n"],
        "estimate": fmt(r["estimate"], nd),
        "ci_low": fmt(r["ci_lower"], nd),
        "ci_high": fmt(r["ci_upper"], nd),
        "p_value": fmt_p(r["p_value"]),
        "global_p": fmt_p(r["global_p"]),
        "scale": r["scale"],
    })
write_csv(OUT / "Table_S19_k3_outcome_sensitivity.csv", s17_fields, s17_rows)

print("QC assertions passed.")
print(f"Wrote {OUT / 'Table_S5_continuous_gradient.csv'}")
print(f"Wrote {OUT / 'Table_S19_k3_outcome_sensitivity.csv'}")
