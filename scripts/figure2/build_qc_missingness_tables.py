# -*- coding: utf-8 -*-
"""Build Table S18 (QC summary) and Table S19 (missingness & per-analysis n)
for Minor Comments 1 and 2 (second-round revision).

Sources:
- data/clean/kneaddata_read_count_table.tsv (per-sample read counts, 325 rows)
- results/tables/section_3_3/section3_6_analysis_data.csv (305 analysis-sample IDs)
- results/tables/section_3_3/missingness_summary.csv
- results/tables/section_3_3/model_complete_case_counts.csv
- round-1 BMI missingness CSVs (per-enterotype BMI completeness)

QC assertions verify computed stats against Methods 2.2 text numbers
(mean raw 4.6e7, range 3.9e7-1.2e8; mean post-filter 4.2e7, range 2.2e7-1.1e8).
"""
import csv
import os
import statistics

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT_DIR = os.path.join(BASE, "results", "tables", "final")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- load 305 analysis sample IDs ----------
id_file = os.path.join(BASE, "results", "tables", "section_3_3", "section3_6_analysis_data.csv")
with open(id_file, encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    id_field = next((name for name in ("SampleID", "sampleid", "Sample_ID")
                     if name in reader.fieldnames), reader.fieldnames[0])
    ids = [row[id_field] for row in reader]
assert len(ids) == 305, f"expected 305 IDs, got {len(ids)}"
id_set = set(ids)

# ---------- load kneaddata table, restrict to 305 ----------
kn_file = os.path.join(BASE, r"data\clean\kneaddata_read_count_table.tsv")
rows = []
with open(kn_file, encoding="utf-8") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        rows.append(parts)
print(f"kneaddata rows total: {len(rows)}")
kn_ids = {r[idx["Sample"]] for r in rows}
missing_in_kn = id_set - kn_ids
assert not missing_in_kn, f"analysis samples missing from kneaddata table: {missing_in_kn}"
rows305 = [r for r in rows if r[idx["Sample"]] in id_set]
print(f"rows restricted to analysis cohort: {len(rows305)}")
assert len(rows305) == 305
extra = len(rows) - 305
print(f"non-cohort rows excluded: {extra}")


def val(r, name):
    return float(r[idx[name]])


raw_totals, final_totals, host_fracs = [], [], []
for r in rows305:
    raw = val(r, "raw pair1") + val(r, "raw pair2")
    final = (val(r, "final pair1") + val(r, "final pair2")
             + val(r, "final orphan1") + val(r, "final orphan2"))
    raw_totals.append(raw)
    final_totals.append(final)
    host_removed = raw - final
    host_fracs.append(host_removed / raw)


def stats(xs):
    return dict(mean=statistics.mean(xs), median=statistics.median(xs),
                mn=min(xs), mx=max(xs))


s_raw = stats(raw_totals)
s_fin = stats(final_totals)
retained = [f_ / r_ for f_, r_ in zip(final_totals, raw_totals)]
s_ret = stats(retained)

print("\n--- computed over 305 analysis samples ---")
print(f"raw reads:   mean {s_raw['mean']:.3e}  median {s_raw['median']:.3e}  "
      f"range {s_raw['mn']:.3e}-{s_raw['mx']:.3e}")
print(f"final reads: mean {s_fin['mean']:.3e}  median {s_fin['median']:.3e}  "
      f"range {s_fin['mn']:.3e}-{s_fin['mx']:.3e}")
print(f"retained fraction: mean {s_ret['mean']:.3f}  range {s_ret['mn']:.3f}-{s_ret['mx']:.3f}")

# ---------- QC assertions against Methods 2.2 text ----------
# text: mean raw 4.6e7 (range 3.9e7-1.2e8); mean post-filter 4.2e7 (range 2.2e7-1.1e8)
assert abs(s_raw["mean"] - 4.6e7) / 4.6e7 < 0.02, f"raw mean mismatch: {s_raw['mean']:.3e}"
assert abs(s_fin["mean"] - 4.2e7) / 4.2e7 < 0.02, f"final mean mismatch: {s_fin['mean']:.3e}"
assert abs(s_raw["mn"] - 3.9e7) / 3.9e7 < 0.05, f"raw min mismatch: {s_raw['mn']:.3e}"
assert abs(s_raw["mx"] - 1.2e8) / 1.2e8 < 0.05, f"raw max mismatch: {s_raw['mx']:.3e}"
assert abs(s_fin["mn"] - 2.2e7) / 2.2e7 < 0.05, f"final min mismatch: {s_fin['mn']:.3e}"
assert abs(s_fin["mx"] - 1.1e8) / 1.1e8 < 0.05, f"final max mismatch: {s_fin['mx']:.3e}"
print("QC assertions vs Methods 2.2 numbers: PASS")


def fmt_e(x):
    """format like 4.6 x 10^7"""
    return f"{x:.2e}".replace("e+0", " x 10^").replace("e+", " x 10^")


# ---------- Table S18: QC summary ----------
s18_path = os.path.join(OUT_DIR, "Table_S1_qc_summary.csv")
with open(s18_path, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["Dataset", "Item", "Value"])
    w.writerow(["A. Shotgun metagenomics", "Samples sequenced and analysed", "305"])
    w.writerow(["A. Shotgun metagenomics", "Platform",
                "Illumina NovaSeq 6000, paired-end 150 bp"])
    w.writerow(["A. Shotgun metagenomics", "Raw reads per sample, mean (range)",
                f"{fmt_e(s_raw['mean'])} ({fmt_e(s_raw['mn'])}-{fmt_e(s_raw['mx'])})"])
    w.writerow(["A. Shotgun metagenomics",
                "High-quality non-host reads per sample after quality filtering and host-DNA removal, mean (range)",
                f"{fmt_e(s_fin['mean'])} ({fmt_e(s_fin['mn'])}-{fmt_e(s_fin['mx'])})"])
    w.writerow(["A. Shotgun metagenomics", "Reads retained after filtering, mean % (range)",
                f"{100*s_ret['mean']:.1f}% ({100*s_ret['mn']:.1f}%-{100*s_ret['mx']:.1f}%)"])
    w.writerow(["A. Shotgun metagenomics", "Quality / host filtering pipeline",
                "KneadData v0.10.0 (bioBakery3 workflow)"])
    w.writerow(["A. Shotgun metagenomics", "Taxonomic / functional profiling",
                "MetaPhlAn3; HUMAnN3 (MetaCyc pathways, KEGG Orthologs, EC numbers)"])
    w.writerow(["A. Shotgun metagenomics",
                "Features quantified before prevalence filtering",
                "509 species-level clades; 474 MetaCyc pathways"])
    w.writerow(["A. Shotgun metagenomics",
                "Prevalence filter before differential-abundance testing",
                "features present in <10% of samples removed"])
    w.writerow(["B. Serum untargeted metabolomics", "Samples profiled", "303 of 305"])
    w.writerow(["B. Serum untargeted metabolomics", "Platform",
                "Vanquish UHPLC coupled to Orbitrap Exploris 120; HILIC (BEH Amide); positive and negative ion modes"])
    w.writerow(["B. Serum untargeted metabolomics", "QC design",
                "Pooled QC sample (equal aliquots of all samples) injected every 10 study samples"])
    w.writerow(["B. Serum untargeted metabolomics", "QC reproducibility filter",
                "Features with RSD > 30% in QC samples removed"])
    w.writerow(["B. Serum untargeted metabolomics", "Missing-value handling",
                "Overall missing rate < 5%; features missing in > 20% of samples excluded; remaining missing values imputed at half the per-feature minimum"])
    w.writerow(["B. Serum untargeted metabolomics", "Normalisation",
                "Isotopically labelled internal standards"])
    w.writerow(["B. Serum untargeted metabolomics", "Batch-effect assessment",
                "PCA of QC samples showed no detectable batch effects or instrumental drift"])
    w.writerow(["B. Serum untargeted metabolomics", "Feature filtering",
                "Unannotated features removed before downstream analysis"])
    w.writerow(["B. Serum untargeted metabolomics", "Metabolite-identification confidence",
                "MSI confidence level 2 (in-house MS2 database, BiotreeDB)"])
    w.writerow(["B. Serum untargeted metabolomics", "Metabolites retained after QC", "829"])
print(f"written: {s18_path}")

# ---------- Table S19: missingness & per-analysis sample sizes ----------
ms_file = os.path.join(BASE, r"results\tables\section_3_3\missingness_summary.csv")
cc_file = os.path.join(BASE, r"results\tables\section_3_3\model_complete_case_counts.csv")
with open(ms_file, encoding="utf-8-sig") as f:
    ms = list(csv.DictReader(f))
with open(cc_file, encoding="utf-8-sig") as f:
    cc = {row["model"]: int(row["n_complete"]) for row in csv.DictReader(f)}

# assertions against known manuscript numbers
msd = {row["variable"]: (int(row["n_nonmissing"]), float(row["pct_nonmissing"])) for row in ms}
assert msd["BMI"][0] == 225 and msd["CRP"][0] == 302
assert cc == {"demo": 305, "bmi": 225, "renal_cr": 305, "renal_egfr": 305,
              "cardiomet": 305, "full_clinical": 223, "bp_extended": 174}

var_labels = [
    ("UA", "Serum uric acid (UA)"),
    ("Hyperuricemia", "Hyperuricemia status (HUA)"),
    ("Enterotype", "Enterotype assignment"),
    ("Gender", "Sex"),
    ("Age", "Age"),
    ("BMI", "Body mass index (BMI)"),
    ("Urea", "Serum urea"),
    ("Cr", "Serum creatinine (Cr)"),
    ("eGFR", "Estimated glomerular filtration rate (eGFR)"),
    ("GLU", "Fasting glucose (GLU)"),
    ("TG", "Triglycerides (TG)"),
    ("HDLC", "HDL cholesterol (HDL-C)"),
    ("LDLC", "LDL cholesterol (LDL-C)"),
    ("CRP", "C-reactive protein (CRP)"),
    ("SBP", "Systolic blood pressure (SBP)"),
    ("DBP", "Diastolic blood pressure (DBP)"),
]

s19_path = os.path.join(OUT_DIR, "Table_S2_missingness_and_analysis_n.csv")
with open(s19_path, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["Part", "Item", "Non-missing n", "Missing n", "Missing %", "Note"])
    for key, label in var_labels:
        n_ok, pct_ok = msd[key]
        note = ""
        if key in ("SBP", "DBP"):
            note = "excluded from the primary robustness analyses owing to substantial missingness"
        w.writerow(["A. Variable-level missingness (N = 305)", label,
                    n_ok, 305 - n_ok, f"{100 - pct_ok:.1f}", note])
    w.writerow(["A. Variable-level missingness (N = 305)",
                "Serum untargeted metabolomics", 303, 2, "0.7",
                "metabolome-dependent analyses use the matched subset"])
    for model, label in [
        ("demo", "M1 Demographic (sex, age)"),
        ("bmi", "M2 +BMI"),
        ("renal_cr", "M3 +Renal (creatinine, urea)"),
        ("cardiomet", "M4 +Cardiometabolic (+GLU, TG, HDL-C, LDL-C)"),
        ("renal_egfr", "S1 +eGFR (replacing creatinine)"),
        ("full_clinical", "S2 Full clinical (+BMI, +CRP)"),
    ]:
        w.writerow(["B. Analysis-specific sample sizes", label, cc[model], "", "", ""])
    w.writerow(["B. Analysis-specific sample sizes",
                "Nested-model refits (fixed populations)", "223 / 225", "", "",
                "full-clinical and BMI complete-case populations"])
    w.writerow(["B. Analysis-specific sample sizes",
                "Subgroup analyses (males / females / non-obese / eGFR >= 90)",
                "81 / 224 / 218 / 271", "", "", ""])
    w.writerow(["B. Analysis-specific sample sizes",
                "Mediation models (Type 1 / Type 4 vs pooled Types 2+3)",
                "264 / 185", "", "", "restricted to participants with matched metabolomics"])
    w.writerow(["B. Analysis-specific sample sizes",
                "Within-Type 4 HUA+ vs HUA- comparison", "39 (12 / 27)", "", "", ""])
print(f"written: {s19_path}")
print("ALL DONE")
