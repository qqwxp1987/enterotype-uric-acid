#!/usr/bin/env python
"""Build numbered Supplementary Tables S1-S6 for the HUA_FCTJ2021 manuscript.

Reads existing result CSVs (read-only) and compiles them into a single
multi-sheet workbook: supplementary/Supplementary_Tables.xlsx
"""
import os
import sys

ROOT = r"c:\Users\qqwxp\OneDrive\Work\HUA_FCTJ2021"
TBL = os.path.join(ROOT, "results", "tables")
OUT_DIR = os.path.join(ROOT, "supplementary")
OUT_XLSX = os.path.join(OUT_DIR, "Supplementary_Tables.xlsx")

try:
    import pandas as pd
except ImportError as e:
    print("MISSING_DEP:pandas", e)
    sys.exit(2)
try:
    import openpyxl  # noqa: F401
except ImportError as e:
    print("MISSING_DEP:openpyxl", e)
    sys.exit(2)


def read_csv(path):
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    if df.empty:
        return df
    return df


# ---------------------------------------------------------------------------
# Source file resolution
# ---------------------------------------------------------------------------
S1_summary = os.path.join(TBL, "section_3_4", "independent_core_summary.csv")
S2_main = os.path.join(TBL, "section_3_6", "mediation_main_results.csv")
S3_binary = os.path.join(TBL, "section_3_6", "mediation_binary_sensitivity.csv")
S4_refgroup = os.path.join(TBL, "section_3_6", "mediation_refgroup_sensitivity.csv")
S5_medsens = os.path.join(TBL, "section_3_6", "mediation_medsens_results.csv")
S6_type4 = os.path.join(TBL, "section_3_6", "type4_heterogeneity_results.csv")

# Per-layer core feature lists for S1 (optional, included if present)
S1_layer_files = []
for layer in ("taxa", "pathway", "ec", "metabolome"):
    for grp in ("E1_only", "E4_only", "Shared"):
        p = os.path.join(TBL, "section_3_4", f"{layer}_{grp}.csv")
        if os.path.isfile(p):
            S1_layer_files.append((layer, grp, p))

missing = []
empty = []

titles = {
    "Table S1": "Independent core multi-omics features by enterotype",
    "Table S2": "Causal mediation analysis - main results (Type 1 and Type 4 models)",
    "Table S3": "Mediation sensitivity analysis - binary hyperuricemia outcome",
    "Table S4": "Mediation sensitivity analysis - reference-group split (Type 2 vs Type 3)",
    "Table S5": "Mediation sensitivity analysis - sequential-ignorability (medsens)",
    "Table S6": "Within-Type 4 heterogeneity - HUA+ vs HUA- comparison",
}
descriptions = {
    "Table S1": ("Counts of enterotype-independent core features per omics layer "
                 "(taxa, pathway, EC, metabolome), with the per-layer feature "
                 "lists for E1-only, E4-only and shared core sets."),
    "Table S2": ("Bootstrapped causal mediation estimates (ACME, ADE, total effect, "
                 "proportion mediated) for each candidate mediator in the Type 1 and "
                 "Type 4 vs reference models; BH-adjusted ACME p-values included."),
    "Table S3": ("Sensitivity analysis re-fitting the mediation models with a binary "
                 "hyperuricemia outcome instead of continuous serum urate."),
    "Table S4": ("Sensitivity analysis splitting the reference group, comparing Type 1 "
                 "and Type 4 against Type 2 (E2) and Type 3 (E3) separately."),
    "Table S5": ("Sequential-ignorability sensitivity analysis (medsens): sensitivity "
                 "parameter rho at which the mediated effect changes sign."),
    "Table S6": ("Within-Type 4 comparison of candidate mediators between HUA-positive "
                 "and HUA-negative individuals (Wilcoxon test, Cliff's delta)."),
}

sheets = {}  # sheet name -> DataFrame

# --- Table S1 -------------------------------------------------------------
df = read_csv(S1_summary)
if df is None:
    missing.append(("Table S1 summary", S1_summary))
    sheets["Table S1"] = pd.DataFrame({"NOTE": ["source CSV missing: " + S1_summary]})
elif df.empty:
    empty.append(("Table S1 summary", S1_summary))
    sheets["Table S1"] = df
else:
    blocks = []
    hdr = pd.DataFrame({df.columns[0]: ["A. Core feature counts by enterotype"]})
    blocks.append((hdr, df))
    # append per-layer lists
    layer_parts = []
    for layer, grp, p in S1_layer_files:
        ldf = read_csv(p)
        if ldf is None:
            continue
        if ldf.empty:
            empty.append((f"S1 {layer}_{grp}", p))
        section = pd.DataFrame(
            {"feature": [f"--- {layer.upper()} core: {grp} (n={len(ldf)}) ---"]}
        )
        # reindex layer df to a common schema
        layer_parts.append(section)
        layer_parts.append(ldf)
    # Build S1 as the summary table; per-layer lists go to separate companion sheets
    sheets["Table S1"] = df
    # store layer lists for separate sheets
    s1_layers = []
    for layer, grp, p in S1_layer_files:
        ldf = read_csv(p)
        if ldf is None or ldf.empty:
            continue
        s1_layers.append((layer, grp, ldf))

# --- Table S2 -------------------------------------------------------------
df = read_csv(S2_main)
if df is None:
    missing.append(("Table S2", S2_main))
    sheets["Table S2"] = pd.DataFrame({"NOTE": ["source CSV missing: " + S2_main]})
else:
    if df.empty:
        empty.append(("Table S2", S2_main))
    sheets["Table S2"] = df

# --- Table S3 -------------------------------------------------------------
df = read_csv(S3_binary)
if df is None:
    missing.append(("Table S3", S3_binary))
    sheets["Table S3"] = pd.DataFrame({"NOTE": ["source CSV missing: " + S3_binary]})
else:
    if df.empty:
        empty.append(("Table S3", S3_binary))
    sheets["Table S3"] = df

# --- Table S4 -------------------------------------------------------------
df = read_csv(S4_refgroup)
if df is None:
    missing.append(("Table S4", S4_refgroup))
    sheets["Table S4"] = pd.DataFrame({"NOTE": ["source CSV missing: " + S4_refgroup]})
else:
    if df.empty:
        empty.append(("Table S4", S4_refgroup))
    sheets["Table S4"] = df

# --- Table S5 -------------------------------------------------------------
df = read_csv(S5_medsens)
if df is None:
    missing.append(("Table S5", S5_medsens))
    sheets["Table S5"] = pd.DataFrame({"NOTE": ["source CSV missing: " + S5_medsens]})
else:
    if df.empty:
        empty.append(("Table S5", S5_medsens))
    sheets["Table S5"] = df

# --- Table S6 -------------------------------------------------------------
df = read_csv(S6_type4)
if df is None:
    missing.append(("Table S6", S6_type4))
    sheets["Table S6"] = pd.DataFrame({"NOTE": ["source CSV missing: " + S6_type4]})
else:
    if df.empty:
        empty.append(("Table S6", S6_type4))
    sheets["Table S6"] = df

# --- Contents sheet -------------------------------------------------------
contents = pd.DataFrame(
    [
        {"Table": t, "Title": titles[t], "Description": descriptions[t]}
        for t in ["Table S1", "Table S2", "Table S3", "Table S4",
                  "Table S5", "Table S6"]
    ]
)

# ---------------------------------------------------------------------------
# Write workbook
# ---------------------------------------------------------------------------
os.makedirs(OUT_DIR, exist_ok=True)
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
    contents.to_excel(xw, sheet_name="Contents", index=False)
    order = ["Table S1", "Table S2", "Table S3", "Table S4",
             "Table S5", "Table S6"]
    for name in order:
        sheets[name].to_excel(xw, sheet_name=name, index=False)
    # per-layer S1 companion sheets
    if 's1_layers' in dir() or True:
        try:
            layers = s1_layers
        except NameError:
            layers = []
        for layer, grp, ldf in layers:
            sn = f"S1_{layer}_{grp}"[:31]
            ldf.to_excel(xw, sheet_name=sn, index=False)

    # autosize columns lightly
    wb = xw.book
    for ws in wb.worksheets:
        for col in ws.columns:
            length = 10
            letter = col[0].column_letter
            for cell in col:
                v = cell.value
                if v is not None:
                    length = max(length, min(len(str(v)) + 2, 80))
            ws.column_dimensions[letter].width = length

print("WROTE:", OUT_XLSX)
print("SHEETS:", ", ".join(["Contents"] + order +
      [f"S1_{l}_{g}" for l, g, _ in (s1_layers if 's1_layers' in dir() else [])]))
print("MISSING:", missing if missing else "none")
print("EMPTY:", empty if empty else "none")
for name in order:
    print(f"  {name}: rows={len(sheets[name])} cols={len(sheets[name].columns)}")
