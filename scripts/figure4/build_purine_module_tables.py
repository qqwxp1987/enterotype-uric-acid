#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build final Tables S3, S9 and S10 from the verified purine-module outputs.

This script repackages frozen statistical results; it does not refit models.
Run after purine_module_analysis.R from any working directory.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CURATED = ROOT / "data" / "curated" / "purine"
RESULTS = ROOT / "results" / "tables" / "section_3_4" / "purine_modules"
OUT = ROOT / "results" / "tables" / "final"
OUT.mkdir(parents=True, exist_ok=True)

# Table S3: composition-based curation of the 34 purine-related pathways.
reg = pd.read_csv(CURATED / "frozen_inclusion_register.csv", encoding="utf-8-sig")
inc = reg[reg["final_action"] != "exclude"].copy()
assert len(reg) == 474 and len(inc) == 34
role_map = {
    "core": "module-defining member",
    "redundant_in_sum": "descriptive annotation (reactions fully contained in another member)",
    "hybrid_excluded": "descriptive annotation (spans modules; excluded from module membership)",
}
module_map = {"de_novo": "De novo biosynthesis", "degradation": "Degradation", "salvage": "Salvage"}
t1 = pd.DataFrame({
    "Pathway ID": inc["pathway_id"],
    "Pathway name": inc["feature"].str.split(": ", n=1).str[1],
    "Module": inc["final_module"].map(module_map),
    "Curation role": inc["final_sum_role"].map(role_map),
    "Number of reactions": inc["n_reactions"],
})
t1.to_csv(OUT / "Table_S3_purine_pathway_curation.csv", index=False, encoding="utf-8-sig")
note = (
    "Note. PWY-6595 is annotated in MetaCyc as 'superpathway of guanosine nucleotides "
    "degradation (plants)'; the '(plants)' suffix denotes the curation source of the "
    "pathway definition, not its ecological occurrence - its constituent reactions are "
    "generic guanosine nucleotide degradation steps found in gut bacteria."
)
(OUT / "Table_S3_purine_pathway_curation_footnote.txt").write_text(note, encoding="utf-8")

# Table S9: group medians, overall GLS tests and six Holm-adjusted contrasts.
gls = pd.read_csv(RESULTS / "03_q2_gls_contrasts.csv")
med = pd.read_csv(RESULTS / "verification" / "v07_direction_check_by_median.csv")
mod_names = {
    "E_de_novo_share": "De novo biosynthesis",
    "E_degradation_share": "Degradation",
    "E_salvage_share": "Salvage",
}
rows = []
for score, name in mod_names.items():
    m = med[med["score"] == score].iloc[0]
    for enterotype in range(1, 5):
        rows.append({"Module": name, "Statistic": f"Median share (%), Type {enterotype}", "Value": m[f"T{enterotype}"]})
    rows.append({
        "Module": name,
        "Statistic": "Overall enterotype effect, P (GLS, marginal ANOVA)",
        "Value": m["overall_P"],
    })
    for _, r in gls[gls["score"] == score].iterrows():
        rows.append({
            "Module": name,
            "Statistic": f"Pairwise contrast {r['contrast']} (log10 scale)",
            "Value": f"estimate={r['estimate']:.4f}, SE={r['SE']:.4f}, Holm P={r['p_holm']:.3g}",
        })
t2 = pd.DataFrame(rows)
assert len(t2) == 33
t2.to_csv(OUT / "Table_S9_purine_module_statistics.csv", index=False, encoding="utf-8-sig")

# Table S10: partial Spearman correlations, with BH adjustment across three modules.
pc = pd.read_csv(RESULTS / "04_q3_module_UA_pcor.csv")
pc = pc[pc["score"].isin(mod_names)].copy()
assert len(pc) == 3 and set(pc["n"]) == {305}
pvals = pc["p"].to_numpy()
order = np.argsort(pvals)
adjusted_sorted = np.minimum.accumulate(
    (pvals[order] * len(pvals) / np.arange(1, len(pvals) + 1))[::-1]
)[::-1]
padj = np.empty(len(pvals))
padj[order] = np.minimum(adjusted_sorted, 1.0)
pc["BH-adjusted P (across the 3 module scores)"] = padj
t3 = pd.DataFrame({
    "Module": pc["score"].map(mod_names),
    "Partial Spearman rho (vs serum UA, adjusted for sex and age)": pc["rho"].round(3),
    "P": pc["p"].map("{:.3g}".format),
    "BH-adjusted P (across the 3 module scores)": pc["BH-adjusted P (across the 3 module scores)"].map("{:.3g}".format),
    "n": pc["n"].astype(int),
})
de_novo = t3[t3["Module"] == "De novo biosynthesis"].iloc[0]
assert abs(float(de_novo["BH-adjusted P (across the 3 module scores)"]) - 0.0365) < 1e-12
t3.to_csv(OUT / "Table_S10_purine_module_UA.csv", index=False, encoding="utf-8-sig")

sys.stdout.reconfigure(encoding="utf-8")
print(f"Wrote S3 ({len(t1)} rows), S9 ({len(t2)} rows), and S10 ({len(t3)} rows) to {OUT}")
print(t3.to_string(index=False))
