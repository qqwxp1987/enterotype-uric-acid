# Validation of Core EC Enzymes in Section 3.3

This document provides the official KEGG annotations for the 36 shared core Enzyme Commission (EC) numbers between Type 1 and Type 4 enterotypes, acting as a reference baseline to verify the descriptive logic in Section 3.3 of the manuscript.

## 1. Verified Core Enzymes Explicitly Named in Section 3.3

A direct comparison against the KEGG Database confirms that all EC enzymes explicitly named in the text are **accurately annotated**:

| EC Number | Named in Manuscript | Official KEGG Name | KEGG Class | Status |
| :--- | :--- | :--- | :--- | :--- |
| **EC 1.4.7.1** | glutamate synthase | glutamate synthase (ferredoxin) | Oxidoreductases | ✅ Accurate |
| **EC 3.4.11.5** | prolyl aminopeptidase | prolyl aminopeptidase | Hydrolases | ✅ Accurate |
| **EC 2.6.1.11** | acetylornithine transaminase | acetylornithine transaminase | Transferases | ✅ Accurate |
| **EC 6.1.1.5** | isoleucine-tRNA ligase | isoleucine---tRNA ligase | Ligases | ✅ Accurate |
| **EC 6.1.1.9** | valine-tRNA ligase | valine---tRNA ligase | Ligases | ✅ Accurate |
| **EC 2.4.2.45** | decaprenyl-phosphate phosphoribosyltransferase | decaprenyl-phosphate phosphoribosyltransferase | Transferases | ✅ Accurate |
| **EC 2.4.2.9** | uracil phosphoribosyltransferase | uracil phosphoribosyltransferase | Transferases | ✅ Accurate |

## 2. Validation of the 14 Convergent Shared ECs

The manuscript states: *"The remaining 14 shared enzymes that showed convergent changes... were predominantly involved in basic carbohydrate processing and one-carbon metabolism."* 

Below is the full list of the 14 convergent enzymes, which validates this claim:

**Carbohydrate Processing (7 enzymes):**
*   **EC 3.2.1.22:** alpha-galactosidase (Hydrolases)
*   **EC 3.2.1.23:** beta-galactosidase (Hydrolases)
*   **EC 3.2.1.55:** non-reducing end alpha-L-arabinofuranosidase (Hydrolases)
*   **EC 4.2.1.8:** mannonate dehydratase (Lyases)
*   **EC 5.1.3.3:** aldose 1-epimerase (Isomerases)
*   **EC 5.4.99.62:** D-ribose pyranase (Isomerases)
*   **EC 5.1.3.n2:** (Unassigned novel epimerase, KEGG unmapped, functionally clustering with carbohydrate epimerization)

**One-Carbon / Folate Metabolism (2 enzymes):**
*   **EC 2.1.1.13:** methionine synthase (Transferases, 1-carbon transfer)
*   **EC 2.7.6.3:** 2-amino-4-hydroxy-6-hydroxymethyldihydropteridine diphosphokinase (Transferases, folate biosynthesis)

**Other Convergent Functions (5 enzymes):**
*   **EC 1.16.3.1:** ferroxidase
*   **EC 2.6.1.1:** aspartate transaminase
*   **EC 2.3.1.129:** acyl-[acyl-carrier-protein]---UDP-N-acetylglucosamine O-acyltransferase
*   **EC 3.1.3.15:** histidinol-phosphatase
*   **EC 4.1.1.37:** uroporphyrinogen decarboxylase

**Conclusion:** The claim holds true. Exactly half (7/14) of the convergent enzymes are heavily involved in plant-glycan and fundamental carbohydrate processing, while another two represent the core backbone of one-carbon/methyl-group transfer. 

## 3. Full List of Remaining Divergent Shared ECs (22 total minus those named in Section 1)
For completeness, these are the remaining divergent ECs not specifically mentioned by name in the text, representing broader pathway divergences.
*   **EC 1.16.3.1:** ferroxidase
*   **EC 1.5.1.7:** saccharopine dehydrogenase (NAD+, L-lysine-forming)
*   **EC 1.8.1.8:** protein-disulfide reductase
*   **EC 2.1.1.148:** thymidylate synthase (FAD)
*   **EC 2.1.1.37:** DNA (cytosine-5-)-methyltransferase
*   **EC 2.7.1.50:** hydroxyethylthiazole kinase
*   **EC 2.7.4.7:** phosphooxymethylpyrimidine kinase
*   **EC 2.7.9.1:** pyruvate, phosphate dikinase
*   **EC 3.2.1.4:** cellulase
*   **EC 3.4.24.55:** pitrilysin
*   **EC 3.6.1.13:** ADP-ribose diphosphatase
*   **EC 3.6.3.14:** Transferred to 7.1.2.2 (proton-translocating ATPase)
*   **EC 4.1.1.41:** Transferred to 7.2.4.3
*   **EC 4.3.3.6:** pyridoxal 5'-phosphate synthase (glutamine hydrolysing)
*   **EC 6.3.1.5:** NAD+ synthase
*   **EC 6.3.4.15:** biotin---[biotin carboxyl-carrier protein] ligase
