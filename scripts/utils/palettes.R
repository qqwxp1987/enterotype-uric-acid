# ==============================================================================
# Global palettes shared across all figures.
# Single source of truth for enterotype & genus colors.
# Source this file from any panel script: source(here("scripts", "utils", "palettes.R"))
# ==============================================================================

# 4-color enterotype palette (NPG-style; color-blind safe against each other)
ENTEROTYPE_COLORS <- c(
  "Type 1" = "#E64B35",  # red
  "Type 2" = "#4DBBD5",  # cyan
  "Type 3" = "#00A087",  # teal
  "Type 4" = "#3C5488"   # blue
)

# Genus palette used by stacked-bar panels (Figure 1C and Supp S1).
# Bacteroides / Prevotella deliberately aligned with Type 1 / Type 4 identity.
GENUS_COLORS <- c(
  "Bacteroides"                  = "#E64B35",
  "Prevotella"                   = "#3C5488",
  "Alistipes"                    = "#F39B7F",
  "Blautia"                      = "#91D1C2",
  "Eubacterium"                  = "#8491B4",
  "Faecalibacterium"             = "#7E6148",
  "Lachnospiraceae (unclassified)" = "#B09C85",
  "Parabacteroides"              = "#00A087",
  "Roseburia"                    = "#4DBBD5",
  "Ruminococcus"                 = "#DC9A6C",
  "Others"                       = "#D9D9D9"
)

# Omics / data-layer palette shared by cross-omics figures (Figure 5 & 6).
# Single source of truth so "Convergent EC" etc. stay consistent across figures.
LAYER_COLORS <- c(
  "Serum metabolite"  = "#A6611A",
  "Microbial pathway" = "#2C7FB8",
  "Convergent EC"     = "#1B9E77",
  "Taxonomy"          = "#7B5EA7"
)
