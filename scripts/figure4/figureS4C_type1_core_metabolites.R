# ==============================================================================
# Figure 3C: Type 1 Core Serum Metabolites — Lollipop Chart
#
# Unified style with Panel B (lollipop chart).
# All 24 core metabolites displayed.
# Mean log2FC = average of E1vE2 and E1vE3 comparisons.
# Metabolite class indicated by colored dots on left margin.
# Colors: Higher = #E64B35, Lower = #3C5488 (same as Panel B).
# ==============================================================================

library(tidyverse)
library(here)

# --- 1. Load data ---
metab_e1 <- read.csv(here("results", "tables", "section_3_4", "metabolome_E1_core.csv"))
cat(sprintf("E1 core metabolites: %d\n", nrow(metab_e1)))

# --- 2. Compute mean LFC & clean names ---
metab_e1 <- metab_e1 %>%
    mutate(
        mean_lfc = (log2FoldChange_E1vE2 + log2FoldChange_E1vE3) / 2,
        ion_mode = str_extract(feature, "^(POS|NEG)"),
        name_raw = str_remove(feature, "^(POS|NEG)\\d+_"),
        direction = ifelse(mean_lfc > 0, "Higher", "Lower")
    )

# Handle duplicate names: append ion mode if duplicated
dup_names <- metab_e1 %>%
    count(name_raw) %>%
    filter(n > 1) %>%
    pull(name_raw)
metab_e1 <- metab_e1 %>%
    mutate(
        name_display = ifelse(name_raw %in% dup_names,
            paste0(name_raw, " [", ion_mode, "]"),
            name_raw
        )
    )

# Assign functional group
metab_e1 <- metab_e1 %>%
    mutate(
        group = case_when(
            str_detect(name_raw, regex("cresol|indol|indox|phenylacetyl|methoxyindole|dihydroindol", ignore_case = TRUE)) ~ "AAA catabolites",
            str_detect(name_raw, regex("cholic|bile|deoxy", ignore_case = TRUE)) ~ "Bile acids",
            TRUE ~ "Other"
        )
    ) %>%
    arrange(mean_lfc) %>%
    mutate(name_display = factor(name_display, levels = name_display))

cat("Groups:\n")
print(table(metab_e1$group))

# --- 3. Colors (matching Panel B exactly) ---
col_higher <- "#E64B35"
col_lower <- "#3C5488"

# Metabolite-class colors aligned with Figure 5 (scripts/figure5/*: class_colors)
group_colors <- c(
    "AAA catabolites" = "#D2691E",
    "Bile acids"      = "#2E8B57",
    "Other"           = "#BDC3C7"
)

# --- 4. Lollipop plot (same style as Panel B) ---
p <- ggplot(metab_e1, aes(x = mean_lfc, y = name_display)) +

    # Zero reference line
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey50", linetype = "solid") +

    # Lollipop segments
    geom_segment(
        aes(
            x = 0, xend = mean_lfc, y = name_display, yend = name_display,
            color = direction
        ),
        linewidth = 0.6
    ) +

    # Lollipop dots
    geom_point(aes(color = direction, fill = direction),
        size = 2.5, shape = 21, stroke = 0.3
    ) +

    # Metabolite class indicator: small squares on right margin
    geom_point(aes(x = max(metab_e1$mean_lfc) + 0.35, fill = group),
        size = 2.2, shape = 22, stroke = 0.3, color = "grey50",
        show.legend = TRUE
    ) +

    # Direction colors
    scale_color_manual(
        values = c("Higher" = col_higher, "Lower" = col_lower),
        name = NULL
    ) +
    scale_fill_manual(
        values = c(
            "Higher" = col_higher, "Lower" = col_lower,
            "AAA catabolites" = "#D2691E", "Bile acids" = "#2E8B57", "Other" = "#BDC3C7"
        ),
        breaks = c("AAA catabolites", "Bile acids", "Other"),
        name = "Metabolite class"
    ) +
    labs(
        x = expression("Mean log"[2] * " fold-change"),
        y = NULL
    ) +

    # Theme: matching Panel B exactly
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
        axis.text.y = element_text(size = 8.5, color = "grey20"),
        axis.text.x = element_text(size = 9, color = "grey35"),
        axis.title.x = element_text(size = 10, margin = margin(t = 6)),
        legend.position = "top",
        legend.box = "vertical",
        legend.justification = "center",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        plot.background = element_rect(fill = "#FAFBFC", color = NA),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        plot.margin = margin(t = 8, r = 15, b = 8, l = 5)
    ) +
    coord_cartesian(clip = "off") +
    guides(
        color = guide_legend(order = 1),
        fill  = guide_legend(order = 2, override.aes = list(shape = 22, size = 3, color = "grey50"))
    )

# --- 5. Save ---
out_dir <- here("results", "figures", "panels")

ggsave(file.path(out_dir, "Figure4C_Type1_Core_Metabolites.pdf"), p,
    width = 7.5, height = 6, units = "in"
)
ggsave(file.path(out_dir, "Figure4C_Type1_Core_Metabolites.png"), p,
    width = 7.5, height = 6, units = "in", dpi = 300
)

cat("\nFigure 4C saved to:", out_dir, "\n")
