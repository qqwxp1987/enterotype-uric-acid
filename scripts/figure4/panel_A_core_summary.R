# ==============================================================================
# Figure 3A: Core Feature Summary — Butterfly / Mirror Bar Chart
#
# Type 1 core features extend LEFT, Type 4 extend RIGHT.
# Shared counts annotated at center axis.
# Each row = one data layer (Taxa, Pathway, EC, Metabolome).
#
# Colors match global enterotype palette:
#   Type 1 = #E64B35, Type 4 = #3C5488
# ==============================================================================

library(tidyverse)
library(here)

# --- 1. Load data ---
summary_df <- read.csv(here("results", "tables", "section_3_3", "independent_core_summary.csv"))
cat("Data loaded:\n")
print(summary_df)

# --- 2. Prepare butterfly layout ---
plot_df <- summary_df %>%
    mutate(
        data_type = factor(data_type,
            levels = rev(c("Taxa", "Pathway", "EC", "Metabolome"))
        ),
        E1_plot = -E1_only,
        E4_plot = E4_only
    )

# --- 3. Global enterotype color palette (matching Figure 1 & 2) ---
col_e1 <- "#E64B35" # Type 1: warm red (ggsci nrc)
col_e4 <- "#3C5488" # Type 4: navy blue (ggsci nrc)
col_shared <- "#FFB347" # amber for shared badges
col_bg <- "#FAFBFC" # subtle background

# --- 4. Build the plot ---
max_val <- max(c(summary_df$E1_only, summary_df$E4_only)) * 1.15

p <- ggplot(plot_df, aes(y = data_type)) +

    # ----- Background panels for left/right -----
    annotate("rect",
        xmin = -max_val, xmax = 0,
        ymin = 0.4, ymax = 4.6,
        fill = col_e1, alpha = 0.03
    ) +
    annotate("rect",
        xmin = 0, xmax = max_val,
        ymin = 0.4, ymax = 4.6,
        fill = col_e4, alpha = 0.03
    ) +

    # ----- Type 1 bars (LEFT) -----
    geom_col(aes(x = E1_plot), fill = col_e1, width = 0.55, color = NA) +

    # ----- Type 4 bars (RIGHT) -----
    geom_col(aes(x = E4_plot), fill = col_e4, width = 0.55, color = NA) +

    # ----- Center axis line -----
    geom_vline(xintercept = 0, linewidth = 0.6, color = "grey40") +

    # ----- Count labels on bars -----
    # Type 1 counts (left, outside bar)
    geom_text(aes(x = E1_plot, label = E1_only),
        hjust = 1.15, size = 3.8, fontface = "bold", color = col_e1
    ) +

    # Type 4 counts (right, outside bar — offset small values to avoid badge overlap)
    geom_text(
        aes(
            x = ifelse(E4_only < 10, pmax(E4_plot, 18), E4_plot),
            label = E4_only
        ),
        hjust = -0.15, size = 3.8, fontface = "bold", color = col_e4
    ) +

    # ----- Shared count badges (center) -----
    geom_point(aes(x = 0),
        size = ifelse(plot_df$Shared > 0, 7, 5),
        shape = 21, fill = ifelse(plot_df$Shared > 0, col_shared, "grey85"),
        color = ifelse(plot_df$Shared > 0, "grey30", "grey70"),
        stroke = 0.5
    ) +
    geom_text(aes(x = 0, label = Shared),
        size = 3.0, fontface = "bold",
        color = ifelse(plot_df$Shared > 0, "grey20", "grey50")
    ) +

    # ----- Axis labels -----
    scale_x_continuous(
        limits = c(-max_val, max_val),
        breaks = seq(-400, 400, by = 100),
        labels = function(x) abs(x),
        expand = expansion(mult = c(0.02, 0.02))
    ) +

    # ----- Header annotations -----
    annotate("text",
        x = -max_val * 0.5, y = 4.85,
        label = "\u2190 Type 1 only", fontface = "bold",
        size = 4, color = col_e1, hjust = 0.5
    ) +
    annotate("text",
        x = max_val * 0.5, y = 4.85,
        label = "Type 4 only \u2192", fontface = "bold",
        size = 4, color = col_e4, hjust = 0.5
    ) +
    annotate("text",
        x = 0, y = 4.85,
        label = "Shared", fontface = "bold.italic",
        size = 3.2, color = "grey40", hjust = 0.5
    ) +
    labs(x = "Number of core features", y = NULL) +

    # ----- Theme -----
    theme_minimal(base_size = 12, base_family = "") +
    theme(
        plot.background = element_rect(fill = col_bg, color = NA),
        panel.background = element_rect(fill = col_bg, color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(
            color = "grey88", linewidth = 0.25,
            linetype = "dotted"
        ),
        axis.text.y = element_text(
            size = 11, face = "bold", color = "grey25",
            margin = margin(r = 5)
        ),
        axis.text.x = element_text(size = 9, color = "grey45"),
        axis.title.x = element_text(size = 10, color = "grey40", margin = margin(t = 8)),
        plot.margin = margin(t = 25, r = 20, b = 10, l = 10),
        plot.title = element_blank()
    ) +
    coord_cartesian(clip = "off", ylim = c(0.5, 4.5))

# --- 5. Save ---
out_dir <- here("results", "figures", "panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "Figure4A_Core_Summary.pdf"), p,
    width = 7, height = 3.2, units = "in"
)
ggsave(file.path(out_dir, "Figure4A_Core_Summary.png"), p,
    width = 7, height = 3.2, units = "in", dpi = 300
)

cat("\nFigure 4A saved to:", out_dir, "\n")
