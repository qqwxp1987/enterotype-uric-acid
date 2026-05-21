# ==============================================================================
# Figure 3E: Shared Core ECs — Diverging Dot Plot
#
# Shows all 36 shared core enzymes with E1 and E4 mean LFC side by side.
# Classified as Convergent (same direction) or Divergent (opposite direction).
# Uses paired dots connected by lines to highlight direction agreement/conflict.
# ==============================================================================

library(tidyverse)
library(here)

# --- 1. Load data ---
ec_dir <- read.csv(here("results", "tables", "section_3_3", "ec_shared_direction.csv"))
cat(sprintf(
    "Shared core ECs: %d (%d convergent, %d divergent)\n",
    nrow(ec_dir),
    sum(ec_dir$direction_type == "Convergent"),
    sum(ec_dir$direction_type == "Divergent")
))

# --- 2. Add enzyme names from KEGG (authoritative source) ---
ec_kegg <- read.csv(here("results", "tables", "section_3_3", "ec_kegg_names.csv"))

# Only these 3 ECs are explicitly numbered in manuscript Section 3.3 for Convergence -> will be bold
manuscript_ecs <- c(
    "2.7.6.3", "2.6.1.1", "2.1.1.13"
)

ec_dir <- ec_dir %>%
    left_join(ec_kegg, by = c("feature" = "ec")) %>%
    mutate(
        # Clean up transferred / NA names
        name = case_when(
            is.na(name) ~ NA_character_,
            str_detect(name, "^Transferred") ~ NA_character_,
            TRUE ~ name
        ),
        # Build label: "EC#  enzyme name" for all with KEGG name, else just EC#
        ec_label = ifelse(!is.na(name),
            paste0(feature, "  ", name),
            feature
        ),
        is_named = feature %in% manuscript_ecs
    )

# --- 3. Order: divergent first (by |E4_mean_lfc - E1_mean_lfc|), then convergent ---
ec_dir <- ec_dir %>%
    mutate(
        delta = abs(E4_mean_lfc - E1_mean_lfc),
        direction_type = factor(direction_type, levels = c("Divergent", "Convergent"))
    ) %>%
    arrange(direction_type, delta) %>%
    mutate(ec_label = factor(ec_label, levels = ec_label))

# --- 4. Colors ---
col_e1 <- "#E64B35" # Type 1
col_e4 <- "#3C5488" # Type 4

# --- 5. Reshape for paired dot plot ---
ec_long <- ec_dir %>%
    select(ec_label, E1_mean_lfc, E4_mean_lfc, direction_type, is_named) %>%
    pivot_longer(
        cols = c(E1_mean_lfc, E4_mean_lfc),
        names_to = "enterotype", values_to = "mean_lfc"
    ) %>%
    mutate(enterotype = ifelse(enterotype == "E1_mean_lfc", "Type 1", "Type 4"))

# --- 6. Plot ---
p <- ggplot() +

    # Zero line
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey50") +

    # Connecting segments between E1 and E4 dots
    geom_segment(
        data = ec_dir,
        aes(
            x = E1_mean_lfc, xend = E4_mean_lfc,
            y = ec_label, yend = ec_label
        ),
        linewidth = 0.4,
        color = ifelse(ec_dir$direction_type == "Divergent", "#E8D5B7", "#D5E8D4"),
        alpha = 0.8
    ) +

    # Dots
    geom_point(
        data = ec_long,
        aes(x = mean_lfc, y = ec_label, color = enterotype),
        size = 2.2, alpha = 0.85
    ) +

    # Colors
    scale_color_manual(
        values = c("Type 1" = col_e1, "Type 4" = col_e4),
        name = NULL
    ) +

    # Facet by direction type
    facet_grid(direction_type ~ .,
        scales = "free_y", space = "free_y",
        switch = "y"
    ) +
    labs(
        x = expression("Mean log"[2] * " fold-change"),
        y = NULL
    ) +

    # Theme matching other panels
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
        axis.text.y = element_text(
            size = 7.5,
            color = "grey20",
            face = ifelse(levels(ec_dir$ec_label) %in%
                ec_dir$ec_label[ec_dir$is_named], "bold", "plain")
        ),
        axis.text.x = element_text(size = 9, color = "grey35"),
        axis.title.x = element_text(size = 10, margin = margin(t = 6)),
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 10, face = "bold", angle = 0),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.background = element_rect(fill = "#FAFBFC", color = NA),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        plot.margin = margin(t = 8, r = 12, b = 8, l = 5)
    )

# --- 7. Save ---
out_dir <- here("results", "figures", "panels")

ggsave(file.path(out_dir, "Figure4E_Shared_EC_Dumbbell.pdf"), p,
    width = 8, height = 9, units = "in"
)
ggsave(file.path(out_dir, "Figure4E_Shared_EC_Dumbbell.png"), p,
    width = 8, height = 9, units = "in", dpi = 300
)

cat("\nPanel E (Dumbbell Plot) saved successfully as Figure4E_Shared_EC_Dumbbell.\n")
