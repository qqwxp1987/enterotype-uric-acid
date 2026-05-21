# ==============================================================================
# Figure 3E: Species Attribution of Key Pathways (Stacked Bar Plot)
#
# Generates a 100% stacked bar plot comparing species contributions 
# to selected key pathways across Enterotype 1 and Enterotype 4.
# Replicates the style of the user-provided screenshot.
# ==============================================================================

library(tidyverse)
library(here)

# --- 1. Load Data ---
res <- readRDS(here("results", "tables", "section_3_3", "attribution", "pathway_contribution_data.rds"))
df <- res$contrib_df

# Define key pathways based on the text (Section 3.3)
key_pathways <- c(
  "PWY-5695"             = "Urate biosyn. / IMP degrad.",
  "PWY-7219"             = "Adenosine ribonucleotides de novo",
  "PWY-6608"             = "Guanosine nucleotides degrad. III",
  "PEPTIDOGLYCANSYN-PWY" = "Peptidoglycan biosynthesis I",
  "PWY-7234"             = "IMP biosynthesis III",
  "ARGSYN-PWY"           = "L-Arg biosynthesis I"
)

# Filter for Type 1 & 4 and the key pathways
plot_df <- df %>%
  filter(Enterotype %in% c("1", "4"), pathway_id %in% names(key_pathways)) %>%
  mutate(
    Pathway = recode(pathway_id, !!!key_pathways),
    Enterotype = paste0("Type ", Enterotype)
  )

# Calculate the mean contribution percentage per species per pathway per enterotype
mean_pct <- plot_df %>%
  group_by(Enterotype, Pathway, species) %>%
  summarise(mean_pct = mean(contribution_pct, na.rm=TRUE), .groups = "drop")

# We want to identify top species across ALL plotted pathways to keep colors consistent
# Also, group unclassified separately
mean_pct <- mean_pct %>%
  mutate(
    # Unclassified often has 'unclassified' in name
    is_unclassified = str_detect(tolower(species), "unclassified")
  )

# --- 2. Custom Species Classification ---
# Identify which specific Bacteroides and Prevotella species exist in this filtered subset
# so we can build dynamic gradient palettes for them.
available_spp <- unique(mean_pct$species[!mean_pct$is_unclassified])

bacteroides_spp <- sort(available_spp[str_detect(available_spp, "^Bacteroides ")])
prevotella_spp <- sort(available_spp[str_detect(available_spp, "^Prevotella ")])

target_spp_3 <- c("Faecalibacterium prausnitzii", "Bifidobacterium pseudocatenulatum", "Alistipes onderdonkii")
actual_target_spp_3 <- target_spp_3[target_spp_3 %in% available_spp]

# The universe of "named" species we care about
focus_species <- c(bacteroides_spp, prevotella_spp, actual_target_spp_3)

# Recode species
agg_df <- mean_pct %>%
  mutate(
    Species = case_when(
      species %in% focus_species ~ species,
      TRUE ~ "Other/Unclassified"
    )
  ) %>%
  group_by(Enterotype, Pathway, Species) %>%
  summarise(mean_pct = sum(mean_pct), .groups = "drop")

# Normalize to exactly 100% per Enterotype/Pathway
agg_df <- agg_df %>%
  group_by(Enterotype, Pathway) %>%
  mutate(pct_normalized = mean_pct / sum(mean_pct)) %>%
  ungroup()

# Ordering for the plot
# Pathways: Order as defined in key_pathways
pathway_order <- rev(unname(key_pathways))
pathway_order <- pathway_order[pathway_order %in% unique(agg_df$Pathway)]
agg_df$Pathway <- factor(agg_df$Pathway, levels = pathway_order)

# Species ordering: 
# 1. P. copri strains at the top (so they plot together)
# 2. Other Prevotella
# 3. Bacteroides
# 4. Target 3
# 5. Other/Unclassified
p_copri_strains <- prevotella_spp[str_detect(prevotella_spp, "Prevotella copri")]
other_prevotella <- setdiff(prevotella_spp, p_copri_strains)

spp_levels <- c(
  sort(p_copri_strains),      # Darkest Oranges/Reds
  sort(other_prevotella),     # Lighter Oranges/Reds
  sort(bacteroides_spp),      # Blues
  sort(actual_target_spp_3),  # Key targets
  "Other/Unclassified"
)
agg_df$Species <- factor(agg_df$Species, levels = spp_levels)

# --- 3. Color Palette Assembly ---
library(RColorBrewer)

# Generates gradients. We reverse them if needed to ensure the first elements (like P. copri) get the darkest colors.
generate_gradient <- function(n, base_palette, reverse = FALSE) {
  if (n == 0) return(character(0))
  if (n == 1) return(brewer.pal(3, base_palette)[3]) # Just pick a decent mid-dark color
  # colorRampPalette generates from light to dark if the base palette is sequential
  pal <- colorRampPalette(brewer.pal(9, base_palette))(n)
  if (reverse) {
    return(rev(pal))
  } else {
    return(pal)
  }
}

# 1. Prevotella (Oranges/Reds) - The levels matrix ensures P. copri is first. We want darkest colors for them.
# So we generate a gradient and Reverse it.
n_prev <- length(prevotella_spp)
color_prevotella <- generate_gradient(n_prev, "OrRd", reverse = TRUE)
names(color_prevotella) <- c(sort(p_copri_strains), sort(other_prevotella))

# 2. Bacteroides (Blues) - Gradient
n_bact <- length(bacteroides_spp)
color_bacteroides <- generate_gradient(n_bact, "Blues", reverse = TRUE)
names(color_bacteroides) <- sort(bacteroides_spp)

# 3. Target 3 (Distinct Solids)
# Assign fixed colors for consistency, but only include if they exist in this cut
fixed_targets <- c(
  "Faecalibacterium prausnitzii" = "#2ca02c",      # Distinct Green
  "Bifidobacterium pseudocatenulatum" = "#9467bd", # Distinct Purple
  "Alistipes onderdonkii" = "#8c564b"              # Distinct Brown/Gold
)
color_targets <- fixed_targets[actual_target_spp_3]

custom_colors <- c(
  color_prevotella,
  color_bacteroides,
  color_targets,
  "Other/Unclassified" = "#D3D3D3" # Light Grey
)

# --- 3. Plotting the Stacked Bar Chart ---
p <- ggplot(agg_df, aes(x = pct_normalized, y = Pathway, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75, color = "white", linewidth = 0.1) +
  facet_wrap(~ Enterotype, ncol = 2) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 22)) +
  labs(
    x = "Mean species contribution to pathway abundance",
    y = NULL,
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    axis.text.y = element_text(size = 10, color = "black", lineheight = 0.8),
    axis.text.x = element_text(size = 10, color = "grey30"),
    axis.title.x = element_text(size = 11, margin = margin(t = 10)),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8, face = "italic"),
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )

# --- 4. Save The Plot ---
out_dir <- here("results", "figures", "panels")

ggsave(file.path(out_dir, "Figure3E_Species_Attribution_Bar.png"), plot = p, width = 12, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure3E_Species_Attribution_Bar.pdf"), plot = p, width = 12, height = 6, device = cairo_pdf)

cat("Panel E (Stacked Bar Chart) saved successfully as Figure3E_Species_Attribution_Bar.\n")
