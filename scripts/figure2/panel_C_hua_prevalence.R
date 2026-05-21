# ==============================================================================
# Figure 2C: Enterotype Distribution Across Sex × HUA Subgroups
# Stacked bar with alluvial flow connections
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(ggalluvial)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")

meta$Enterotype <- factor(
  meta$Dirichlettype_Genus,
  levels = c("4", "3", "2", "1"),
  labels = c("Type 4", "Type 3", "Type 2", "Type 1")
)
meta$Gender <- factor(meta$Gender, levels = c("Male", "Female"))
meta$Hyperuricemia <- factor(meta$Hyperuricemia, levels = c("Yes", "No"))

# --- 2. Color Palette ---
enterotype_colors <- c(
  "Type 1" = "#E64B35",
  "Type 2" = "#4DBBD5",
  "Type 3" = "#00A087",
  "Type 4" = "#3C5488"
)

# --- 3. Compute Proportions ---
# Create Sex × HUA groups
meta$Group <- interaction(meta$Gender, meta$Hyperuricemia, sep = " ")
meta$Group <- factor(meta$Group, levels = c(
  "Male Yes", "Male No", "Female Yes", "Female No"
))

# Count table
count_data <- meta %>%
  count(Group, Enterotype) %>%
  group_by(Group) %>%
  mutate(
    pct = n / sum(n) * 100,
    total_n = sum(n)
  ) %>%
  ungroup()

# Display labels for x-axis
group_labels <- c(
  "Male Yes"   = "Male\nHUA+",
  "Male No"    = "Male\nHUA\u2212",
  "Female Yes" = "Female\nHUA+",
  "Female No"  = "Female\nHUA\u2212"
)

count_data$Group_label <- group_labels[as.character(count_data$Group)]
count_data$Group_label <- factor(count_data$Group_label,
  levels = group_labels
)

cat("=== Enterotype distribution by Sex × HUA ===\n")
print(count_data %>% select(Group, Enterotype, n, pct) %>% as.data.frame())

# --- 4. Statistical Tests ---
cat("\n=== Statistical Tests ===\n")

# Overall: Enterotype × HUA (ignoring sex)
hua_tab_all <- table(meta$Dirichlettype_Genus, meta$Hyperuricemia)

# Male subgroup
meta_male <- meta %>% filter(Gender == "Male")
male_tab <- table(meta_male$Enterotype, meta_male$Hyperuricemia)
set.seed(42)
male_fisher <- fisher.test(male_tab, simulate.p.value = TRUE, B = 10000)
cat(sprintf(
  "Male:   Enterotype × HUA  Fisher's P = %.4f (N=%d)\n",
  male_fisher$p.value, nrow(meta_male)
))

# Female subgroup
meta_female <- meta %>% filter(Gender == "Female")
female_tab <- table(meta_female$Enterotype, meta_female$Hyperuricemia)
set.seed(42)
female_fisher <- fisher.test(female_tab, simulate.p.value = TRUE, B = 10000)
cat(sprintf(
  "Female: Enterotype × HUA  Fisher's P = %.4f (N=%d)\n",
  female_fisher$p.value, nrow(meta_female)
))

# Overall chi-squared
overall_tab <- table(meta$Enterotype, meta$Hyperuricemia)
overall_chi <- chisq.test(overall_tab)
cat(sprintf(
  "Overall: χ² = %.2f, P = %.2e\n",
  overall_chi$statistic, overall_chi$p.value
))

# CMH test (Enterotype × HUA, stratified by Sex)
cmh <- mantelhaen.test(
  table(meta$Enterotype, meta$Hyperuricemia, meta$Gender)
)
cat(sprintf(
  "CMH test (stratified by sex): χ² = %.2f, P = %.4f\n",
  cmh$statistic, cmh$p.value
))

# --- 6. Plot ---

# For alluvial flow, use a lode-form data
alluvial_data <- count_data %>%
  select(Group_label, Enterotype, n) %>%
  rename(x = Group_label, stratum = Enterotype, freq = n)

p <- ggplot(
  count_data,
  aes(
    x = Group_label, y = pct, fill = Enterotype,
    stratum = Enterotype, alluvium = Enterotype
  )
) +
  # Alluvial flow between bars
  geom_flow(
    aes(y = pct),
    alpha = 0.25, width = 0.45,
    curve_type = "sigmoid"
  ) +
  # Stacked bar
  geom_stratum(
    aes(y = pct),
    width = 0.45, color = "white", linewidth = 0.4
  ) +
  # Percentage labels
  geom_text(
    stat = "stratum",
    aes(label = ifelse(pct >= 5, sprintf("%.1f%%", pct), "")),
    size = 3, fontface = "bold", color = "black"
  ) +
  # Statistical annotations
  annotate("text",
    x = 1.5, y = 105,
    label = sprintf("italic(P) == %.3f", male_fisher$p.value),
    parse = TRUE, size = 3.2, color = "grey30"
  ) +
  annotate("segment",
    x = 1, xend = 2, y = 103, yend = 103,
    linewidth = 0.3, color = "grey50"
  ) +
  annotate("text",
    x = 3.5, y = 105,
    label = sprintf("italic(P) == %.3f", female_fisher$p.value),
    parse = TRUE, size = 3.2, color = "grey30"
  ) +
  annotate("segment",
    x = 3, xend = 4, y = 103, yend = 103,
    linewidth = 0.3, color = "grey50"
  ) +
  # N labels
  annotate("text",
    x = 1:4,
    y = rep(-4, 4),
    label = count_data %>%
      group_by(Group_label) %>%
      summarise(N = first(total_n), .groups = "drop") %>%
      mutate(lab = paste0("n=", N)) %>%
      pull(lab),
    size = 3, color = "grey40"
  ) +
  # CMH test annotation (Enterotype × HUA, stratified by sex)
  annotate("text",
    x = 2.5, y = 116,
    label = sprintf(
      "CMH (sex-stratified): italic(P) == %.3f",
      cmh$p.value
    ),
    parse = TRUE, size = 3.2, fontface = "bold", color = "grey20"
  ) +
  # Scales
  scale_fill_manual(values = enterotype_colors, name = "Enterotype") +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x),
    expand = expansion(mult = c(0.08, 0.16))
  ) +
  labs(x = NULL, y = "Percentage (%)") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 10, color = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10),
    plot.margin = margin(22, 10, 15, 10)
  )

# --- 7. Save ---
out_dir <- here("results", "figures", "panels")

ggsave(file.path(out_dir, "Figure2C_HUA_Prevalence.pdf"),
  p,
  width = 6, height = 5.5, device = cairo_pdf
)
ggsave(file.path(out_dir, "Figure2C_HUA_Prevalence.png"),
  p,
  width = 6, height = 5.5, dpi = 300
)

cat("\nFigure 2C saved.\n")
