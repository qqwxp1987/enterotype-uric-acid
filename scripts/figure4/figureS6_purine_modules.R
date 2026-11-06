# Supplementary Figure S6: purine module scores (EC-union share) by enterotype
# Input : results/tables/section_3_4/purine_modules/02_module_scores_by_sample.csv
# Output: supplementary/figures/FigureS6_module_scores.{pdf,png}
# Colors: global enterotype palette (scripts/utils/palettes.R), consistent with all main figures
suppressPackageStartupMessages({ library(ggplot2); library(tidyr); library(dplyr); library(here) })
set.seed(42)

ET_COLORS <- c("Type 1" = "#E64B35", "Type 2" = "#4DBBD5",
               "Type 3" = "#00A087", "Type 4" = "#3C5488")

df <- read.csv(here("results", "tables", "section_3_4", "purine_modules", "02_module_scores_by_sample.csv"))
df$Enterotype <- factor(df$Enterotype, levels = paste("Type", 1:4))

long <- df %>%
  select(Enterotype, E_de_novo_share, E_degradation_share, E_salvage_share) %>%
  pivot_longer(-Enterotype, names_to = "module", values_to = "share") %>%
  mutate(module = factor(module,
    levels = c("E_de_novo_share", "E_degradation_share", "E_salvage_share"),
    labels = c("De novo biosynthesis", "Degradation", "Salvage")))

p <- ggplot(long, aes(Enterotype, share, colour = Enterotype)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = NA, linewidth = 0.6) +
  geom_jitter(width = 0.16, size = 0.9, alpha = 0.5) +
  scale_colour_manual(values = ET_COLORS, guide = "none") +
  facet_wrap(~ module, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Module score (% of total EC abundance)") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13))

ggsave(here("supplementary", "figures", "FigureS6_module_scores.pdf"), p, width = 9.5, height = 3.6)
ggsave(here("supplementary", "figures", "FigureS6_module_scores.png"), p, width = 9.5, height = 3.6, dpi = 600)
cat("saved FigureS6_module_scores.pdf / .png\n")
