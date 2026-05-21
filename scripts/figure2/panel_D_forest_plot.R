# ==============================================================================
# Figure 2D: Dual Forest Plot
# Top: WLS regression β coefficients + 95% CI (Serum UA)
# Bottom: Logistic regression OR + 95% CI (Hyperuricemia)
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(nlme)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")

meta$Enterotype <- factor(
  meta$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)
meta$Enterotype <- relevel(meta$Enterotype, ref = "Type 1")
meta$Gender <- factor(meta$Gender, levels = c("Female", "Male"))
meta$Age    <- as.numeric(meta$Age)
meta$UA     <- as.numeric(meta$UA)
meta$Hyperuricemia <- factor(meta$Hyperuricemia, levels = c("No", "Yes"))

# --- 2. WLS Model ---
wls <- gls(UA ~ Enterotype + Gender + Age, data = meta,
           weights = varIdent(form = ~ 1 | Enterotype),
           method = "REML", na.action = na.omit)

wls_sum <- summary(wls)$tTable
wls_anova <- anova(wls, type = "marginal")

cat("=== WLS Coefficients ===\n")
print(round(wls_sum, 4))
cat("\n=== WLS Marginal ANOVA ===\n")
print(wls_anova)

# Extract WLS β and 95% CI
wls_terms <- rownames(wls_sum)
wls_df <- data.frame(
  term = wls_terms,
  beta = wls_sum[, "Value"],
  se = wls_sum[, "Std.Error"],
  p = wls_sum[, "p-value"],
  stringsAsFactors = FALSE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    ci_low = beta - 1.96 * se,
    ci_high = beta + 1.96 * se,
    model = "Serum UA (GLS \u03b2)"
  )

# Clean term labels
wls_df$label <- case_when(
  wls_df$term == "EnterotypeType 2" ~ "Type 2 vs Type 1",
  wls_df$term == "EnterotypeType 3" ~ "Type 3 vs Type 1",
  wls_df$term == "EnterotypeType 4" ~ "Type 4 vs Type 1",
  wls_df$term == "GenderMale" ~ "Sex (Male vs Female)",
  wls_df$term == "Age" ~ "Age (per year)",
  TRUE ~ wls_df$term
)

# Add overall Enterotype p-value as annotation
wls_enterotype_p <- wls_anova["Enterotype", "p-value"]

# --- 3. Logistic Model ---
logit <- glm(Hyperuricemia ~ Enterotype + Gender + Age, data = meta,
             family = binomial)

logit_sum <- summary(logit)$coefficients
logit_or <- exp(coef(logit))
logit_ci <- exp(confint.default(logit))

cat("\n=== Logistic Regression ===\n")
cat("OR and 95% CI:\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

# Likelihood ratio test
logit_null <- glm(Hyperuricemia ~ Gender + Age, data = meta, family = binomial)
lr_test <- anova(logit_null, logit, test = "Chisq")
logit_enterotype_p <- lr_test$`Pr(>Chi)`[2]
cat(sprintf("\nLR test for Enterotype: P = %.4e\n", logit_enterotype_p))

# Build logistic data frame
logit_terms <- rownames(logit_sum)[-1]  # exclude intercept
logit_df <- data.frame(
  term = logit_terms,
  beta = logit_or[logit_terms],  # OR
  ci_low = logit_ci[logit_terms, 1],
  ci_high = logit_ci[logit_terms, 2],
  p = logit_sum[logit_terms, "Pr(>|z|)"],
  stringsAsFactors = FALSE
) %>%
  mutate(
    se = NA,
    model = "HUA (Logistic OR)"
  )

logit_df$label <- case_when(
  logit_df$term == "EnterotypeType 2" ~ "Type 2 vs Type 1",
  logit_df$term == "EnterotypeType 3" ~ "Type 3 vs Type 1",
  logit_df$term == "EnterotypeType 4" ~ "Type 4 vs Type 1",
  logit_df$term == "GenderMale" ~ "Sex (Male vs Female)",
  logit_df$term == "Age" ~ "Age (per year)",
  TRUE ~ logit_df$term
)

# --- 4. Combined Data ---
# For WLS: reference line at 0, values are β
# For Logistic: reference line at 1, values are OR
# We'll use separate panels (facets) with free x-scales

# Format p-value labels
fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
  ifelse(p < 0.01, sprintf("%.3f", p),
  sprintf("%.2f", p)))
}

wls_plot_df <- wls_df %>%
  select(label, beta, ci_low, ci_high, p, model) %>%
  rename(estimate = beta) %>%
  mutate(
    p_label = fmt_p(p),
    sig = ifelse(p < 0.05, "sig", "ns"),
    ref_line = 0
  )

logit_plot_df <- logit_df %>%
  select(label, beta, ci_low, ci_high, p, model) %>%
  rename(estimate = beta) %>%
  mutate(
    p_label = fmt_p(p),
    sig = ifelse(p < 0.05, "sig", "ns"),
    ref_line = 1
  )

# Reverse order for plotting (top-to-bottom reading)
var_order <- rev(c("Type 2 vs Type 1", "Type 3 vs Type 1", "Type 4 vs Type 1",
                    "Sex (Male vs Female)", "Age (per year)"))

wls_plot_df$label <- factor(wls_plot_df$label, levels = var_order)
logit_plot_df$label <- factor(logit_plot_df$label, levels = var_order)

# --- 5. Plot WLS Panel ---
p_wls <- ggplot(wls_plot_df, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = sig),
                 height = 0.2, linewidth = 0.7) +
  geom_point(aes(color = sig, shape = sig), size = 3) +
  geom_text(aes(label = p_label, x = ci_high),
            hjust = -0.2, size = 3, color = "grey30") +
  scale_color_manual(values = c("sig" = "#B2182B", "ns" = "grey50"), guide = "none") +
  scale_shape_manual(values = c("sig" = 16, "ns" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    x = expression(beta ~ "(95% CI)"),
    y = NULL,
    subtitle = bquote("Serum UA (GLS regression)" ~~
                        "Enterotype" ~ italic(P) == .(sprintf("%.3f", wls_enterotype_p)))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 9),
    plot.subtitle = element_text(size = 10, face = "bold", color = "grey20"),
    plot.margin = margin(5, 15, 5, 5)
  )

# --- 6. Plot Logistic Panel ---
p_logit <- ggplot(logit_plot_df, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = sig),
                 height = 0.2, linewidth = 0.7) +
  geom_point(aes(color = sig, shape = sig), size = 3) +
  geom_text(aes(label = p_label, x = ci_high),
            hjust = -0.2, size = 3, color = "grey30") +
  scale_color_manual(values = c("sig" = "#B2182B", "ns" = "grey50"), guide = "none") +
  scale_shape_manual(values = c("sig" = 16, "ns" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = NULL,
    subtitle = bquote("Hyperuricemia (Logistic regression)" ~~
                        "Enterotype LR" ~ italic(P) == .(sprintf("%.3f", logit_enterotype_p)))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 9),
    plot.subtitle = element_text(size = 10, face = "bold", color = "grey20"),
    plot.margin = margin(5, 15, 5, 5)
  )

# --- 7. Combine ---
library(patchwork)

p_combined <- p_wls / p_logit +
  plot_layout(heights = c(1, 1))

# --- 8. Save ---
out_dir <- here("results", "figures", "panels")

ggsave(file.path(out_dir, "Figure2D_Forest_Plot.pdf"),
       p_combined, width = 6.5, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure2D_Forest_Plot.png"),
       p_combined, width = 6.5, height = 6, dpi = 300)

cat("\nFigure 2D saved.\n")

# --- 9. Save regression stats ---
table_dir <- here("results", "tables", "section_3_2")

wls_out <- wls_df %>% select(label, beta, se, ci_low, ci_high, p)
write.csv(wls_out, file.path(table_dir, "wls_regression_coefficients.csv"),
          row.names = FALSE)

logit_out <- logit_df %>% select(label, beta, ci_low, ci_high, p) %>%
  rename(OR = beta)
write.csv(logit_out, file.path(table_dir, "logistic_regression_OR.csv"),
          row.names = FALSE)

cat("Regression tables saved.\n")

# --- 10. All pairwise contrasts from WLS via emmeans (adjusted for Sex + Age) ---
library(emmeans)

emm <- emmeans(wls, ~ Enterotype)
pw_contrasts <- pairs(emm, adjust = "holm")
pw_confint <- confint(pw_contrasts)
pw_summary <- summary(pw_contrasts)

contrasts_df <- as.data.frame(pw_summary) %>%
  left_join(
    as.data.frame(pw_confint) %>% select(contrast, lower.CL, upper.CL),
    by = "contrast"
  ) %>%
  mutate(
    significant = ifelse(p.value < 0.05, "Yes", "No"),
    signif_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  arrange(p.value)

cat("\n=== WLS emmeans: All pairwise contrasts (Holm-corrected) ===\n")
print(as.data.frame(contrasts_df))

write.csv(contrasts_df,
          file.path(table_dir, "wls_emmeans_pairwise_UA.csv"),
          row.names = FALSE)
cat("WLS emmeans pairwise contrasts saved.\n")
