# ==============================================================================
# R1.2 / R2.1 canonical enterotype robustness analysis
#
# Reviewer repair, 2026-08-28.
#
# This script deliberately separates:
#   1) assignment confidence conditional on k = 4;
#   2) component-wise recovery under bootstrap perturbation;
#   3) uncertainty in the number of DMM components;
#   4) UA/HUA sensitivity to a k = 3 representation; and
#   5) P/B-ratio sensitivity to zero replacement.
#
# The outcome analyses are NOT used to select k. All seeds and unfavorable
# results are retained. See results/r1_2_r2_1/analysis_protocol.md.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(phyloseq)
  library(microbiome)
  library(DirichletMultinomial)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(clue)
  library(mgcv)
  library(nlme)
})

analysis_seed <- 20260828L
set.seed(analysis_seed)

out_audit <- here("results", "r1_2_r2_1")
out_tables <- here("results", "tables", "supplementary")
out_figures <- here("results", "figures", "panels")
out_cache <- here("results", "cache")
dir.create(out_audit, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(out_cache, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_audit, "execution_log.txt")
sink(log_path, split = TRUE)
on.exit({
  while (sink.number() > 0) sink()
}, add = TRUE)

cat("R1.2 / R2.1 canonical robustness analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Analysis seed:", analysis_seed, "\n\n")

# ---- Figure style -----------------------------------------------------------
ET_COLORS <- c(
  "Type 1" = "#E64B35",
  "Type 2" = "#4DBBD5",
  "Type 3" = "#00A087",
  "Type 4" = "#3C5488"
)

theme_jtm <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, colour = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7.5),
      plot.margin = margin(5, 8, 5, 5),
      strip.background = element_blank(),
      strip.text = element_text(size = 8)
    )
}

save_plot_pair <- function(plot, stem, width, height) {
  ggsave(file.path(out_figures, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out_figures, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 300, device = "png", type = "cairo")
}

# ---- Helpers ----------------------------------------------------------------
compute_centroids <- function(count_mat, assignments, k) {
  rs <- rowSums(count_mat)
  rs[rs == 0] <- 1
  rel <- sweep(count_mat, 1, rs, "/")
  ans <- matrix(0, nrow = k, ncol = ncol(count_mat),
                dimnames = list(seq_len(k), colnames(count_mat)))
  for (j in seq_len(k)) {
    idx <- which(assignments == j)
    if (length(idx) > 0) ans[j, ] <- colMeans(rel[idx, , drop = FALSE])
  }
  ans
}

match_centroids <- function(reference, candidate) {
  k <- nrow(reference)
  cost <- matrix(1, nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      rij <- suppressWarnings(cor(reference[i, ], candidate[j, ],
                                  method = "spearman"))
      if (is.finite(rij)) cost[i, j] <- 1 - rij
    }
  }
  sol <- as.integer(solve_LSAP(cost))
  candidate_to_reference <- integer(k)
  for (i in seq_len(k)) candidate_to_reference[sol[i]] <- i
  candidate_to_reference
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  ta <- rowSums(tab)
  tb <- colSums(tab)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  sum_ij <- sum(choose(tab, 2))
  sum_a <- sum(choose(ta, 2))
  sum_b <- sum(choose(tb, 2))
  expected <- sum_a * sum_b / choose(n, 2)
  maximum <- (sum_a + sum_b) / 2
  if (maximum == expected) return(NA_real_)
  (sum_ij - expected) / (maximum - expected)
}

collapse_duplicate_assignments <- function(ids, assignments) {
  tmp <- split(assignments, ids)
  vapply(tmp, function(z) {
    counts <- table(z)
    as.integer(names(counts)[which.max(counts)])
  }, integer(1))
}

partial_spearman <- function(x, y, covs) {
  d <- data.frame(x = x, y = y, covs, check.names = FALSE)
  d <- d[complete.cases(d), , drop = FALSE]
  rx <- rank(d$x, ties.method = "average")
  ry <- rank(d$y, ties.method = "average")
  rcov <- as.data.frame(lapply(d[, -(1:2), drop = FALSE], function(z) {
    if (is.factor(z) || is.character(z)) z <- as.numeric(factor(z)) - 1
    rank(z, ties.method = "average")
  }))
  fit_x <- lm(rx ~ ., data = rcov)
  fit_y <- lm(ry ~ ., data = rcov)
  r <- cor(residuals(fit_x), residuals(fit_y))
  df <- nrow(d) - 2 - ncol(rcov)
  stat <- r * sqrt(df / (1 - r^2))
  list(rho = r, p = 2 * pt(abs(stat), df = df, lower.tail = FALSE),
       n = nrow(d), df = df)
}

fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "fg")

# ---- Canonical data ---------------------------------------------------------
ps_raw <- readRDS(here("data", "clean", "UA_Analysis_20230811.rds"))
ps_genus <- microbiome::aggregate_taxa(ps_raw, level = "Genus")
otu_mat <- t(as(otu_table(ps_genus), "matrix"))
otu_mat <- matrix(as.integer(round(otu_mat)), nrow = nrow(otu_mat),
                  dimnames = dimnames(otu_mat))

ps_labeled <- readRDS(here("data", "clean",
                           "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta_s4 <- as.data.frame(sample_data(ps_labeled))
# phyloseq::sample_data inherits from data.frame but nlme::gls rejects the
# subclass. Rebuild a plain data.frame without changing values or row names.
meta <- data.frame(lapply(meta_s4, identity), row.names = rownames(meta_s4),
                   check.names = FALSE)
meta <- meta[rownames(otu_mat), , drop = FALSE]
stopifnot(identical(rownames(meta), rownames(otu_mat)))

orig_assign <- as.integer(meta$Dirichlettype_Genus)
names(orig_assign) <- rownames(meta)
meta$Enterotype <- factor(orig_assign, levels = 1:4,
                         labels = paste("Type", 1:4))
meta$Gender <- factor(meta$Gender)
meta$HUA <- ifelse(meta$Hyperuricemia == "Yes", 1, 0)
numeric_vars <- intersect(c("UA", "Age", "BMI", "Cr", "Urea", "GLU",
                            "TG", "HDLC", "LDLC"), names(meta))
for (v in numeric_vars) meta[[v]] <- as.numeric(meta[[v]])

cat("Canonical data:", nrow(otu_mat), "participants x", ncol(otu_mat),
    "genera\n")
print(table(meta$Enterotype))

# =============================================================================
# 1. Posterior assignment confidence, conditional on k=4
# =============================================================================
cat("\n[1] Posterior confidence conditional on k=4\n")
dmm_cache <- readRDS(here("results", "cache", "dmm_model_fits.rds"))
best_fit <- dmm_cache$dmm_fits[[4]]
post <- mixture(best_fit, assign = FALSE)
post_sorted <- t(apply(post, 1, sort, decreasing = TRUE))
posterior <- data.frame(
  sample_id = rownames(post),
  assigned_type = apply(post, 1, which.max),
  top1 = post_sorted[, 1],
  top2 = post_sorted[, 2],
  margin = post_sorted[, 1] - post_sorted[, 2],
  stringsAsFactors = FALSE
)

posterior_summary <- data.frame(
  metric = c("n", "mean_top1", "median_top1", "min_top1",
             "n_top1_below_0.7", "pct_top1_below_0.7",
             "min_top1_top2_margin"),
  value = c(nrow(posterior), mean(posterior$top1), median(posterior$top1),
            min(posterior$top1), sum(posterior$top1 < 0.7),
            100 * mean(posterior$top1 < 0.7), min(posterior$margin))
)
print(posterior_summary)
cat("Least-certain assignment:\n")
print(posterior[which.min(posterior$margin), ])

# Frozen baseline reproduction
stopifnot(abs(mean(posterior$top1) - 0.987740089348957) < 1e-10)

write.csv(posterior, file.path(out_tables, "posterior_assignment_details.csv"),
          row.names = FALSE)

p_post <- ggplot(posterior, aes(x = top1)) +
  geom_histogram(binwidth = 0.025, boundary = 0,
                 fill = "grey70", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0.7, linetype = "dashed",
             colour = "#3C5488", linewidth = 0.55) +
  annotate("text", x = 0.72, y = Inf, vjust = 1.5, hjust = 0,
           label = "0.70", colour = "#3C5488", size = 2.7) +
  annotate("text", x = 0.51, y = Inf, vjust = 3.0, hjust = 0,
           label = sprintf("mean = %.3f\n%d/%d below 0.70",
                           mean(posterior$top1),
                           sum(posterior$top1 < 0.7), nrow(posterior)),
           size = 2.7) +
  coord_cartesian(xlim = c(0.48, 1.005), clip = "off") +
  labs(x = "Maximum posterior assignment probability",
       y = "Number of participants") +
  theme_jtm()
save_plot_pair(p_post, "FigureS_stability_posterior", 3.45, 2.6)

# =============================================================================
# 2. DMM initialization sensitivity and global multi-start minima
# =============================================================================
cat("\n[2] DMM initialization sensitivity\n")
seed_results <- readRDS(here("results", "cache", "seed_search_results.rds"))
seed_results$gap_k4_minus_k3 <- seed_results$lp4 - seed_results$lp3
winner_counts <- as.data.frame(table(seed_results$best_k),
                               stringsAsFactors = FALSE)
names(winner_counts) <- c("k", "n_starts")
winner_counts$k <- as.integer(as.character(winner_counts$k))
winner_counts$pct <- 100 * winner_counts$n_starts / nrow(seed_results)

minima <- data.frame(
  k = 3:5,
  min_laplace = c(min(seed_results$lp3), min(seed_results$lp4),
                  min(seed_results$lp5)),
  seed_at_min = c(seed_results$seed[which.min(seed_results$lp3)],
                  seed_results$seed[which.min(seed_results$lp4)],
                  seed_results$seed[which.min(seed_results$lp5)])
)
minima$delta_from_best <- minima$min_laplace - min(minima$min_laplace)

gap_summary <- data.frame(
  metric = c("n_starts", "mean_k4_minus_k3", "sd_k4_minus_k3",
             "n_k4_better_pairwise", "n_k3_better_pairwise",
             "analysis_seed10_gap"),
  value = c(nrow(seed_results), mean(seed_results$gap_k4_minus_k3),
            sd(seed_results$gap_k4_minus_k3),
            sum(seed_results$gap_k4_minus_k3 < 0),
            sum(seed_results$gap_k4_minus_k3 > 0),
            seed_results$gap_k4_minus_k3[seed_results$seed == 10])
)
print(winner_counts)
print(minima)
print(gap_summary)
stopifnot(identical(winner_counts$n_starts, c(79L, 18L, 3L)))

k_summary <- bind_rows(
  winner_counts %>% transmute(section = "winner_frequency",
                              metric = paste0("k", k), value = n_starts,
                              detail = paste0(fmt(pct, 3), "%")),
  minima %>% transmute(section = "multi_start_minimum",
                       metric = paste0("k", k), value = min_laplace,
                       detail = paste0("seed=", seed_at_min,
                                       "; delta=", fmt(delta_from_best, 5))),
  gap_summary %>% transmute(section = "k4_minus_k3_gap", metric,
                            value, detail = "")
)
write.csv(k_summary, file.path(out_tables, "k_sensitivity_summary.csv"),
          row.names = FALSE)

p_seed_a <- ggplot(winner_counts, aes(x = factor(k), y = n_starts,
                                      fill = factor(k))) +
  geom_col(width = 0.65, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = n_starts), vjust = -0.35, size = 2.8) +
  scale_fill_manual(values = c("3" = "grey65", "4" = "#E64B35",
                               "5" = "grey82"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Laplace-optimal number of components (k)",
       y = "Number of initializations") +
  theme_jtm()

p_seed_b <- ggplot(seed_results, aes(x = gap_k4_minus_k3)) +
  geom_histogram(binwidth = 25, boundary = 0, fill = "grey70",
                 colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#3C5488", linewidth = 0.5) +
  geom_vline(xintercept = mean(seed_results$gap_k4_minus_k3),
             colour = "#E64B35", linewidth = 0.55) +
  annotate("text", x = -165, y = Inf, vjust = 1.4, hjust = 0,
           label = "k=4 favoured\n18/100", size = 2.6,
           colour = "#3C5488") +
  annotate("text", x = 295, y = Inf, vjust = 1.4, hjust = 1,
           label = "k=3 favoured\n82/100", size = 2.6,
           colour = "#3C5488") +
  labs(x = "Laplace(k=4) - Laplace(k=3)\n(negative favours k=4)",
       y = "Number of initializations") +
  theme_jtm()

p_seed <- ggarrange(p_seed_a, p_seed_b, ncol = 2,
                    labels = c("a", "b"), font.label = list(size = 10,
                                                             face = "bold"),
                    widths = c(0.85, 1.2))
save_plot_pair(p_seed, "FigureS_stability_seed", 7.1, 2.75)

# =============================================================================
# 3. Bootstrap component recovery: all 100 seeds, unique in-bag participants
# =============================================================================
cat("\n[3] Bootstrap component recovery conditional on k=4\n")
n_boot <- as.integer(Sys.getenv("R12_N_BOOT", unset = "100"))
force_boot <- identical(tolower(Sys.getenv("R12_FORCE_BOOT", unset = "false")),
                        "true")
boot_cache_path <- file.path(out_cache,
                             paste0("r12_bootstrap_cluster_metrics_", n_boot,
                                    ".rds"))
orig_centroids <- compute_centroids(otu_mat, orig_assign, 4)

if (file.exists(boot_cache_path) && !force_boot) {
  boot_obj <- readRDS(boot_cache_path)
  cat("Loaded bootstrap cache:", boot_cache_path, "\n")
} else {
  iteration_rows <- list()
  cluster_rows <- list()
  failures <- data.frame(iteration = integer(), reason = character())

  for (b in seq_len(n_boot)) {
    set.seed(1000L + b)
    boot_idx <- sample(seq_len(nrow(otu_mat)), replace = TRUE)
    boot_counts <- otu_mat[boot_idx, , drop = FALSE]
    fit_b <- tryCatch(
      dmn(boot_counts, k = 4, seed = b, verbose = FALSE),
      error = function(e) e
    )
    if (inherits(fit_b, "error")) {
      failures <- rbind(failures,
                        data.frame(iteration = b,
                                   reason = conditionMessage(fit_b)))
      cat("bootstrap", b, "failed:", conditionMessage(fit_b), "\n")
      next
    }

    boot_raw <- apply(mixture(fit_b, assign = FALSE), 1, which.max)
    boot_centroids <- compute_centroids(boot_counts, boot_raw, 4)
    relabel_map <- match_centroids(orig_centroids, boot_centroids)
    boot_relabel_rows <- relabel_map[boot_raw]

    ids <- rownames(otu_mat)[boot_idx]
    boot_unique_raw <- collapse_duplicate_assignments(ids, boot_raw)
    unique_ids <- names(boot_unique_raw)
    boot_unique <- relabel_map[boot_unique_raw]
    orig_unique <- orig_assign[unique_ids]

    iteration_rows[[length(iteration_rows) + 1]] <- data.frame(
      iteration = b,
      resample_seed = 1000L + b,
      dmm_seed = b,
      n_unique_inbag = length(unique_ids),
      row_weighted_agreement = mean(boot_relabel_rows == orig_assign[ids]),
      unique_agreement = mean(boot_unique == orig_unique),
      ari = adjusted_rand_index(orig_unique, boot_unique)
    )

    for (cl in 1:4) {
      original_set <- unique_ids[orig_unique == cl]
      recovered_set <- unique_ids[boot_unique == cl]
      intersection_n <- length(intersect(original_set, recovered_set))
      union_n <- length(union(original_set, recovered_set))
      cluster_rows[[length(cluster_rows) + 1]] <- data.frame(
        iteration = b,
        enterotype = paste("Type", cl),
        original_inbag_n = length(original_set),
        recovered_n = length(recovered_set),
        intersection_n = intersection_n,
        jaccard = ifelse(union_n > 0, intersection_n / union_n, NA_real_),
        recall = ifelse(length(original_set) > 0,
                        intersection_n / length(original_set), NA_real_),
        precision = ifelse(length(recovered_set) > 0,
                           intersection_n / length(recovered_set), NA_real_)
      )
    }
    cat(sprintf("bootstrap %03d/%03d: unique agreement %.3f; ARI %.3f\n",
                b, n_boot,
                iteration_rows[[length(iteration_rows)]]$unique_agreement,
                iteration_rows[[length(iteration_rows)]]$ari))
  }

  boot_obj <- list(
    iteration = bind_rows(iteration_rows),
    cluster = bind_rows(cluster_rows),
    failures = failures,
    n_boot_requested = n_boot,
    created = Sys.time()
  )
  saveRDS(boot_obj, boot_cache_path)
}

boot_iterations <- boot_obj$iteration
boot_clusters <- boot_obj$cluster
boot_failures <- boot_obj$failures
stopifnot(nrow(boot_iterations) + nrow(boot_failures) == n_boot)

write.csv(boot_iterations,
          file.path(out_tables, "bootstrap_cluster_stability_iterations.csv"),
          row.names = FALSE)
write.csv(boot_clusters,
          file.path(out_tables, "bootstrap_cluster_stability_by_type.csv"),
          row.names = FALSE)
write.csv(boot_failures,
          file.path(out_audit, "bootstrap_failures.csv"), row.names = FALSE)

summarise_metric <- function(data, metric_name, group = NULL) {
  if (is.null(group)) {
    x <- data[[metric_name]]
    data.frame(group = "overall", metric = metric_name,
               n = sum(is.finite(x)), mean = mean(x, na.rm = TRUE),
               sd = sd(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
               q025 = quantile(x, 0.025, na.rm = TRUE),
               q975 = quantile(x, 0.975, na.rm = TRUE),
               min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
  } else {
    data %>%
      group_by(.data[[group]]) %>%
      summarise(
        metric = metric_name,
        n = sum(is.finite(.data[[metric_name]])),
        mean = mean(.data[[metric_name]], na.rm = TRUE),
        sd = sd(.data[[metric_name]], na.rm = TRUE),
        median = median(.data[[metric_name]], na.rm = TRUE),
        q025 = quantile(.data[[metric_name]], 0.025, na.rm = TRUE),
        q975 = quantile(.data[[metric_name]], 0.975, na.rm = TRUE),
        min = min(.data[[metric_name]], na.rm = TRUE),
        max = max(.data[[metric_name]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(group = all_of(group))
  }
}

boot_summary <- bind_rows(
  summarise_metric(boot_iterations, "row_weighted_agreement"),
  summarise_metric(boot_iterations, "unique_agreement"),
  summarise_metric(boot_iterations, "ari"),
  summarise_metric(boot_clusters, "jaccard", "enterotype"),
  summarise_metric(boot_clusters, "recall", "enterotype"),
  summarise_metric(boot_clusters, "precision", "enterotype")
)
write.csv(boot_summary,
          file.path(out_tables, "bootstrap_cluster_stability_summary.csv"),
          row.names = FALSE)
print(boot_summary)

if (n_boot == 100) {
  stopifnot(nrow(boot_failures) == 0)
  # Reproduce the legacy row-weighted baseline before replacing it. Four of
  # 100 local-optimum fits differ from the 2026-08-23 cache in the current R
  # package build; the mean differs by 0.39% relative, within the protocol's
  # <=1% baseline-reproduction tolerance. Record every differing iteration.
  legacy_rates_path <- file.path(out_cache,
                                 "bootstrap_stability_100_rates.rds")
  if (file.exists(legacy_rates_path)) {
    legacy_rates <- readRDS(legacy_rates_path)
    baseline_compare <- data.frame(
      iteration = seq_len(100),
      legacy_row_agreement = legacy_rates,
      current_row_agreement = boot_iterations$row_weighted_agreement
    )
    baseline_compare$difference <- baseline_compare$current_row_agreement -
      baseline_compare$legacy_row_agreement
    write.csv(baseline_compare,
              file.path(out_audit,
                        "legacy_bootstrap_baseline_comparison.csv"),
              row.names = FALSE)
    legacy_mean <- mean(legacy_rates)
    current_mean <- mean(boot_iterations$row_weighted_agreement)
    relative_difference <- abs(current_mean - legacy_mean) / legacy_mean
    cat(sprintf("Legacy bootstrap mean %.6f; current %.6f; relative difference %.3f%%; differing iterations %d/100\n",
                legacy_mean, current_mean, 100 * relative_difference,
                sum(abs(baseline_compare$difference) > 1e-12)))
    stopifnot(relative_difference <= 0.01)
  }
}

# R1.2 bootstrap recovery: give each cited metric a matching visual. Panel (a)
# shows unique-participant agreement across the 100 refits, panel (b) the ARI,
# and panel (c) the component-specific Jaccard similarity.
boot_agree_lab <- sprintf("mean = %.1f%%\nSD = %.1f%%",
                          100 * mean(boot_iterations$unique_agreement),
                          100 * sd(boot_iterations$unique_agreement))
p_boot_agree <- ggplot(boot_iterations, aes(x = unique_agreement)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "grey70",
                 colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = mean(boot_iterations$unique_agreement),
             colour = "#3C5488", linewidth = 0.55) +
  annotate("text", x = min(boot_iterations$unique_agreement), y = Inf,
           hjust = 0, vjust = 1.6, label = boot_agree_lab,
           colour = "#3C5488", size = 2.4) +
  coord_cartesian(clip = "off") +
  labs(x = "Unique-participant agreement", y = "Number of refits") +
  theme_jtm()

boot_ari_lab <- sprintf("median = %.3f",
                        median(boot_iterations$ari))
p_boot_ari <- ggplot(boot_iterations, aes(x = ari)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "grey70",
                 colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = median(boot_iterations$ari),
             colour = "#3C5488", linewidth = 0.55) +
  annotate("text", x = min(boot_iterations$ari), y = Inf,
           hjust = 0, vjust = 1.6, label = boot_ari_lab,
           colour = "#3C5488", size = 2.4) +
  coord_cartesian(clip = "off") +
  labs(x = "Adjusted Rand index (ARI)", y = "Number of refits") +
  theme_jtm()

jaccard_labels <- boot_summary %>%
  filter(metric == "jaccard") %>%
  transmute(enterotype = factor(group, levels = paste("Type", 1:4)),
            label = sprintf("median %.2f", median), y = 1.015)

p_boot_jaccard <- ggplot(boot_clusters,
                 aes(x = factor(enterotype, levels = paste("Type", 1:4)),
                     y = jaccard, fill = enterotype)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", colour = "grey45",
             linewidth = 0.45) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = "grey25",
               linewidth = 0.35) +
  geom_jitter(width = 0.13, height = 0, size = 0.65, alpha = 0.32,
              colour = "grey15") +
  geom_text(data = jaccard_labels,
            aes(x = enterotype, y = y, label = label), inherit.aes = FALSE,
            size = 2.45, vjust = 0) +
  scale_fill_manual(values = ET_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Bootstrap Jaccard similarity") +
  theme_jtm()

p_boot <- ggarrange(p_boot_agree, p_boot_ari, p_boot_jaccard, ncol = 3,
                    widths = c(1, 1, 1.3))
save_plot_pair(p_boot, "FigureS_stability_bootstrap", 6.8, 2.6)

# =============================================================================
# 4. k=3 sensitivity of UA/HUA outcomes (not model selection)
# =============================================================================
cat("\n[4] UA/HUA sensitivity under the best-attained k=3 fit\n")
fit_k3 <- dmn(otu_mat, k = 3, seed = 27, verbose = FALSE)
laplace_k3 <- as.numeric(fit_k3@goodnessOfFit["Laplace"])
stopifnot(abs(laplace_k3 - minima$min_laplace[minima$k == 3]) < 1e-6)
assign_k3 <- apply(mixture(fit_k3, assign = FALSE), 1, which.max)
names(assign_k3) <- rownames(otu_mat)
meta$k3 <- factor(assign_k3[rownames(meta)], levels = 1:3,
                  labels = c("C1", "C2", "C3"))
meta$k3 <- relevel(meta$k3, ref = "C2")

k4_k3_cross <- as.data.frame.matrix(table(meta$Enterotype, meta$k3))
k4_k3_cross$k4_enterotype <- rownames(k4_k3_cross)
k4_k3_cross <- k4_k3_cross %>%
  select(k4_enterotype, everything())
write.csv(k4_k3_cross,
          file.path(out_tables, "k4_k3_membership_crosswalk.csv"),
          row.names = FALSE)
cat("k=4 by k=3 crosswalk:\n")
print(table(meta$Enterotype, meta$k3))

run_k3_models <- function(label, covariates) {
  needed <- c("UA", "HUA", "k3", covariates)
  d <- meta[complete.cases(meta[, needed]), , drop = FALSE]
  rhs <- paste(c("k3", covariates), collapse = " + ")

  gls_fit <- gls(as.formula(paste("UA ~", rhs)), data = d,
                 weights = varIdent(form = ~ 1 | k3), method = "REML")
  gls_tab <- summary(gls_fit)$tTable
  gls_ci <- intervals(gls_fit, which = "coef")$coef
  gls_global <- anova(gls_fit, type = "marginal")["k3", "p-value"]

  logit_fit <- glm(as.formula(paste("HUA ~", rhs)), data = d,
                   family = binomial)
  logit_reduced <- update(logit_fit, . ~ . - k3)
  logit_global <- anova(logit_reduced, logit_fit, test = "LRT")$`Pr(>Chi)`[2]
  logit_tab <- summary(logit_fit)$coefficients
  logit_ci <- confint.default(logit_fit)

  rows <- list()
  for (term in c("k3C1", "k3C3")) {
    contrast <- ifelse(term == "k3C1", "C1_vs_C2", "C3_vs_C2")
    rows[[length(rows) + 1]] <- data.frame(
      adjustment = label, outcome = "UA", contrast = contrast,
      n = nrow(d), estimate = gls_tab[term, "Value"],
      ci_lower = gls_ci[term, "lower"], ci_upper = gls_ci[term, "upper"],
      p_value = gls_tab[term, "p-value"], global_p = gls_global,
      scale = "beta_umol_L"
    )
    rows[[length(rows) + 1]] <- data.frame(
      adjustment = label, outcome = "HUA", contrast = contrast,
      n = nrow(d), estimate = exp(logit_tab[term, "Estimate"]),
      ci_lower = exp(logit_ci[term, "2.5 %"]),
      ci_upper = exp(logit_ci[term, "97.5 %"]),
      p_value = logit_tab[term, "Pr(>|z|)"], global_p = logit_global,
      scale = "odds_ratio"
    )
  }
  bind_rows(rows)
}

k3_outcomes <- bind_rows(
  run_k3_models("sex_age", c("Gender", "Age")),
  run_k3_models("renal_cardiometabolic",
                c("Gender", "Age", "Cr", "Urea", "GLU", "TG", "HDLC",
                  "LDLC"))
)
k3_desc <- meta %>%
  group_by(k3) %>%
  summarise(n = n(), mean_UA = mean(UA), sd_UA = sd(UA),
            HUA_n = sum(HUA), HUA_pct = 100 * mean(HUA), .groups = "drop")
write.csv(k3_outcomes,
          file.path(out_tables, "k3_outcome_sensitivity.csv"), row.names = FALSE)
write.csv(k3_desc,
          file.path(out_tables, "k3_outcome_descriptives.csv"), row.names = FALSE)
print(k3_desc)
print(k3_outcomes)

# =============================================================================
# 5. Continuous P/B ratio, zero handling, and model comparison
# =============================================================================
cat("\n[5] Continuous P/B ratio sensitivity\n")
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
rel_mat <- t(as(otu_table(ps_rel), "matrix"))
if (!identical(rownames(rel_mat), rownames(meta))) {
  rel_mat <- rel_mat[rownames(meta), , drop = FALSE]
}
tax <- as.data.frame(tax_table(ps_rel))
genus_names <- gsub("^g__", "", as.character(tax$Genus))
prev <- rowSums(rel_mat[, genus_names == "Prevotella", drop = FALSE])
bact <- rowSums(rel_mat[, genus_names == "Bacteroides", drop = FALSE])
stopifnot(all(bact > 0))

half_min_prev <- min(prev[prev > 0]) / 2
ratio_specs <- list(
  pc_1e_6 = list(pc = 1e-6, detected_only = FALSE),
  pc_1e_5 = list(pc = 1e-5, detected_only = FALSE),
  pc_1e_4 = list(pc = 1e-4, detected_only = FALSE),
  pc_1e_3 = list(pc = 1e-3, detected_only = FALSE),
  half_min_positive = list(pc = half_min_prev, detected_only = FALSE),
  prevotella_detected = list(pc = 1e-6, detected_only = TRUE)
)

run_pb_spec <- function(name, spec) {
  ratio <- log2(prev + spec$pc) - log2(bact + spec$pc)
  keep <- if (spec$detected_only) prev > 0 else rep(TRUE, length(prev))
  d0 <- meta[keep, , drop = FALSE]
  d0$pb <- ratio[keep]
  d_ua <- d0[complete.cases(d0[, c("UA", "pb", "Age", "Gender",
                                    "Enterotype")]), , drop = FALSE]
  d_hua <- d0[complete.cases(d0[, c("HUA", "pb", "Age", "Gender",
                                     "BMI", "Cr")]), , drop = FALSE]

  ps <- partial_spearman(d_ua$pb, d_ua$UA,
                         d_ua[, c("Age", "Gender")])

  # Primary nonlinear check: fixed-df natural splines nested in the same
  # Gaussian linear model. This directly tests added curvature degrees of
  # freedom and avoids interpreting a penalized smooth's overall test as a
  # test of nonlinearity.
  linear <- lm(UA ~ pb + Age + Gender, data = d_ua)
  spline3 <- lm(UA ~ splines::ns(pb, df = 3) + Age + Gender, data = d_ua)
  spline4 <- lm(UA ~ splines::ns(pb, df = 4) + Age + Gender, data = d_ua)
  nl3 <- anova(linear, spline3)
  nl4 <- anova(linear, spline4)

  # Retired diagnostic retained for auditability. The previous workflow used
  # anova.gam(linear, penalized smooth) as a nonlinearity test. The fitted
  # smooth has edf approximately 1, so that P-value does not support curvature.
  legacy_linear <- gam(UA ~ pb + Age + Gender, data = d_ua, method = "ML")
  legacy_smooth <- gam(UA ~ s(pb, bs = "cr", k = 10) + Age + Gender,
                       data = d_ua, method = "ML")
  legacy_cmp <- anova(legacy_linear, legacy_smooth, test = "LRT")

  logit <- glm(HUA ~ pb + Age + Gender + BMI + Cr,
               data = d_hua, family = binomial)
  et_model <- lm(UA ~ Enterotype + Age + Gender, data = d_ua)
  legacy_st <- summary(legacy_smooth)$s.table

  data.frame(
    specification = name,
    pseudocount = spec$pc,
    detected_only = spec$detected_only,
    n_ua = nrow(d_ua),
    n_hua = nrow(d_hua),
    prevotella_zero_n = sum(prev[keep] == 0),
    partial_rho = ps$rho,
    partial_p = ps$p,
    nonlinear_ns3_p = nl3$`Pr(>F)`[2],
    nonlinear_ns4_p = nl4$`Pr(>F)`[2],
    legacy_gam_edf = legacy_st[1, "edf"],
    legacy_gam_comparison_p = legacy_cmp$`Pr(>Chi)`[2],
    hua_or = exp(coef(logit)["pb"]),
    hua_ci_lower = exp(confint.default(logit)["pb", 1]),
    hua_ci_upper = exp(confint.default(logit)["pb", 2]),
    hua_p = summary(logit)$coefficients["pb", "Pr(>|z|)"],
    aic_enterotype = AIC(et_model),
    aic_pb_linear = AIC(linear),
    aic_pb_spline3 = AIC(spline3),
    aic_pb_spline4 = AIC(spline4),
    stringsAsFactors = FALSE
  )
}

pb_sensitivity <- bind_rows(lapply(names(ratio_specs), function(nm) {
  run_pb_spec(nm, ratio_specs[[nm]])
}))
write.csv(pb_sensitivity,
          file.path(out_tables, "pb_ratio_sensitivity.csv"), row.names = FALSE)
print(pb_sensitivity)

main_pb <- pb_sensitivity[pb_sensitivity$specification == "pc_1e_6", ]
stopifnot(abs(main_pb$partial_rho - (-0.151825421762325)) < 1e-10)
stopifnot(abs(main_pb$partial_p - 0.0081149509560736) < 1e-10)
stopifnot(abs(main_pb$aic_enterotype - 3515.87234428002) < 1e-8)
stopifnot(abs(main_pb$aic_pb_linear - 3516.12582305862) < 1e-8)

# Adjusted linear and fixed-df spline fits for the primary specification; raw
# points are shown only to display observed data density. Predictions hold
# covariates constant.
main_ratio <- log2(prev + 1e-6) - log2(bact + 1e-6)
plot_df <- meta
plot_df$pb <- main_ratio
linear_plot <- lm(UA ~ pb + Age + Gender, data = plot_df)
spline_plot <- lm(UA ~ splines::ns(pb, df = 3) + Age + Gender,
                  data = plot_df)
modal_gender <- names(which.max(table(plot_df$Gender)))
grid <- data.frame(
  pb = seq(min(plot_df$pb), max(plot_df$pb), length.out = 300),
  Age = median(plot_df$Age, na.rm = TRUE),
  Gender = factor(modal_gender, levels = levels(plot_df$Gender))
)
pred_linear <- predict(linear_plot, newdata = grid, se.fit = TRUE)
pred_spline <- predict(spline_plot, newdata = grid, se.fit = TRUE)
grid$linear <- as.numeric(pred_linear$fit)
grid$linear_lower <- grid$linear - 1.96 * as.numeric(pred_linear$se.fit)
grid$linear_upper <- grid$linear + 1.96 * as.numeric(pred_linear$se.fit)
grid$spline <- as.numeric(pred_spline$fit)

p_pb <- ggplot(plot_df, aes(x = pb, y = UA)) +
  geom_point(aes(colour = Enterotype), alpha = 0.38, size = 0.85) +
  geom_ribbon(data = grid,
              aes(x = pb, ymin = linear_lower, ymax = linear_upper),
              inherit.aes = FALSE, fill = "#3C5488", alpha = 0.14) +
  geom_line(data = grid, aes(x = pb, y = linear), inherit.aes = FALSE,
            colour = "#3C5488", linewidth = 0.65) +
  geom_line(data = grid, aes(x = pb, y = spline), inherit.aes = FALSE,
            colour = "#E64B35", linewidth = 0.55, linetype = "dashed") +
  annotate("text", x = min(plot_df$pb), y = max(plot_df$UA),
           hjust = 0, vjust = 1,
           label = sprintf("partial Spearman rho = %.3f, P = %.3f\nnonlinearity (natural spline, df=3) P = %.3f; Prevotella zeros = %d/%d",
                           main_pb$partial_rho, main_pb$partial_p,
                           main_pb$nonlinear_ns3_p, sum(prev == 0), length(prev)),
           size = 2.45) +
  scale_colour_manual(values = ET_COLORS) +
  labs(x = expression(log[2] * "(" * italic(Prevotella) * "/" *
                        italic(Bacteroides) * ")"),
       y = expression("Serum UA (" * mu * "mol/L)"), colour = "Enterotype") +
  theme_jtm() +
  theme(legend.position = "right")
save_plot_pair(p_pb, "FigureS_PB_ratio_UA", 5.2, 3.15)

# Complete Supplementary Figure S4 composite. Individual panel files are kept
# for production editing, while this file prevents A-C integration omissions.
# Layout: row 1 = A (posterior) + C (k-sensitivity, two unlabelled sub-panels),
# row 2 = B (bootstrap, three unlabelled sub-panels); each row holds three plots.
p_seed_clean <- ggarrange(p_seed_a, p_seed_b, ncol = 2,
                          widths = c(0.85, 1.2))
p_s4_top <- ggarrange(p_post, p_seed_clean, ncol = 2,
                      labels = c("A", "C"),
                      font.label = list(size = 10, face = "bold"),
                      widths = c(0.75, 1.25))
p_s4 <- ggarrange(p_s4_top, p_boot, ncol = 1,
                  labels = c("", "B"),
                  font.label = list(size = 10, face = "bold"),
                  heights = c(1.0, 1.0))
save_plot_pair(p_s4, "FigureS1_enterotype_robustness", 9.0, 6.0)

# =============================================================================
# 6. Consolidated reviewer-facing sensitivity table and environment
# =============================================================================
cat("\n[6] Consolidating evidence\n")
reviewer_table <- bind_rows(
  posterior_summary %>%
    transmute(domain = "posterior_k4", metric, group = "overall", value,
              secondary = NA_real_, note = "conditional_on_k4"),
  boot_summary %>%
    transmute(domain = "bootstrap_k4", metric, group,
              value = median, secondary = mean,
              note = paste0("median; mean=", fmt(mean), "; sd=", fmt(sd))),
  k_summary %>%
    transmute(domain = "k_sensitivity", metric,
              group = section, value, secondary = NA_real_, note = detail),
  k3_outcomes %>%
    transmute(domain = "k3_outcome", metric = paste(outcome, contrast, sep = ":"),
              group = adjustment, value = estimate, secondary = p_value,
              note = paste0("95% CI [", fmt(ci_lower), ", ", fmt(ci_upper),
                            "]; global P=", fmt(global_p), "; ", scale)),
  pb_sensitivity %>%
    transmute(domain = "pb_sensitivity", metric = specification,
              group = ifelse(detected_only, "detected_only", "all_samples"),
              value = partial_rho, secondary = nonlinear_ns3_p,
              note = paste0("n_UA=", n_ua, "; n_HUA=", n_hua,
                            "; partial P=", fmt(partial_p),
                            "; HUA OR=", fmt(hua_or), ", P=", fmt(hua_p),
                            "; AIC ET/linear/ns3=", fmt(aic_enterotype), "/",
                            fmt(aic_pb_linear), "/", fmt(aic_pb_spline3)))
)
write.csv(reviewer_table,
          file.path(out_tables, "enterotype_robustness_sensitivity.csv"),
          row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_audit, "session_info.txt"))
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Canonical outputs written to:\n  ", out_tables, "\n  ", out_figures,
    "\n  ", out_audit, "\n", sep = "")

sink()
