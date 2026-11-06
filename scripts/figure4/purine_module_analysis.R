# Purine-module analysis
#
# Uses the frozen pathway register and exact pathway composition to calculate
# EC-union module scores, enterotype contrasts and module-UA partial correlations.

suppressPackageStartupMessages({
  library(phyloseq); library(nlme); library(emmeans); library(ppcor); library(here)
})

RES <- here("results", "tables", "section_3_4", "purine_modules")
VER <- file.path(RES, "verification")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(VER, showWarnings = FALSE, recursive = TRUE)

runlog <- file.path(VER, "run_log.txt")
con <- file(runlog, open = "wt"); sink(con, split = TRUE)
cat("Purine module analysis run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# ---- QC assertion harness -------------------------------------------------
qc <- data.frame(id = character(), pass = logical(), detail = character(),
                 stringsAsFactors = FALSE)
check <- function(id, pass, detail = "") {
  qc <<- rbind(qc, data.frame(id = id, pass = isTRUE(pass), detail = detail))
  cat(sprintf("[QC] %-28s %s  %s\n", id, ifelse(isTRUE(pass), "PASS", "FAIL"), detail))
  if (!isTRUE(pass)) {
    sink(); close(con)
    stop("QC assertion failed: ", id, " — ", detail, call. = FALSE)
  }
}

# ---- inputs ----------------------------------------------------------------
pw <- readRDS(here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"))
ec <- readRDS(here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"))
reg <- read.csv(here("data", "curated", "purine", "frozen_inclusion_register.csv"),
                fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
comp <- read.csv(here("data", "curated", "purine", "pathway_composition.csv"),
                 fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
linda <- read.csv(here("results", "tables", "section_3_4", "linda_pathway_all_comparisons.csv"),
                  stringsAsFactors = FALSE)

# ---- P1: register <-> phyloseq, ID-exact matching --------------------------
pwy_id <- sub(":.*$", "", taxa_names(pw))
check("P1a_474_id_match", identical(sort(pwy_id), sort(reg$pathway_id)),
      "474 phyloseq pathway IDs == 474 register pathway_ids")
inc <- reg[reg$final_action != "exclude", ]
check("P1b_purine_set_34", nrow(inc) == 34, "frozen purine set = 34")
cores <- reg[reg$final_sum_role == "core", ]
check("P1c_cores_15", nrow(cores) == 15 &&
        identical(as.numeric(table(cores$final_module)[c("de_novo", "degradation", "salvage")]),
                  c(3, 11, 1)), "cores = 15 (de_novo 3, degradation 11, salvage 1)")

# ---- abundance matrices ----------------------------------------------------
P <- as(otu_table(pw), "matrix"); if (!taxa_are_rows(pw)) P <- t(P)
E <- as(otu_table(ec), "matrix"); if (!taxa_are_rows(ec)) E <- t(E)
rownames(P) <- pwy_id
check("P2b_structure", nrow(P) == 474 && ncol(P) == 305, "pathway matrix 474 x 305")

sdat <- as.data.frame(sample_data(pw))
sdat$Enterotype <- factor(paste0("Type ", sdat$Dirichlettype_Genus),
                          levels = paste("Type", 1:4))
check("P2b_enterotype_n",
      identical(as.integer(table(sdat$Enterotype)), c(120L, 98L, 48L, 39L)),
      "enterotype counts 120/98/48/39")
check("P2b_UA_complete", !any(is.na(sdat$UA)), "UA available for all 305 samples")

# ---- denominator provenance (分母来源检验) ----------------------------------
# phyloseq per-sample totals must equal totals recomputed from the HUMAnN tsv
tsv <- read.delim(here("data", "clean", "pathabundance_stratified.tsv"),
                  check.names = FALSE, comment.char = "")
names(tsv)[1] <- "feature"
tsv <- tsv[!grepl("\\|", tsv$feature), ]                    # unstratified rows only
tsv <- tsv[!tsv$feature %in% c("UNMAPPED", "UNINTEGRATED"), ]
tsv_id <- sub(":.*$", "", tsv$feature)
tsv_tot <- colSums(tsv[, -1, drop = FALSE])
phy_tot <- colSums(P)
common <- intersect(names(tsv_tot), colnames(P))
check("P2e_denominator", length(common) == 305 &&
        max(abs(tsv_tot[common] - phy_tot[common]) / phy_tot[common]) < 1e-6,
      "phyloseq totals == tsv pathway totals (rel diff < 1e-6)")
den_df <- data.frame(sample = names(phy_tot), total_pathway = as.numeric(phy_tot),
                     Enterotype = sdat$Enterotype)
den_by_et <- aggregate(total_pathway ~ Enterotype, den_df,
                       function(x) c(n = length(x), median = median(x), IQR = IQR(x)))
write.csv(data.frame(Enterotype = den_by_et$Enterotype,
                     n = den_by_et$total_pathway[, "n"],
                     median_total = den_by_et$total_pathway[, "median"],
                     IQR_total = den_by_et$total_pathway[, "IQR"]),
          file.path(VER, "v06_denominator_by_enterotype.csv"), row.names = FALSE)

# ---- module membership from the frozen register ----------------------------
sco <- reg[reg$in_pathway_module_score %in% c("primary", "yes", "sensitivity"), ]
mod_pwy <- list(
  de_novo_primary     = sco$pathway_id[sco$final_module == "de_novo" & sco$in_pathway_module_score == "primary"],
  de_novo_sensitivity = sco$pathway_id[sco$final_module == "de_novo" & sco$in_pathway_module_score == "sensitivity"],
  degradation         = sco$pathway_id[sco$final_module == "degradation" & sco$in_pathway_module_score == "yes"],
  salvage             = sco$pathway_id[sco$final_module == "salvage" & sco$in_pathway_module_score == "yes"]
)
check("P1d_score_membership",
      identical(lengths(mod_pwy), c(de_novo_primary = 1L, de_novo_sensitivity = 1L,
                                    degradation = 11L, salvage = 1L)),
      "score membership: de_novo 1+1, degradation 11, salvage 1")

# EC union per module over core members (composition table, ID-exact)
core_ids <- list(de_novo = cores$pathway_id[cores$final_module == "de_novo"],
                 degradation = cores$pathway_id[cores$final_module == "degradation"],
                 salvage = cores$pathway_id[cores$final_module == "salvage"])
mod_ec <- lapply(core_ids, function(ids)
  unique(unlist(strsplit(comp$ecs[match(ids, comp$pathway_id)], "\\s+"))))
mod_ec <- lapply(mod_ec, function(x) x[!is.na(x) & x != ""])
ec_present <- rownames(E)
ec_absent <- lapply(mod_ec, setdiff, ec_present)
mod_ec <- lapply(mod_ec, intersect, ec_present)
write.csv(do.call(rbind, lapply(names(mod_ec), function(m)
  data.frame(module = m, EC = mod_ec[[m]]))),
  file.path(RES, "05_ec_union_membership.csv"), row.names = FALSE)
cat("EC union sizes:", paste(names(mod_ec), lengths(mod_ec), sep = "=", collapse = "; "),
    " | absent from EC table:", paste(names(ec_absent), lengths(ec_absent), sep = "=", collapse = "; "), "\n")
check("P1e_ec_union_nonempty", all(lengths(mod_ec) > 0), "every module EC union nonempty")

# ---- module scores: DUAL IMPLEMENTATION (P2a) ------------------------------
# impl 1: vectorised colSums; impl 2: explicit per-sample double loop
sum_rows <- function(M, ids) if (length(ids) == 1) M[ids, ] else colSums(M[ids, , drop = FALSE])
sum_rows_loop <- function(M, ids) {
  out <- numeric(ncol(M))
  for (j in seq_len(ncol(M))) for (i in ids) out[j] <- out[j] + M[i, j]
  out
}
scores_abs <- list()
for (nm in names(mod_pwy)) scores_abs[[paste0("P_", nm)]] <- sum_rows(P, mod_pwy[[nm]])
for (nm in names(mod_ec))  scores_abs[[paste0("E_", nm)]] <- sum_rows(E, mod_ec[[nm]])
for (nm in names(mod_pwy)) {
  a <- scores_abs[[paste0("P_", nm)]]; b <- sum_rows_loop(P, mod_pwy[[nm]])
  check(paste0("P2a_dual_", nm), max(abs(a - b)) < 1e-10, sprintf("max diff %.3e", max(abs(a - b))))
}
for (nm in names(mod_ec)) {
  a <- scores_abs[[paste0("E_", nm)]]; b <- sum_rows_loop(E, mod_ec[[nm]])
  check(paste0("P2a_dual_EC_", nm), max(abs(a - b)) < 1e-10, sprintf("max diff %.3e", max(abs(a - b))))
}
tot_P <- colSums(P); tot_E <- colSums(E)
scores <- list()
for (nm in names(scores_abs)) {
  den <- if (startsWith(nm, "P_")) tot_P else tot_E
  scores[[paste0(nm, "_share")]] <- 100 * scores_abs[[nm]] / den   # percent
  scores[[paste0(nm, "_abs")]]   <- scores_abs[[nm]]
}
S <- as.data.frame(scores)
check("P2c_share_range", all(S$P_degradation_share >= 0 & S$P_degradation_share <= 100) &&
        all(S$E_degradation_share >= 0 & S$E_degradation_share <= 100), "shares within [0,100]")

adf <- cbind(data.frame(sample = rownames(sdat), Enterotype = sdat$Enterotype,
                        Gender = sdat$Gender, Age = sdat$Age, UA = sdat$UA,
                        Hyperuricemia = sdat$Hyperuricemia,
                        total_pathway = tot_P, total_ec = tot_E), S)
write.csv(adf, file.path(RES, "02_module_scores_by_sample.csv"), row.names = FALSE)

# ---- Q2: GLS per score variant ---------------------------------------------
gls_one <- function(y, pseudo_mult = 1, label) {
  pseudo <- min(y[y > 0], na.rm = TRUE) / 2 * pseudo_mult
  d <- data.frame(ly = log10(y + pseudo), Enterotype = sdat$Enterotype,
                  Gender = sdat$Gender, Age = sdat$Age)
  fit <- gls(ly ~ Enterotype + Gender + Age, data = d,
             weights = varIdent(form = ~ 1 | Enterotype))
  overall_p <- anova(fit, type = "marginal")["Enterotype", "p-value"]
  # dual implementation of the 6 pairwise contrasts: emmeans vs manual
  em <- emmeans(fit, ~ Enterotype, data = d)
  ce <- as.data.frame(pairs(em, adjust = "holm"))
  beta <- coef(fit); V <- vcov(fit)
  tix <- grep("^Enterotype", names(beta))
  Z <- rbind(0, diag(3))                       # 3 betas -> 4 EMMs
  C <- combn(4, 2, function(ix) { v <- numeric(4); v[ix] <- c(1, -1); v })  # 4x6
  M <- t(C) %*% Z                              # 6x3 contrast map on betas
  est_m <- as.numeric(M %*% beta[tix])
  se_m  <- sqrt(diag(M %*% V[tix, tix] %*% t(M)))
  lv <- levels(d$Enterotype)
  pairs_nm <- apply(combn(4, 2), 2, function(ix) paste(lv[ix], collapse = " - "))
  ord <- match(pairs_nm, ce$contrast)
  check(paste0("P2d_gls_dual_", label),
        max(abs(ce$estimate[ord] - est_m)) < 1e-8 && max(abs(ce$SE[ord] - se_m)) < 1e-6,
        sprintf("emmeans vs manual: max d_est %.2e, d_SE %.2e",
                max(abs(ce$estimate[ord] - est_m)), max(abs(ce$SE[ord] - se_m))))
  check(paste0("P2c_prange_", label), overall_p >= 0 && overall_p <= 1, "P in [0,1]")
  data.frame(score = label, overall_enterotype_p = overall_p, pseudocount = pseudo,
             contrast = ce$contrast, estimate = ce$estimate, SE = ce$SE,
             t_ratio = ce$`t.ratio`, p_holm = ce$p.value, row.names = NULL)
}

variants <- grep("_share$", names(S), value = TRUE)             # primary: shares
variants <- c(variants, grep("_abs$", names(S), value = TRUE))  # sensitivity: absolute
gls_all <- do.call(rbind, lapply(variants, function(v) gls_one(S[[v]], 1, v)))
write.csv(gls_all, file.path(RES, "03_q2_gls_contrasts.csv"), row.names = FALSE)

# pseudo-count sensitivity on the primary (share) variants: pseudo x0.5 / x2
sens <- do.call(rbind, lapply(grep("_share$", names(S), value = TRUE), function(v) {
  r5 <- gls_one(S[[v]], 0.5, paste0(v, " x0.5")); r2 <- gls_one(S[[v]], 2, paste0(v, " x2"))
  data.frame(score = v, pseudo_mult = c(0.5, 2),
             overall_enterotype_p = c(r5$overall_enterotype_p[1], r2$overall_enterotype_p[1]),
             min_p_contrast = c(min(r5$p_holm), min(r2$p_holm)))
}))
write.csv(sens, file.path(RES, "06_pseudo_sensitivity.csv"), row.names = FALSE)

# ---- Q3: partial Spearman vs UA --------------------------------------------
cor_rows <- lapply(variants, function(v) {
  r <- pcor.test(S[[v]], sdat$UA,
                 data.frame(male = as.numeric(sdat$Gender == "Male"), age = sdat$Age),
                 method = "spearman")
  data.frame(score = v, rho = unname(r$estimate), p = unname(r$p.value), n = r$n)
})
cor_df <- do.call(rbind, cor_rows)
cor_df$padj <- p.adjust(cor_df$p, method = "BH")
check("P2c_pcor_n", all(cor_df$n == 305), "pcor n=305 for all variants")
write.csv(cor_df, file.path(RES, "04_q3_module_UA_pcor.csv"), row.names = FALSE)

# ---- Q1: curated extraction from existing LinDA table -----------------------
linda$pathway_id <- sub(":.*$", "", linda$feature)
q1 <- linda[linda$pathway_id %in% inc$pathway_id, ]
# P3 fidelity: re-merge by feature+comparison and confirm identical values
m <- merge(q1, linda, by = c("feature", "comparison"), suffixes = c("", ".src"))
check("P3_linda_fidelity", nrow(m) == nrow(q1) &&
        all(m$log2FoldChange == m$log2FoldChange.src) && all(m$padj == m$padj.src),
      "extracted log2FC/padj identical to LinDA source table")
tested <- unique(q1$pathway_id)
untested <- setdiff(inc$pathway_id, tested)
cat("Q1: ", length(tested), "/", nrow(inc), " frozen pathways tested in LinDA; untested: ",
    paste(untested, collapse = ", "), "\n")
q1_out <- merge(inc[, c("pathway_id", "feature", "final_module", "final_sum_role",
                        "in_pathway_module_score")],
                q1[, c("pathway_id", "comparison", "log2FoldChange", "pvalue", "padj", "sig")],
                by = "pathway_id", all.x = TRUE)
q1_out$tested_in_linda <- !is.na(q1_out$comparison)
write.csv(q1_out, file.path(RES, "01_q1_pathway_linda_frozen34.csv"), row.names = FALSE)

# ---- P4: provenance manifest -------------------------------------------------
inputs <- c(here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"),
            here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"),
            here("data", "clean", "pathabundance_stratified.tsv"),
            here("data", "curated", "purine", "frozen_inclusion_register.csv"),
            here("data", "curated", "purine", "pathway_composition.csv"),
            here("results", "tables", "section_3_4", "linda_pathway_all_comparisons.csv"),
            here("scripts", "figure4", "purine_module_analysis.R"))
man <- data.frame(file = inputs, md5 = tools::md5sum(inputs), row.names = NULL)
write.csv(man, file.path(VER, "05_run_manifest_md5.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(VER, "session_info.txt"))

write.csv(qc, file.path(VER, "06_in_run_qc_assertions.csv"), row.names = FALSE)
cat("\nAll", nrow(qc), "QC assertions passed.\n")
sink(); close(con)

# status gate: DRAFT until the independent verification script flips it
writeLines(c("Purine module analysis", paste("main run:", Sys.time()),
             sprintf("QC assertions passed: %d/%d", nrow(qc), nrow(qc)),
             "STATUS: DRAFT (awaiting 99_independent_verification.R)"),
           file.path(VER, "v05_STATUS_DRAFT.txt"))
cat("DONE. Results in", RES, "\n")
