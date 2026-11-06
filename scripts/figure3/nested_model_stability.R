# Nested-model comparison and Firth stability analysis
#
# Refits the prespecified clinical models in fixed complete-case populations to
# separate sample-missingness and covariate-adjustment contributions. Reports
# events per variable and Firth-penalized HUA models with QC assertions.

suppressPackageStartupMessages({
  library(nlme)
  library(logistf)
  library(here)
})

RES <- here("results", "tables", "section_3_3", "nested_model_stability")
VER <- file.path(RES, "verification")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(VER, showWarnings = FALSE, recursive = TRUE)

runlog <- file.path(VER, "run_log.txt")
con <- file(runlog, open = "wt"); sink(con, split = TRUE)
cat("Nested-model and Firth stability analysis run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# ---- QC assertion harness --------------------------------------------------
qc <- data.frame(id = character(), pass = logical(), detail = character(),
                 stringsAsFactors = FALSE)
check <- function(id, pass, detail = "") {
  qc <<- rbind(qc, data.frame(id = id, pass = isTRUE(pass), detail = detail))
  cat(sprintf("[QC] %-34s %s  %s\n", id, ifelse(isTRUE(pass), "PASS", "FAIL"), detail))
  if (!isTRUE(pass)) {
    sink(); close(con)
    stop("QC assertion failed: ", id, " — ", detail, call. = FALSE)
  }
}

# ---- inputs ----------------------------------------------------------------
DATA_CSV <- here("results", "tables", "section_3_3", "section3_6_analysis_data.csv")
PUB_UA   <- here("results", "tables", "section_3_3", "ua_model_comparison.csv")
PUB_HUA  <- here("results", "tables", "section_3_3", "hua_model_comparison.csv")

raw <- read.csv(DATA_CSV, stringsAsFactors = FALSE)
raw$Enterotype    <- factor(raw$Enterotype, levels = paste("Type", 1:4))
raw$Gender        <- factor(raw$Gender, levels = c("Female", "Male"))
raw$Hyperuricemia <- factor(raw$Hyperuricemia, levels = c("No", "Yes"))

# ---- P1: model register (explicit, machine-readable) ------------------------
model_specs <- list(
  M1_Demographic = list(
    ua  = UA ~ Enterotype + Gender + Age,
    hua = Hyperuricemia ~ Enterotype + Gender + Age,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age"),
    npred = 5L),
  M2_BMI = list(
    ua  = UA ~ Enterotype + Gender + Age + BMI,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + BMI,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI"),
    npred = 6L),
  M3_Renal_Cr = list(
    ua  = UA ~ Enterotype + Gender + Age + Cr + Urea,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + Cr + Urea,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea"),
    npred = 7L),
  M4_Cardiomet = list(
    ua  = UA ~ Enterotype + Gender + Age + Cr + Urea + GLU + TG + HDLC + LDLC,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + Cr + Urea + GLU + TG + HDLC + LDLC,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC"),
    npred = 10L),
  S1_Renal_eGFR = list(
    ua  = UA ~ Enterotype + Gender + Age + eGFR + Urea,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + eGFR + Urea,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "eGFR", "Urea"),
    npred = 7L),
  S2_Full_Clinical = list(
    ua  = UA ~ Enterotype + Gender + Age + BMI + Cr + Urea + GLU + TG + HDLC + LDLC + CRP,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + BMI + Cr + Urea + GLU + TG + HDLC + LDLC + CRP,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC", "CRP"),
    npred = 13L)
)
register <- data.frame(
  model = names(model_specs),
  ua_formula  = vapply(model_specs, function(s) paste(deparse(s$ua), collapse = ""), character(1)),
  hua_formula = vapply(model_specs, function(s) paste(deparse(s$hua), collapse = ""), character(1)),
  n_predictors = vapply(model_specs, function(s) s$npred, integer(1)),
  stringsAsFactors = FALSE)

check("P1_register_6_models", nrow(register) == 6, "6 models registered")
check("P1_predictor_counts", all(register$n_predictors == c(5, 6, 7, 10, 7, 13)),
      "non-intercept predictor counts 5/6/7/10/7/13")

# ---- structural checks on the analysis dataset ------------------------------
check("P2b_n_rows", nrow(raw) == 305, sprintf("n = %d", nrow(raw)))
check("P2b_enterotype_counts",
      identical(as.integer(table(raw$Enterotype)), c(120L, 98L, 48L, 39L)),
      paste(as.integer(table(raw$Enterotype)), collapse = "/"))
check("P2b_hua_events", sum(raw$Hyperuricemia == "Yes") == 64, "64 HUA events")
check("P2b_no_missing_core",
      all(complete.cases(raw[, c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age",
                                 "Cr", "Urea", "eGFR", "GLU", "TG", "HDLC", "LDLC")])),
      "core vars complete for all 305")

cc_bmi  <- complete.cases(raw[, model_specs$M2_BMI$vars])
cc_full <- complete.cases(raw[, model_specs$S2_Full_Clinical$vars])
check("P2b_cc_bmi225",  sum(cc_bmi) == 225,  sprintf("BMI-complete n = %d", sum(cc_bmi)))
check("P2b_cc_full223", sum(cc_full) == 223, sprintf("full-clinical-complete n = %d", sum(cc_full)))
check("P2b_et_counts_223",
      identical(as.integer(table(raw$Enterotype[cc_full])), c(88L, 77L, 36L, 22L)),
      paste(as.integer(table(raw$Enterotype[cc_full])), collapse = "/"))

# ---- model fitting engines --------------------------------------------------
fit_ua <- function(form, df) {
  gls(form, data = df, weights = varIdent(form = ~ 1 | Enterotype),
      method = "REML", na.action = na.omit)
}
ua_overall_p <- function(fit) {
  anova(fit, type = "marginal")["Enterotype", "p-value"]
}
# dual: overall P must be invariant to reference-level choice. Tolerance 1e-5
# accommodates GLS REML optimizer noise, which grows with model complexity
# (max observed 1.4e-6 at S2; agreement to >= 6 decimals, far below the
# 3-decimal reporting precision).
ua_overall_p_relevel <- function(form, df) {
  df2 <- df; df2$Enterotype <- relevel(df2$Enterotype, ref = "Type 3")
  anova(fit_ua(form, df2), type = "marginal")["Enterotype", "p-value"]
}

fit_hua <- function(form, df) glm(form, data = df, family = binomial)
hua_overall_p <- function(form, df) {
  fit  <- fit_hua(form, df)
  null <- glm(update(form, . ~ . - Enterotype), data = df, family = binomial)
  anova(null, fit, test = "Chisq")$`Pr(>Chi)`[2]
}
# dual: manual LRT from log-likelihoods (base R path)
hua_overall_p_manual <- function(form, df) {
  fit  <- fit_hua(form, df)
  null <- glm(update(form, . ~ . - Enterotype), data = df, family = binomial)
  lr <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(null)))
  pchisq(lr, df = 3, lower.tail = FALSE)
}

extract_et_terms <- function(est, se) {
  data.frame(term = c("Type2_vs_Type1", "Type3_vs_Type1", "Type4_vs_Type1"),
             estimate = est, se = se, stringsAsFactors = FALSE)
}

ua_et_estimates <- function(fit) {
  tt <- summary(fit)$tTable
  i <- match(paste0("EnterotypeType ", 2:4), rownames(tt))
  extract_et_terms(tt[i, "Value"], tt[i, "Std.Error"])
}
hua_et_estimates <- function(fit) {
  cf <- summary(fit)$coefficients
  i <- match(paste0("EnterotypeType ", 2:4), rownames(cf))
  data.frame(term = c("Type2_vs_Type1", "Type3_vs_Type1", "Type4_vs_Type1"),
             OR = exp(cf[i, "Estimate"]),
             ci_low = exp(cf[i, "Estimate"] - 1.96 * cf[i, "Std. Error"]),
             ci_high = exp(cf[i, "Estimate"] + 1.96 * cf[i, "Std. Error"]),
             p = cf[i, "Pr(>|z|)"], stringsAsFactors = FALSE)
}

run_one_model <- function(spec, df, model_label, panel_label) {
  fu <- fit_ua(spec$ua, df)
  fh <- fit_hua(spec$hua, df)
  list(
    summary = data.frame(
      panel = panel_label, model = model_label, n = nrow(df),
      ua_enterotype_p = ua_overall_p(fu),
      ua_enterotype_p_dual = ua_overall_p_relevel(spec$ua, df),
      hua_enterotype_p = hua_overall_p(spec$hua, df),
      hua_enterotype_p_dual = hua_overall_p_manual(spec$hua, df),
      hua_events = sum(df$Hyperuricemia == "Yes"),
      stringsAsFactors = FALSE),
    ua_terms = cbind(panel = panel_label, model = model_label, ua_et_estimates(fu)),
    hua_terms = cbind(panel = panel_label, model = model_label, hua_et_estimates(fh))
  )
}

# ============================== PART A: fidelity rerun =======================
cat("---- Part A: fidelity rerun of the six published models ----\n")
pub_ua  <- read.csv(PUB_UA, stringsAsFactors = FALSE)
pub_hua <- read.csv(PUB_HUA, stringsAsFactors = FALSE)

partA <- lapply(names(model_specs), function(mn) {
  spec <- model_specs[[mn]]
  df <- raw[complete.cases(raw[, spec$vars, drop = FALSE]), ]
  run_one_model(spec, df, mn, "original_n")
})
partA_sum <- do.call(rbind, lapply(partA, `[[`, "summary"))

tol <- 1e-10
pub_ua_p  <- tapply(pub_ua$overall_enterotype_p, pub_ua$model, unique)[partA_sum$model]
pub_hua_p <- tapply(pub_hua$overall_enterotype_p, pub_hua$model, unique)[partA_sum$model]
check("P3_fidelity_ua_overall_p",
      all(abs(partA_sum$ua_enterotype_p - as.numeric(pub_ua_p)) < tol),
      paste(sprintf("%s: %.3e vs %.3e", partA_sum$model, partA_sum$ua_enterotype_p,
                    as.numeric(pub_ua_p)), collapse = " | "))
check("P3_fidelity_hua_overall_p",
      all(abs(partA_sum$hua_enterotype_p - as.numeric(pub_hua_p)) < tol),
      paste(sprintf("%s: %.3e vs %.3e", partA_sum$model, partA_sum$hua_enterotype_p,
                    as.numeric(pub_hua_p)), collapse = " | "))
pub_n <- tapply(pub_ua$n, pub_ua$model, unique)[partA_sum$model]
check("P3_fidelity_n", all(partA_sum$n == as.numeric(pub_n)),
      paste(partA_sum$model, partA_sum$n, collapse = "; "))
pub_ua_t4 <- pub_ua[pub_ua$label == "Type 4 vs Type 1", ]
pub_ua_t4 <- pub_ua_t4[match(partA_sum$model, pub_ua_t4$model), ]
partA_ua_terms <- do.call(rbind, lapply(partA, `[[`, "ua_terms"))
t4_rows <- partA_ua_terms[partA_ua_terms$term == "Type4_vs_Type1", ]
t4_rows <- t4_rows[match(partA_sum$model, t4_rows$model), ]
check("P3_fidelity_ua_t4_beta",
      all(abs(t4_rows$estimate - pub_ua_t4$estimate) < tol),
      paste(sprintf("%s: %.4f vs %.4f", t4_rows$model, t4_rows$estimate,
                    pub_ua_t4$estimate), collapse = " | "))

check("P2a_dual_ua_overall_p",
      all(abs(partA_sum$ua_enterotype_p - partA_sum$ua_enterotype_p_dual) < 1e-5),
      sprintf("relevel-invariance of marginal Wald P; max diff = %.2e",
              max(abs(partA_sum$ua_enterotype_p - partA_sum$ua_enterotype_p_dual))))
check("P2a_dual_hua_overall_p",
      all(abs(partA_sum$hua_enterotype_p - partA_sum$hua_enterotype_p_dual) < 1e-12),
      "anova LRT vs manual logLik LRT")

# ============================== PART B: nested models ========================
cat("\n---- Part B: nested models within fixed complete-case populations ----\n")
df_bmi  <- raw[cc_bmi, , drop = FALSE]
df_full <- raw[cc_full, , drop = FALSE]

panel_bmi <- lapply(c("M1_Demographic", "M2_BMI"), function(mn)
  run_one_model(model_specs[[mn]], df_bmi, mn, "BMI225"))
panel_full <- lapply(c("M1_Demographic", "M2_BMI", "M3_Renal_Cr", "M4_Cardiomet",
                       "S2_Full_Clinical"), function(mn)
  run_one_model(model_specs[[mn]], df_full, mn, "FULL223"))

partB <- c(panel_bmi, panel_full)
partB_sum <- do.call(rbind, lapply(partB, `[[`, "summary"))
partB_ua  <- do.call(rbind, lapply(partB, `[[`, "ua_terms"))
partB_hua <- do.call(rbind, lapply(partB, `[[`, "hua_terms"))

check("P2a_dualB_ua_overall_p",
      all(abs(partB_sum$ua_enterotype_p - partB_sum$ua_enterotype_p_dual) < 1e-5),
      sprintf("relevel-invariance in nested panels; max diff = %.2e",
              max(abs(partB_sum$ua_enterotype_p - partB_sum$ua_enterotype_p_dual))))
check("P2a_dualB_hua_overall_p",
      all(abs(partB_sum$hua_enterotype_p - partB_sum$hua_enterotype_p_dual) < 1e-12),
      "manual LRT in nested panels")
check("P2c_p_range",
      all(partB_sum$ua_enterotype_p >= 0 & partB_sum$ua_enterotype_p <= 1 &
          partB_sum$hua_enterotype_p >= 0 & partB_sum$hua_enterotype_p <= 1),
      "all P in [0,1]")
check("P2b_panel_n",
      all(partB_sum$n[partB_sum$panel == "BMI225"] == 225) &&
      all(partB_sum$n[partB_sum$panel == "FULL223"] == 223),
      "panel sizes fixed at 225 / 223")

# ============================== PART C: stability ============================
cat("\n---- Part C: EPV and Firth penalized logistic models ----\n")
plrt_enterotype <- function(form, df) {
  ff <- logistf(form, data = df)
  f0 <- logistf(update(form, . ~ . - Enterotype), data = df)
  stat <- 2 * (unname(ff$loglik["full"]) - unname(f0$loglik["full"]))
  pchisq(stat, df = 3, lower.tail = FALSE)
}
plrt_enterotype_relevel <- function(form, df) {
  df2 <- df; df2$Enterotype <- relevel(df2$Enterotype, ref = "Type 3")
  plrt_enterotype(form, df2)
}
epv_rows <- list()
firth_rows <- list()
for (mn in names(model_specs)) {
  spec <- model_specs[[mn]]
  df <- raw[complete.cases(raw[, spec$vars, drop = FALSE]), ]
  events <- sum(df$Hyperuricemia == "Yes")
  epv_rows[[mn]] <- data.frame(model = mn, n = nrow(df), events = events,
                               predictors = spec$npred, epv = events / spec$npred,
                               stringsAsFactors = FALSE)
  ff <- logistf(spec$hua, data = df)
  cf <- ff$coefficients
  i <- match(paste0("EnterotypeType ", 2:4), names(cf))
  ci <- confint(ff)
  firth_rows[[mn]] <- data.frame(
    model = mn, n = nrow(df), events = events,
    term = c("Type2_vs_Type1", "Type3_vs_Type1", "Type4_vs_Type1"),
    OR = exp(cf[i]), ci_low = exp(ci[i, 1]), ci_high = exp(ci[i, 2]),
    enterotype_plrt_p = plrt_enterotype(spec$hua, df),
    enterotype_plrt_p_dual = plrt_enterotype_relevel(spec$hua, df),
    enterotype_glm_lrt_p = hua_overall_p(spec$hua, df),
    stringsAsFactors = FALSE)
}
epv_df <- do.call(rbind, epv_rows)
firth_df <- do.call(rbind, firth_rows)

check("P2a_dual_firth_plrt",
      all(abs(firth_df$enterotype_plrt_p - firth_df$enterotype_plrt_p_dual) < 1e-10),
      sprintf("penalized LRT relevel-invariance; max diff = %.2e",
              max(abs(firth_df$enterotype_plrt_p - firth_df$enterotype_plrt_p_dual))))
check("P2c_epv_positive", all(epv_df$epv > 0), "EPV computed for all models")
check("P2b_events_sum",
      sum(epv_df$events[epv_df$model == "M1_Demographic"]) == 64, "M1 events = 64")
check("P2c_firth_or_range",
      all(firth_df$OR > 0 & firth_df$ci_low > 0 &
          firth_df$ci_low <= firth_df$OR & firth_df$OR <= firth_df$ci_high),
      "Firth ORs positive and within profile CIs")

# ============================== write outputs ================================
write.csv(register, file.path(RES, "00_model_register.csv"), row.names = FALSE)
write.csv(partA_sum, file.path(RES, "01_fidelity_rerun_summary.csv"), row.names = FALSE)
write.csv(partB_sum, file.path(RES, "02_nested_models_summary.csv"), row.names = FALSE)
write.csv(partB_ua, file.path(RES, "03_nested_models_ua_contrasts.csv"), row.names = FALSE)
write.csv(partB_hua, file.path(RES, "04_nested_models_hua_contrasts.csv"), row.names = FALSE)
write.csv(epv_df, file.path(RES, "05_epv_by_model.csv"), row.names = FALSE)
write.csv(firth_df, file.path(RES, "06_firth_hua_models.csv"), row.names = FALSE)
write.csv(qc, file.path(VER, "06_in_run_qc_assertions.csv"), row.names = FALSE)

# ============================== P4: provenance ===============================
md5 <- tools::md5sum(c(DATA_CSV, PUB_UA, PUB_HUA, here("scripts", "figure3", "nested_model_stability.R")))
manifest <- data.frame(file = names(md5), md5 = unname(md5), stringsAsFactors = FALSE)
write.csv(manifest, file.path(VER, "05_run_manifest_md5.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(VER, "session_info.txt"))
writeLines(paste("STATUS: DRAFT —", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
           file.path(VER, "v05_STATUS_DRAFT.txt"))

cat("\n---- Part B headline numbers ----\n")
print(partB_sum[, c("panel", "model", "n", "ua_enterotype_p", "hua_enterotype_p", "hua_events")])
cat("\n---- EPV ----\n")
print(epv_df)
cat("\n---- Firth Type 4 vs Type 1 ----\n")
print(firth_df[firth_df$term == "Type4_vs_Type1",
               c("model", "n", "OR", "ci_low", "ci_high", "enterotype_plrt_p")])
cat("\nAll QC assertions passed:", nrow(qc), "\n")
sink(); close(con)
