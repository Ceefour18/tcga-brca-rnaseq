################################################################################
# 05_survival.R
#
# Project : TCGA-BRCA RNA-seq — Basal-like vs Luminal A
# Purpose : Prognostic analysis of subtype, run IDENTICALLY on two endpoints so
#           they can be compared on equal footing:
#             PFI (primary)  and  OS (secondary)
#           For EACH endpoint: Kaplan-Meier + log-rank, unadjusted and
#           age-adjusted Cox, proportional-hazards check, forest plot.
#           Then a side-by-side comparison + a data-derived summary.
#
# Author  : Clement Akinsola (Ceefor Analytic Hub)
#
#
# Design: ADJUST_FOR = c("age") applied to both endpoints (parallel). Age is
# near-complete; stage (missingness), gender (near-constant) and ER/PR/HER2
# (collinear with subtype) are excluded.
#
# Input   : data/processed/dds_basal_vs_luma.rds
# Outputs : results/figures/km_{pfi,os}_basal_vs_luma.png
#           results/figures/cox_forest_{pfi,os}.png
#           results/tables/cox_all_models.csv
#           results/tables/subtype_hr_comparison.csv
#           results/tables/survival_endpoint_summary.csv
#           results/tables/cox_ph_{pfi,os}.txt
#           results/tables/survival_summary.txt          (data-derived)
#           data/processed/survival_results.rds
#           results/tables/survival_session_info.txt
################################################################################


# ---- 0. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2); library(survival); library(survminer)
  library(dplyr); library(here)
})


# ---- 1. Parameters & paths --------------------------------------------------
ADJUST_FOR          <- c("age")
FOLLOWUP_CAP_MONTHS <- NULL
DAYS_PER_MONTH      <- 30.44
PALETTE             <- c("#0072B2", "#D55E00")   # LumA blue, Basal orange

dir_proc    <- here("data", "processed")
dir_figures <- here("results", "figures")
dir_tables  <- here("results", "tables")
for (d in c(dir_figures, dir_tables)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ---- Helpers ----------------------------------------------------------------
median_followup <- function(time, event) {
  ok <- !is.na(time) & !is.na(event)
  if (sum(ok) == 0) return(NA_real_)
  unname(summary(survfit(Surv(time[ok], 1 - event[ok]) ~ 1))$table["median"])
}
condition_stats <- function(fit) {
  s <- summary(fit); i <- grep("^condition", rownames(s$coefficients))[1]
  list(HR = unname(s$conf.int[i, "exp(coef)"]),
       lo = unname(s$conf.int[i, "lower .95"]),
       hi = unname(s$conf.int[i, "upper .95"]),
       p  = unname(s$coefficients[i, "Pr(>|z|)"]))
}
tidy_cox <- function(fit, label) {
  s <- summary(fit)
  data.frame(model = label, term = rownames(s$coefficients),
             HR = round(s$conf.int[, "exp(coef)"], 3),
             CI_low = round(s$conf.int[, "lower .95"], 3),
             CI_high = round(s$conf.int[, "upper .95"], 3),
             p_value = signif(s$coefficients[, "Pr(>|z|)"], 3), row.names = NULL)
}
apply_cap <- function(time, event, cap) {
  if (is.null(cap)) return(data.frame(time = time, event = event))
  data.frame(time = pmin(time, cap), event = ifelse(time > cap, 0L, event))
}


# ---- 2. Load cohort ---------------------------------------------------------
# Cohort = the exact DE cohort frozen in the dds. Times converted days -> months.
dds <- readRDS(file.path(dir_proc, "dds_basal_vs_luma.rds"))
surv_df <- as.data.frame(colData(dds)) %>%
  transmute(
    condition,
    PFI = as.integer(PFI), PFI.time = as.numeric(PFI.time) / DAYS_PER_MONTH,
    OS  = as.integer(OS),  OS.time  = as.numeric(OS.time)  / DAYS_PER_MONTH,
    age = as.numeric(age)
  )
surv_df$condition <- factor(surv_df$condition, levels = c("LumA", "Basal"),
                            labels = c("Luminal A", "Basal-like"))
message("Cohort: ", nrow(surv_df), " tumours (",
        paste(names(table(surv_df$condition)), table(surv_df$condition),
              sep = "=", collapse = ", "), ")")


# ---- 3. One analysis, applied to each endpoint ------------------------------
# Identical KM + Cox + PH + forest for whichever endpoint is passed in.
analyse_endpoint <- function(df, event_col, time_col, label, ylab,
                             adjust_for, cap, km_path, forest_path, ph_path) {
  d <- df
  d$event <- as.integer(df[[event_col]]); d$time <- as.numeric(df[[time_col]])
  d <- d[!is.na(d$event) & !is.na(d$time), , drop = FALSE]
  cp <- apply_cap(d$time, d$event, cap); d$time <- cp$time; d$event <- cp$event

  # KM + log-rank
  fit_km <- survfit(Surv(time, event) ~ condition, data = d)
  lr  <- survdiff(Surv(time, event) ~ condition, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  # Cox (unadjusted + adjusted) on the same complete-case set
  covar_ok <- if (is.null(adjust_for)) rep(TRUE, nrow(d)) else
              complete.cases(d[, adjust_for, drop = FALSE])
  cd <- d[covar_ok, ]
  bf <- function(av) {
    t <- c("condition", av); t <- t[!is.null(t) & nzchar(t)]
    as.formula(paste("Surv(time, event) ~", paste(t, collapse = " + ")))
  }
  cox_u <- coxph(bf(NULL), data = cd)
  cox_a <- coxph(bf(adjust_for), data = cd)
  ph <- cox.zph(cox_a); ph_p <- unname(ph$table["GLOBAL", "p"])
  writeLines(capture.output(print(ph)), ph_path)

  # KM plot
  km <- ggsurvplot(fit_km, data = d, pval = TRUE, conf.int = TRUE,
                   risk.table = TRUE, palette = PALETTE,
                   legend.labs = c("Luminal A", "Basal-like"),
                   legend.title = "PAM50", xlab = "Time (months)", ylab = ylab,
                   title = paste0(label, " by subtype — TCGA-BRCA"),
                   break.time.by = 24, risk.table.height = 0.28,
                   ggtheme = theme_bw())
  png(km_path, width = 8, height = 7, units = "in", res = 300); print(km); dev.off()

  # Forest (adjusted model)
  fp <- ggforest(cox_a, data = cd,
                 main = paste0("Cox ", label, " (age-adjusted) — Basal vs LumA"))
  ggsave(forest_path, fp, width = 8, height = 4, dpi = 300)

  rows <- rbind(
    cbind(endpoint = label, tidy_cox(cox_u, "Unadjusted (subtype only)")),
    cbind(endpoint = label, tidy_cox(cox_a,
          paste0("Adjusted (+ ", paste(adjust_for, collapse = ", "), ")")))
  )
  list(label = label, fit_km = fit_km, lr_p = lr_p, cox_u = cox_u, cox_a = cox_a,
       ph = ph, ph_p = ph_p, rows = rows,
       cs_u = condition_stats(cox_u), cs_a = condition_stats(cox_a),
       n = nrow(d), events = sum(d$event == 1), n_cox = nrow(cd),
       median_fu = median_followup(d$time, d$event))
}

message("\nAnalysing PFI ...")
res_pfi <- analyse_endpoint(
  surv_df, "PFI", "PFI.time", "PFI", "Progression-free interval probability",
  ADJUST_FOR, FOLLOWUP_CAP_MONTHS,
  file.path(dir_figures, "km_pfi_basal_vs_luma.png"),
  file.path(dir_figures, "cox_forest_pfi.png"),
  file.path(dir_tables, "cox_ph_pfi.txt"))

message("Analysing OS ...")
res_os <- analyse_endpoint(
  surv_df, "OS", "OS.time", "OS", "Overall survival probability",
  ADJUST_FOR, FOLLOWUP_CAP_MONTHS,
  file.path(dir_figures, "km_os_basal_vs_luma.png"),
  file.path(dir_figures, "cox_forest_os.png"),
  file.path(dir_tables, "cox_ph_os.txt"))


# ---- 4. Combined tables -----------------------------------------------------
cox_all <- rbind(res_pfi$rows, res_os$rows)
write.csv(cox_all, file.path(dir_tables, "cox_all_models.csv"), row.names = FALSE)
print(cox_all)

# Subtype (Basal vs LumA) effect across endpoints x adjustment, with log-rank
# and PH verdict alongside — the core comparison.
mk <- function(r) data.frame(
  endpoint = r$label,
  model = c("Unadjusted", "Age-adjusted"),
  HR = round(c(r$cs_u$HR, r$cs_a$HR), 3),
  CI_low = round(c(r$cs_u$lo, r$cs_a$lo), 3),
  CI_high = round(c(r$cs_u$hi, r$cs_a$hi), 3),
  CI_width = round(c(r$cs_u$hi - r$cs_u$lo, r$cs_a$hi - r$cs_a$lo), 3),
  cox_p = signif(c(r$cs_u$p, r$cs_a$p), 3),
  logrank_p = signif(r$lr_p, 3),
  PH_global_p = signif(c(NA, r$ph_p), 3))
subtype_cmp <- rbind(mk(res_pfi), mk(res_os))
write.csv(subtype_cmp, file.path(dir_tables, "subtype_hr_comparison.csv"),
          row.names = FALSE)
print(subtype_cmp)

endpoint_summary <- data.frame(
  endpoint = c("PFI (primary)", "OS (secondary)"),
  n = c(res_pfi$n, res_os$n),
  events = c(res_pfi$events, res_os$events),
  event_rate = round(c(res_pfi$events / res_pfi$n, res_os$events / res_os$n), 3),
  median_fu_months = round(c(res_pfi$median_fu, res_os$median_fu), 1))
write.csv(endpoint_summary, file.path(dir_tables, "survival_endpoint_summary.csv"),
          row.names = FALSE)
print(endpoint_summary)


# ---- 5. Data-derived comparison summary -------------------------------------
verdict <- function(p) ifelse(p < 0.05, "significant", "not significant")
ph_word <- function(p) ifelse(p < 0.05, "violated", "supported")
tighter <- ifelse(res_pfi$cs_a$hi - res_pfi$cs_a$lo <
                  res_os$cs_a$hi - res_os$cs_a$lo, "PFI", "OS")
agree   <- ifelse((res_pfi$cs_a$HR > 1) == (res_os$cs_a$HR > 1),
                  "agree on direction", "disagree on direction")

summary_lines <- c(
  sprintf("Cohort: %d tumours (%s).", nrow(surv_df),
          paste(names(table(surv_df$condition)), table(surv_df$condition),
                sep = "=", collapse = ", ")),
  sprintf("PFI: %d/%d events (%.1f%%), median FU %.1f mo | log-rank p=%.3g (%s) | adj HR=%.2f (%.2f-%.2f) p=%.3g (%s) | PH p=%.3g (%s).",
          res_pfi$events, res_pfi$n, 100 * res_pfi$events / res_pfi$n,
          res_pfi$median_fu, res_pfi$lr_p, verdict(res_pfi$lr_p),
          res_pfi$cs_a$HR, res_pfi$cs_a$lo, res_pfi$cs_a$hi, res_pfi$cs_a$p,
          verdict(res_pfi$cs_a$p), res_pfi$ph_p, ph_word(res_pfi$ph_p)),
  sprintf("OS : %d/%d events (%.1f%%), median FU %.1f mo | log-rank p=%.3g (%s) | adj HR=%.2f (%.2f-%.2f) p=%.3g (%s) | PH p=%.3g (%s).",
          res_os$events, res_os$n, 100 * res_os$events / res_os$n,
          res_os$median_fu, res_os$lr_p, verdict(res_os$lr_p),
          res_os$cs_a$HR, res_os$cs_a$lo, res_os$cs_a$hi, res_os$cs_a$p,
          verdict(res_os$cs_a$p), res_os$ph_p, ph_word(res_os$ph_p)),
  sprintf("Comparison: endpoints %s; adjusted-HR 95%% CI is tighter for %s (PFI width %.2f vs OS width %.2f).",
          agree, tighter, res_pfi$cs_a$hi - res_pfi$cs_a$lo,
          res_os$cs_a$hi - res_os$cs_a$lo)
)
writeLines(summary_lines, file.path(dir_tables, "survival_summary.txt"))
message("\n----- data-derived endpoint comparison -----")
message(paste(summary_lines, collapse = "\n"))


# ---- 6. Save + session ------------------------------------------------------
saveRDS(list(pfi = res_pfi, os = res_os, subtype_cmp = subtype_cmp,
             endpoint_summary = endpoint_summary, summary_lines = summary_lines),
        file.path(dir_proc, "survival_results.rds"))
writeLines(capture.output(sessionInfo()),
           con = file.path(dir_tables, "survival_session_info.txt"))
message("\nSaved survival_results.rds. Pipeline (01-05) complete.")
