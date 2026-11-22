# ==============================================================================
# Helper, Survival Modelling, and Export Functions
# ==============================================================================

# standardise demographic, DHB region, age group, and ethnicity variables
process_demo <- function(dat,
                         age1_upper = 1, age1name = "under 1y",
                         age2_lower = 1, age2_upper = 5, age2name = "1-4y",
                         age3_lower = 5, age3_upper = 15, age3name = "5-14y",
                         age4_lower = 15, age4_upper = 20, age4name = "15-19y",
                         age5_lower = 20, age5_upper = 30, age5name = "20-29y") {
  
  use_age5 <- !(is.null(age5_lower) || is.null(age5_upper) || is.null(age5name))
  age_levels <- if (use_age5) {
    c(age1name, age2name, age3name, age4name, age5name)
  } else {
    c(age1name, age2name, age3name, age4name)
  }
  
  dat %>%
    mutate(
      dob_firstevstdate = floor(interval(start = date_of_birth, end = firstevstdate) / dyears(1)),
      dhbregion = fct_collapse(
        dhb_2015,
        NLD_AKL_WKO = dhbCodes$three_digit_code[1:5],
        BOP_LKS     = dhbCodes$three_digit_code[6:10],
        MID_WTG     = dhbCodes$three_digit_code[11:15],
        STH_IS      = dhbCodes$three_digit_code[16:20]
      ),
      .after = dhb_2015
    ) %>%
    mutate(
      dhb_2015 = factor(dhb_2015, levels = dhbCodes$three_digit_code, labels = dhbCodes$name),
      age_group = case_when(
        dob_firstevstdate < age1_upper ~ age1name,
        dob_firstevstdate >= age2_lower & dob_firstevstdate < age2_upper ~ age2name,
        dob_firstevstdate >= age3_lower & dob_firstevstdate < age3_upper ~ age3name,
        dob_firstevstdate >= age4_lower & dob_firstevstdate < age4_upper ~ age4name,
        use_age5 & dob_firstevstdate >= age5_lower & dob_firstevstdate < age5_upper ~ age5name,
        TRUE ~ NA_character_
      ),
      age_group = factor(age_group, levels = age_levels),
      # Total response ethnicity indicators across recorded ethnicity fields
      is_maori = if_any(starts_with("ethnic_code_"), ~ !is.na(.x) & .x == "21"),
      is_pacific = if_any(starts_with("ethnic_code_"), ~ !is.na(.x) & str_detect(.x, "^3[0-7]")),
      anyMaori_other = factor(if_else(is_maori, "Maori", "non-Maori"),
                              levels = c("non-Maori", "Maori")),
      anyPacific_other = factor(if_else(is_pacific, "Pacific", "non-Pacific"),
                                levels = c("non-Pacific", "Pacific")),
      anyMP_other = factor(if_else(is_maori | is_pacific, "Maori/Pacific", "non-Maori/non-Pacific"),
                           levels = c("non-Maori/non-Pacific", "Maori/Pacific"))
    ) %>%
    select(-is_maori, -is_pacific)
}

# categorise event dates relative to index measles hospitalisation or equivalent (firstevstdate)
categorise_followup_period <- function(event_date, index_date) {
  diff_days <- as.numeric(difftime(event_date, index_date, units = "days"))
  case_when(
    diff_days < -365                     ~ "pre",
    diff_days >= -365 & diff_days < 0    ~ "base",
    diff_days >= 0    & diff_days <= 30  ~ "infn",
    diff_days > 30                       ~ "post",
    TRUE                                 ~ NA_character_
  )
}

# structure post-acute hospitalisation events into survival intervals (ttf or recurrent)
build_hosp_survival_data <- function(events_df,
                                     cohort_df,
                                     visitnumber = NULL,
                                     tcensor_days,
                                     cutoff_date = as.Date("2018-12-31")) {
  
  evsurv <- events_df %>%
    filter(base_post == "post") %>%
    arrange(cluster, hcu_id, evstdate) %>%
    mutate(status = 1)
  
  cc_subs_event <- evsurv %>%
    mutate(
      censor_date = pmin(cutoff_date, firstevstdate + days(tcensor_days)),
      status = if_else(tstop < censor_date, 1, 0),
      tstop = pmin(tstop, censor_date),
      tstart = if_else(enum == 1, firstevstdate + days(30), tstart),
      lagint = as.numeric(tstop - tstart)
    ) %>%
    select(-censor_date)
  
  # exclude same-day transfers or events outside the valid follow-up window for recurrent models
  if (is.null(visitnumber)) {
    cc_subs_event <- cc_subs_event %>%
      mutate(censor_date = pmin(cutoff_date, firstevstdate + days(tcensor_days))) %>%
      filter(tstop >= tstart, evstdate <= censor_date) %>%
      select(-censor_date)
  }
  
  # add terminal censored interval at end of follow-up for individuals with >=1 event
  final_intervals <- cc_subs_event %>%
    group_by(hcu_id) %>%
    filter(enum == max(enum)) %>%
    mutate(
      event_id = NA_character_,
      evstdate = as.Date(NA),
      evendate = as.Date(NA),
      tstart = tstop,
      tstop = pmin(cutoff_date, firstevstdate + days(tcensor_days)),
      lagint = as.numeric(tstop - tstart),
      status = 0,
      enum = enum + 1
    ) %>%
    ungroup()
  
  cc_subs_event <- bind_rows(cc_subs_event, final_intervals)
  
  # create single censored interval for individuals with zero post-acute events
  cc_null_subs <- cohort_df %>%
    filter(!hcu_id %in% cc_subs_event$hcu_id) %>%
    mutate(
      event_id = NA_character_,
      evstdate = as.Date(NA),
      evendate = as.Date(NA),
      enum = 1,
      tstart = firstevstdate + days(30),
      tstop = pmin(cutoff_date, firstevstdate + days(tcensor_days)),
      lagint = as.numeric(tstop - tstart),
      evint = lagint,
      status = 0
    )
  
  out <- bind_rows(cc_subs_event, cc_null_subs) %>%
    standardise_survival_covariates(date_col = "evstdate")
  
  if (!is.null(visitnumber)) {
    out <- out %>% filter(enum == visitnumber)
  }
  
  out
}

# structure post-acute pharmaceutical dispensings into survival intervals
build_pharm_survival_data <- function(phh_df,
                                      cohort_df,
                                      visitnumber = NULL,
                                      tcensor_days,
                                      cutoff_date = as.Date("2018-12-31")) {
  
  evsurv <- phh_df %>%
    arrange(cluster, hcu_id, date_dispensed) %>%
    group_by(hcu_id) %>%
    mutate(
      enum = row_number(),
      tstart = if_else(enum == 1, firstevstdate + days(30), lag(date_dispensed)),
      tstop = date_dispensed,
      lagint = as.numeric(tstop - tstart)
    ) %>%
    ungroup()
  
  cc_subs_disp <- evsurv %>%
    mutate(
      censor_date = pmin(cutoff_date, firstevstdate + days(tcensor_days)),
      status = if_else(tstop < censor_date, 1, 0),
      tstop = pmin(tstop, censor_date),
      lagint = as.numeric(tstop - tstart)
    ) %>%
    select(-censor_date)
  
  if (is.null(visitnumber)) {
    cc_subs_disp <- cc_subs_disp %>%
      filter(tstop >= tstart)
  }
  
  final_intervals <- cc_subs_disp %>%
    group_by(hcu_id) %>%
    filter(enum == max(enum), lagint < tcensor_days) %>%
    mutate(
      dispence_event_id = NA_character_,
      tstart = tstop,
      tstop = pmin(cutoff_date, firstevstdate + days(tcensor_days)),
      lagint = as.numeric(tstop - tstart),
      status = 0,
      enum = enum + 1
    ) %>%
    ungroup()
  
  cc_subs_disp <- bind_rows(cc_subs_disp, final_intervals)
  
  cc_null_disp <- cohort_df %>%
    filter(!hcu_id %in% cc_subs_disp$hcu_id) %>%
    mutate(
      date_dispensed = as.Date(NA),
      dispence_event_id = NA_character_,
      enum = 1,
      status = 0,
      tstart = firstevstdate + days(30),
      tstop = pmin(cutoff_date, firstevstdate + days(tcensor_days)),
      lagint = as.numeric(tstop - tstart)
    )
  
  out <- bind_rows(cc_subs_disp, cc_null_disp) %>%
    standardise_survival_covariates(date_col = "date_dispensed")
  
  if (!is.null(visitnumber)) {
    out <- out %>% filter(enum == visitnumber)
  }
  
  out
}

# shared covariate formatting and reference-level ordering for Cox models
standardise_survival_covariates <- function(df, date_col) {
  df %>%
    group_by(hcu_id) %>%
    arrange(.data[[date_col]], .by_group = TRUE) %>%
    mutate(enum = row_number()) %>%
    ungroup() %>%
    mutate(
      dobmy_firstevstdate = floor(interval(start = dob_monthyear, end = firstevstdate) / dyears(1)),
      age_group = case_when(
        dobmy_firstevstdate < 1 ~ "under 1y",
        dobmy_firstevstdate >= 1 & dobmy_firstevstdate < 5 ~ "1-4y",
        dobmy_firstevstdate >= 5 & dobmy_firstevstdate < 15 ~ "5-14y",
        dobmy_firstevstdate >= 15 & dobmy_firstevstdate < 20 ~ "15-19y",
        dobmy_firstevstdate >= 20 ~ "20-30y"
      ),
      Sex = fct_drop(gender_code),
      Age = dobmy_firstevstdate,
      Deprivation_Quintile = fct_relevel(as.character(dep13_q), c("1", "2", "3", "4", "5")),
      Age_Group = fct_relevel(age_group, c("15-19y", "under 1y", "1-4y", "5-14y", "20-30y")),
      cc_id = factor(cc_id, levels = c("control", "case")),
      anyMP_other = fct_relevel(anyMP_other, c("non-Maori/non-Pacific", "Maori/Pacific")),
      anyMaori_other = fct_relevel(anyMaori_other, c("non-Maori", "Maori")),
      anyPacific_other = fct_relevel(anyPacific_other, c("non-Pacific", "Pacific"))
    ) %>%
    arrange(firstevstdate, cc_id, cluster, enum)
}

# fit ttf or recurrent-event models
# note that frailty models were explored, but ultimately we chose to use PWPGT for recurrent events
fit_cox_model <- function(dat,
                          ttf_or_rec = c("ttf", "rec"),
                          cov_choice,
                          covmodel_formula_list,
                          type = NULL) {
  
  ttf_or_rec <- match.arg(ttf_or_rec)
  if (!cov_choice %in% names(covmodel_formula_list)) {
    stop("Unknown cov_choice: ", cov_choice)
  }
  
  covs_str <- covmodel_formula_list[[cov_choice]]
  additional_args <- list()
  
  if (ttf_or_rec == "ttf") {
    formula <- as.formula(paste("Surv(lagint, status) ~", covs_str))
  } else {
    formula <- switch(
      type,
      "AG"    = { additional_args <- list(method = "breslow", robust = TRUE); as.formula(paste("Surv(lagint, status) ~", covs_str)) },
      "MMR"   = as.formula(paste("Surv(lagint, status) ~", covs_str, "+ cluster(hcu_id)")),
      "PWPTT" = as.formula(paste("Surv(tstop - firstevstdate, status) ~", covs_str, "+ cluster(hcu_id) + strata(enum)")),
      "PWPGT" = as.formula(paste("Surv(lagint, status) ~", covs_str, "+ cluster(hcu_id) + strata(enum)")),
      stop("Unknown recurrent event model type: ", type)
    )
  }
  
  coxph_obj <- do.call(coxph, c(list(formula = formula, data = dat), additional_args))
  
  list(
    survfit = survfit(coxph_obj),
    coxph = coxph_obj,
    assumptions = cox.zph(coxph_obj)
  )
}

# export Cox model coefficients, confidence intervals, and Schoenfeld assumption tests
export_cox_results <- function(res, model_kind, rq, followup, cov_choice, type, output_root) {
  target_dir <- file.path(output_root, rq)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  clean_covs <- str_replace_all(cov_choice, " ", "_")
  type_tag <- if (!is.na(type)) paste0(type, "_") else ""
  prefix <- paste0(model_kind, "_", rq, "_", followup, "_", clean_covs, "_", type_tag)
  
  if (is.null(res)) {
    writeLines(paste("Model failed:", prefix), file.path(target_dir, paste0(prefix, "FAILURE.txt")))
    return(invisible(FALSE))
  }
  
  sm <- summary(res$coxph)
  writeLines(capture.output(sm), file.path(target_dir, paste0(prefix, "coxph_summary.txt")))
  
  coef_tbl <- as_tibble(sm$coefficients, rownames = "variable")
  conf_try <- try(confint(res$coxph), silent = TRUE)
  
  if (!inherits(conf_try, "try-error")) {
    conf_tbl <- as_tibble(conf_try, rownames = "variable")
    coef_tbl <- left_join(coef_tbl, conf_tbl, by = "variable")
  }
  
  write_csv(coef_tbl, file.path(target_dir, paste0(prefix, "coxph.csv")))
  
  if (!is.null(res$assumptions$table)) {
    assump_tbl <- as_tibble(res$assumptions$table, rownames = "variable")
    write_csv(assump_tbl, file.path(target_dir, paste0(prefix, "assumptions.csv")))
  }
  
  invisible(TRUE)
}

# generate Kaplan-Meier cumulative incidence summary at specified days
km_summary_table <- function(dat, times = c(365.25, 1826.25)) {
  fit <- survfit(Surv(lagint, status) ~ cc_id, data = dat)
  km <- summary(fit, times = times, extend = TRUE)
  
  tibble(
    group = sub("^cc_id=", "", km$strata),
    time_days = km$time,
    event_probability = (1 - km$surv) * 100,
    lower95 = (1 - km$upper) * 100,
    upper95 = (1 - km$lower) * 100,
    n_risk = km$n.risk,
    n_event = km$n.event
  ) %>%
    mutate(timepoint = case_when(
      abs(time_days - 365.25) < 30  ~ "1 year",
      abs(time_days - 1826.25) < 30 ~ "5 years",
      TRUE ~ paste0(round(time_days), " days")
    ))
}

# ==============================================================================
# Forest Plot Extraction & Formatting Functions
# ==============================================================================

# extract adjusted case vs control HR, 95% CI, p-value, and event counts
extract_forest_summary <- function(results_list,
                                   frames_list,
                                   outcome_labels,
                                   model_kind = c("ttf", "rec"),
                                   followup = "5y",
                                   cov_suffix = "MPvnMnP",
                                   rec_type = "PWPGT") {
  model_kind <- match.arg(model_kind)
  frame_key  <- paste(model_kind, followup)
  
  imap_dfr(outcome_labels, function(display_label, rq_key) {
    clean_rq <- str_replace_all(rq_key, " ", "_")
    model_key <- if (model_kind == "ttf") {
      paste(clean_rq, model_kind, followup, cov_suffix, sep = "_")
    } else {
      paste(clean_rq, model_kind, followup, cov_suffix, rec_type, sep = "_")
    }
    
    res   <- results_list[[model_key]]
    frame <- frames_list[[rq_key]][[frame_key]]
    if (is.null(res) || is.null(frame)) return(NULL)
    
    sm <- summary(res$coxph)
    ci <- exp(confint(res$coxph))
    
    # locate the cc_idcase row
    row_idx <- grep("^cc_id", rownames(sm$coefficients))[1]
    p_col   <- intersect(c("Pr(>|z|)", "Pr(>|W|)"), colnames(sm$coefficients))[1]
    
    n_cases <- frame %>% filter(cc_id == "case", status == 1) %>% nrow()
    n_ctrls <- frame %>% filter(cc_id == "control", status == 1) %>% nrow()
    
    tibble(
      model           = display_label,
      n_events_cases  = n_cases,
      n_events_controls = n_ctrls,
      exp_coef        = sm$coefficients[row_idx, "exp(coef)"],
      lower_95        = ci[row_idx, 1],
      upper_95        = ci[row_idx, 2],
      p_val           = sm$coefficients[row_idx, p_col]
    )
  }) %>%
    mutate(
      `HR (95% CI)` = sprintf("%.2f (%.2f to %.2f)", exp_coef, lower_95, upper_95),
      pval_str = case_when(
        is.na(p_val)  ~ "",
        p_val < 0.001 ~ "<0.001",
        TRUE          ~ formatC(p_val, digits = 3, format = "f")
      ),
      stars = case_when(
        !is.na(p_val) & p_val < 0.001 ~ "***",
        !is.na(p_val) & p_val < 0.01  ~ "**",
        !is.na(p_val) & p_val < 0.05  ~ "*",
        TRUE                          ~ ""
      ),
      `P value` = str_trim(paste(pval_str, stars)),
      `    `    = paste(rep(" ", 20), collapse = "    "),
      `  `      = " "
    ) %>%
    rename(Cases = n_events_cases, Controls = n_events_controls)
}

# shared forestploter theme
manuscript_forest_theme <- function() {
  forest_theme(
    base_size   = 12,
    core        = list(bg_params = list(fill = c("grey95", "white"))),
    ci_pch      = 15,
    ci_col      = "grey10",
    ci_fill     = "black",
    ci_alpha    = 1,
    ci_lty      = 1,
    ci_lwd      = 1,
    ci_Theight  = 0.2,
    refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
    vertline_gp = gpar(lwd = 1, lty = "dashed", col = "grey20"),
    summary_fill = "grey",
    summary_col  = "grey",
    footnote_gp = gpar(fontsize = 11, fontface = "bold", col = "black")
  )
}

# build a single formatted forest plot panel
build_forest_panel <- function(summary_df,
                               domain = c("hosp", "pharm"),
                               title,
                               xlim,
                               ticks_at,
                               show_pvalue = FALSE,
                               footnote = "Total n cases = 215; Total n controls = 2150") {
  domain <- match.arg(domain)
  
  cat_header <- if (domain == "hosp") {
    "Hospitalisation event category          n events:    "
  } else {
    "Dispensing event category               n events:    "
  }
  
  plot_tbl <- summary_df %>%
    select(model, Cases, Controls, `    `, `  `, `HR (95% CI)`, any_of(if (show_pvalue) "P value")) %>%
    rename(!!cat_header := model)
  
  p <- forest(
    plot_tbl,
    est       = summary_df$exp_coef,
    lower     = summary_df$lower_95,
    upper     = summary_df$upper_95,
    sizes     = 1,
    ci_column = 4,
    ref_line  = 1,
    title     = title,
    arrow_lab = c("Lesser Risk", "Greater Risk"),
    xlim      = xlim,
    ticks_at  = ticks_at,
    footnote  = footnote,
    theme     = manuscript_forest_theme()
  )
  
  p$heights <- rep(unit(10, "mm"), nrow(p))
  p
}