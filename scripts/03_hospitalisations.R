# ==============================================================================
# 03: Baseline & Post-Acute Hospitalisation Analyses
# ==============================================================================

# 1. extract all NMDS hospitalisation events (2005–2018) for matched cases & controls
extract_cohort_events <- function(cohort_subset) {
  yearly_list <- imap(yearly_list[as.character(active_years)], function(index_yr, yr_name) {
    index_ids <- cohort_subset %>%
      filter(year(firstevstdate) == index_yr) %>%
      select(hcu_id)
    
    map(study_years, function(obs_yr) {
      pus_events_df[[as.character(obs_yr)]] %>%
        filter(hcu_id %in% index_ids$hcu_id) %>%
        select(-c(nz_res, dob_monthyear)) %>%
        left_join(cohort_subset, by = "hcu_id", suffix = c(".pus_events", ".demo")) %>%
        rename(
          dhb_2015_code.pus_events = dhb_2015.pus_events,
          dhb_2015_name.demo = dhb_2015.demo
        ) %>%
        mutate(base_post = categorise_followup_period(evstdate, firstevstdate), .after = evyear)
    }) %>% set_names(study_years)
  })
  
  all_events <- yearly_list %>% list_flatten() %>% bind_rows()
  
  baseline_counts <- all_events %>%
    filter(base_post == "base") %>%
    count(hcu_id, name = "baseline_visits") %>%
    right_join(cohort_subset, by = "hcu_id") %>%
    mutate(baseline_visits = replace_na(baseline_visits, 0L))
  
  list(
    yearly_events_list = yearly_list,
    all_events_frame = all_events,
    baseline_events_count = baseline_counts
  )
}

cases_hosp <- extract_cohort_events(cases_demo)
ctrls_hosp <- extract_cohort_events(ctrls_demo)

# 2. baseline hospitalisation rate comparison (1 year pre-index; GEE Poisson)
baseline_visits_df <- bind_rows(
  cases_hosp$baseline_events_count,
  ctrls_hosp$baseline_events_count
) %>%
  mutate(gender_code = fct_drop(gender_code))

baseline_hosp_use_model <- geeglm(
  baseline_visits ~ cc_id + gender_code + age_group + dep13_q,
  id = cluster,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  data = baseline_visits_df
)

b_hosp  <- coef(summary(baseline_hosp_use_model))["cc_id", "Estimate"]
se_hosp <- coef(summary(baseline_hosp_use_model))["cc_id", "Std.err"]

writeLines(
  paste(
    c("IRR_hosp", "CI95_lower", "CI95_upper", "p_hosp"),
    c(exp(b_hosp), exp(b_hosp - 1.96 * se_hosp), exp(b_hosp + 1.96 * se_hosp),
      coef(summary(baseline_hosp_use_model))["cc_id", "Pr(>|W|)"]),
    sep = "="
  ),
  "Output/Demographics/baseline-hosp-use-stats.txt"
)

baseline_visits_df %>%
  group_by(cc_id) %>%
  summarise(
    mean = mean(baseline_visits), sd = sd(baseline_visits),
    median = median(baseline_visits), iqr = IQR(baseline_visits),
    .groups = "drop"
  ) %>%
  write_csv("Output/Demographics/baseline-hosp-use-parameters.csv")

# 3. join ICD-10 diagnosis codes (clin_sys == 10) to hospitalisation events
attach_diagnosis_codes <- function(datlist) {
  map(datlist, function(year_events) {
    imap(year_events, function(events, yr_char) {
      left_join(
        events,
        filter(pus_diags_df[[yr_char]], clin_sys == 10),
        by = "event_id"
      ) %>%
        select(hcu_id, dob, event_id, base_post, diag_typ, clin_cd, clin_sys, evstdate, evendate) %>%
        distinct(hcu_id, event_id, clin_cd, .keep_all = TRUE) %>%
        left_join(final_matched_data_10cc, by = "hcu_id") %>%
        mutate(evint = time_length(interval(firstevstdate + days(30), evstdate, tz = "Pacific/Auckland"), "day"))
    })
  })
}

cases_yearly_events_diagsall <- attach_diagnosis_codes(cases_hosp$yearly_events_list)
ctrls_yearly_events_diagsall <- attach_diagnosis_codes(ctrls_hosp$yearly_events_list)

# 4. ICD-10-AM outcome definitions
cds_asthma       <- c("^J4[5-6]")
cds_cellulitis   <- c("^H000", "^H010", "^J340", "^L0[1-4,8]", "^L980")
cds_constipation <- c("^K590")
cds_dental       <- c("^K0[2,4-5]")
cds_derm_eczm    <- c("^L2", "^L30")
cds_derm         <- c("^L303")
cds_eczm         <- c("^L2")
cds_gast_dehy    <- c("^A0[2-9]", "^R11", "^K529")
cds_gord         <- c("^K21")
cds_nutdef_anaem <- c("^D5[0-3]", "^E4[0-6]", "^E5[0-9]", "^E6[0-4]")
cds_pneumonia    <- c("^J1[3-8]")
cds_urti_enti    <- c("^J0[0-4,6]", "^H6[5-7]")
cds_vp_main      <- c("^A3[3-7]", "^A403", "^A80", "^B1[6,8]")
cds_vp_mmr       <- c("^B0[5-6]", "^B26", "^M014")

cds_meningococcal <- c("^A39")
cds_pneumococcalreport <- c("^G00[1-2,9]", "^A40[3,9]", "^A49[1,9]", "^B953", "^J13", "^A419", "^M001")
cds_allcausepneumo <- c("^J1[2-8]", "^J100", "^J110")
cds_pneumobac    <- c("^J1[3-8]")
cds_spe          <- c("^A811")
cds_j            <- c("^J")
cds_allcause     <- c(".*")

RI_codes     <- c(cds_urti_enti, cds_pneumonia)
ASH_codes    <- c(cds_asthma, cds_cellulitis, cds_constipation, cds_dental, cds_derm_eczm,
                  cds_gast_dehy, cds_gord, cds_nutdef_anaem, cds_pneumonia, cds_urti_enti,
                  cds_vp_main, cds_vp_mmr)
AI_codes     <- c(cds_cellulitis, cds_derm, cds_gast_dehy, cds_pneumonia, cds_urti_enti,
                  cds_vp_main, cds_vp_mmr)
INVASH_codes <- c(cds_constipation, cds_dental, cds_eczm, cds_gord, cds_nutdef_anaem)

codes_list <- list(
  "all cause"                          = cds_allcause,
  "subacute sclerosing panencephalitis" = cds_spe,
  "all J"                              = cds_j,
  "respiratory infection"              = RI_codes,
  "ASH infection"                      = AI_codes,
  "ASH"                                = ASH_codes,
  "inverse ASH"                        = INVASH_codes,
  "meningococcal"                      = cds_meningococcal,
  "pneumococcal report codes"          = cds_pneumococcalreport,
  "all pneumonia"                      = cds_allcausepneumo,
  "bacteria pneumonia"                 = cds_pneumobac
)

# filter events by ICD-10 regex and collapse to distinct hospitalisation episodes per person
extract_coded_visits <- function(events_list, regex_vec) {
  pattern <- paste(regex_vec, collapse = "|")
  
  events_list %>%
    list_flatten() %>%
    bind_rows() %>%
    filter(str_detect(clin_cd, pattern)) %>%
    select(-c(diag_typ, clin_cd, clin_sys)) %>%
    distinct(hcu_id, event_id, .keep_all = TRUE) %>%
    group_by(hcu_id) %>%
    arrange(evstdate, .by_group = TRUE) %>%
    mutate(enum = row_number()) %>%
    ungroup() %>%
    arrange(firstevstdate, hcu_id, enum) %>%
    mutate(
      tstart = if_else(enum == 1, firstevstdate, lag(evstdate)),
      tstop = evstdate,
      lagint = as.numeric(tstop - tstart)
    )
}

cc_coded_visits_list <- map(codes_list, function(regex_vec) {
  bind_rows(
    extract_coded_visits(cases_yearly_events_diagsall, regex_vec),
    extract_coded_visits(ctrls_yearly_events_diagsall, regex_vec)
  )
})

# 5. construct survival dataframes (5-year = 1827d; 10-year = 3653d; ttf & recurrent)
param_combinations <- tibble(
  param_name   = c("ttf 5y", "ttf 10y", "rec 5y", "rec 10y"),
  visitnumber  = list(1, 1, NULL, NULL),
  tcensor_days = c(1827, 3653, 1827, 3653)
)

hospframes_for_surv_analysis <- map(cc_coded_visits_list, function(visits_df) {
  pmap(param_combinations, function(param_name, visitnumber, tcensor_days) {
    tryCatch(
      build_hosp_survival_data(visits_df, final_matched_data_10cc, visitnumber, tcensor_days),
      error = function(e) { message("Skipping frame (", param_name, "): ", e$message); NULL }
    )
  }) %>% set_names(param_combinations$param_name)
})

# 6. define covariate specifications and fit models
covmodel_formula_list <- list(
  "vanilla"               = "cc_id + Sex + Age + Deprivation_Quintile",
  "Age Group"             = "cc_id + Sex + Age_Group + Deprivation_Quintile",
  "prioritised ethnicity" = "cc_id + Sex + Age_Group + Deprivation_Quintile + priority_ethnic_code",
  "MvnM"                  = "cc_id + Sex + Age_Group + Deprivation_Quintile + anyMaori_other",
  "PvnP"                  = "cc_id + Sex + Age_Group + Deprivation_Quintile + anyPacific_other",
  "MPvnMnP"               = "cc_id + Sex + Age_Group + Deprivation_Quintile + anyMP_other",
  "unadjusted"            = "cc_id"
)

model_specs <- tibble(
  model_kind  = c("ttf", "ttf", "rec", "rec"),
  followup    = c("5y", "10y", "5y", "10y"),
  type        = c(NA_character_, NA_character_, "PWPGT", "PWPGT"),
  frame_index = c(1L, 2L, 3L, 4L)
)

hosp_plan <- expand_grid(
  rq = names(codes_list),
  cov_choice = names(covmodel_formula_list)
) %>%
  crossing(model_specs) %>%
  mutate(
    model_name = pmap_chr(list(rq, model_kind, followup, cov_choice, type), function(rq, mk, fu, cv, ty) {
      parts <- c(rq, mk, fu, cv, if (!is.na(ty)) ty)
      str_replace_all(paste(parts, collapse = "_"), " ", "_")
    })
  )

pb <- progress_bar$new(total = nrow(hosp_plan), format = "  Fitting hosp Cox models [:bar] :percent in :elapsed")

hosp_results_list <- pmap(hosp_plan, function(rq, cov_choice, model_kind, followup, type, frame_index, ...) {
  on.exit(pb$tick())
  frame <- hospframes_for_surv_analysis[[rq]][[frame_index]]
  if (is.null(frame) || nrow(frame) == 0) return(NULL)
  
  tryCatch(
    fit_cox_model(frame, model_kind, cov_choice, covmodel_formula_list, type),
    error = function(e) NULL
  )
}) %>% set_names(hosp_plan$model_name)

# 7. export hospitalisation model results
hosp_out_root <- file.path("Output", format(Sys.time(), "%Y%m%d_%H%M"), "CoxResults_hospitalisations")

pwalk(
  list(
    res = hosp_results_list,
    model_kind = hosp_plan$model_kind,
    rq = hosp_plan$rq,
    followup = hosp_plan$followup,
    cov_choice = hosp_plan$cov_choice,
    type = hosp_plan$type
  ),
  export_cox_results,
  output_root = hosp_out_root
)