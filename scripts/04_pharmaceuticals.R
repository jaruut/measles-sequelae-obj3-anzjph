# ==============================================================================
# 04: Baseline & Post-Acute Community Pharmaceutical Dispensing Analyses
# ==============================================================================

# 1. load pharmaceutical formulation reference table
phh_dfps <- read_delim(
  "Data/text files/phh1119_dimformpacksub.txt",
  delim = "|", escape_double = FALSE, trim_ws = TRUE
) %>% clean_names()

# 2. load monthly PHARMS dispensing files (2005–2018) filtered to study cohort
month_char <- sprintf("%02d", 1:12)
phh_yearmonth <- expand_grid(yr = as.character(2005:2018), mo = month_char) %>%
  unite("ym", yr, mo, sep = "") %>%
  pull(ym)

all_phh_list <- map(phh_yearmonth, function(ym) {
  filename <- list.files(path = "Data/text files", pattern = paste0("^phh_", ym), full.names = TRUE)[1]
  
  read_delim(
    filename,
    delim = "|", escape_double = FALSE,
    col_types = cols(
      DATE_DISPENSED = col_date(format = "%d/%m/%Y"),
      DISPENCE_EVENT_ID = col_character(),
      .default = col_guess()
    ),
    trim_ws = TRUE
  ) %>%
    clean_names() %>%
    rename(hcu_id = new_master_encrypted_hcu_id) %>%
    select(hcu_id, date_dispensed, dispence_event_id, dim_form_pack_subsidy_key, dispensings_prescribed, days_supply) %>%
    filter(hcu_id %in% final_matched_data_10cc$hcu_id) %>%
    inner_join(final_matched_data_10cc, by = "hcu_id") %>%
    filter(interval(date_of_birth, date_dispensed) < years(30)) %>%
    mutate(base_post = categorise_followup_period(date_dispensed, firstevstdate))
})

phh_base_all <- map_dfr(all_phh_list, ~ filter(.x, base_post == "base"))
phh_post_all <- map_dfr(all_phh_list, ~ filter(.x, base_post == "post"))

# 3. baseline pharmaceutical dispensing comparison (GEE Poisson)
baseline_counts_pharms <- phh_base_all %>%
  count(hcu_id, name = "baseline_disp") %>%
  right_join(final_matched_data_10cc, by = "hcu_id") %>%
  mutate(
    baseline_disp = replace_na(baseline_disp, 0L),
    gender_code = fct_drop(gender_code)
  )

baseline_pharm_use_model <- geeglm(
  baseline_disp ~ cc_id + gender_code + age_group + dep13_q,
  id = cluster,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  data = baseline_counts_pharms
)

b_pharm  <- coef(summary(baseline_pharm_use_model))["cc_id", "Estimate"]
se_pharm <- coef(summary(baseline_pharm_use_model))["cc_id", "Std.err"]

writeLines(
  paste(
    c("IRR_pharms", "CI95_lower", "CI95_upper", "p_pharms"),
    c(exp(b_pharm), exp(b_pharm - 1.96 * se_pharm), exp(b_pharm + 1.96 * se_pharm),
      coef(summary(baseline_pharm_use_model))["cc_id", "Pr(>|W|)"]),
    sep = "="
  ),
  "Output/Demographics/baseline-pharm-use-stats.txt"
)

baseline_counts_pharms %>%
  group_by(cc_id) %>%
  summarise(
    mean = mean(baseline_disp), sd = sd(baseline_disp),
    median = median(baseline_disp), iqr = IQR(baseline_disp),
    .groups = "drop"
  ) %>%
  write_csv("Output/Demographics/baseline-pharm-use-parameters.csv")

# 4. therapeutic group definitions (PHARMAC Therapeutic Group IDs & Chemical IDs)
dfpsk_finder <- function(level1 = NULL, level2 = NULL, level3 = NULL, chemid = NULL) {
  res <- phh_dfps
  if (!is.null(level1)) res <- res %>% filter(tg_level1_id %in% level1)
  if (!is.null(level2)) res <- res %>% filter(tg_level2_id %in% level2)
  if (!is.null(level3)) res <- res %>% filter(tg_level3_id %in% level3)
  if (!is.null(chemid)) res <- res %>% filter(chemical_id %in% chemid)
  res %>% distinct(dim_form_pack_subsidy_key) %>% pull(dim_form_pack_subsidy_key)
}

keys_topical_abs      <- dfpsk_finder(10, 4)
keys_steroids_preds   <- dfpsk_finder(14, 7, NULL, c(2034, 2036, 2038))
keys_steroids_nasal   <- dfpsk_finder(28, 51, NULL, c(1108, 1168, 3903))
keys_steroids_inhal   <- dfpsk_finder(28, c(10, 7), NULL, c(1065, 1108, 1168))
keys_uti              <- c(dfpsk_finder(16, 23), dfpsk_finder(16, 3, 13, 2300))
keys_abs_pneumo       <- dfpsk_finder(16, 3, 10, 2529)
keys_abs_csu          <- dfpsk_finder(16, 3, 13, 1361)
keys_abs_macrol       <- dfpsk_finder(16, 3, 4)
keys_abs_cephal       <- dfpsk_finder(16, 3, 1, c(1228, 1234))
keys_abs_pen_benzbenz <- dfpsk_finder(16, 3, 7, c(1512, 1070, 1072, 1511, 1970, 1120))

phh_topical_abs_post      <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_topical_abs)
phh_steroids_preds_post   <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_steroids_preds, dispensings_prescribed < 2, days_supply < 7)
phh_steroids_nasal_post   <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_steroids_nasal)
phh_steroids_inhal_post   <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_steroids_inhal)
phh_uti_post              <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_uti)
phh_abs_pneumo_post       <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_abs_pneumo, dispensings_prescribed < 2, days_supply < 7)
phh_abs_csu_post          <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_abs_csu)
phh_abs_macrol_post       <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_abs_macrol)
phh_abs_cephal_post       <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_abs_cephal)
phh_abs_pen_benzbenz_post <- phh_post_all %>% filter(dim_form_pack_subsidy_key %in% keys_abs_pen_benzbenz)

phh_df_combos_list <- list(
  "oralrespinfbucket"  = bind_rows(phh_abs_pneumo_post, phh_uti_post, phh_abs_csu_post, phh_abs_macrol_post, phh_abs_cephal_post, phh_abs_pen_benzbenz_post),
  "urinarybucket"      = bind_rows(phh_uti_post, phh_abs_csu_post),
  "steroidsrespbucket" = bind_rows(phh_steroids_nasal_post, phh_steroids_inhal_post, phh_steroids_preds_post),
  "anyabsbucket"       = bind_rows(phh_topical_abs_post, phh_abs_pneumo_post, phh_abs_csu_post, phh_abs_macrol_post, phh_abs_cephal_post, phh_abs_pen_benzbenz_post, phh_uti_post)
)

# 5. build survival frames & fit models
phh_dataframes <- map(phh_df_combos_list, function(df) {
  list(
    ttf_5y  = build_pharm_survival_data(df, final_matched_data_10cc, visitnumber = 1, tcensor_days = 1827),
    ttf_10y = build_pharm_survival_data(df, final_matched_data_10cc, visitnumber = 1, tcensor_days = 3653),
    rec_5y  = build_pharm_survival_data(df, final_matched_data_10cc, visitnumber = NULL, tcensor_days = 1827),
    rec_10y = build_pharm_survival_data(df, final_matched_data_10cc, visitnumber = NULL, tcensor_days = 3653)
  )
})

phh_plan <- expand_grid(
  rq = names(phh_df_combos_list),
  cov_choice = names(covmodel_formula_list)
) %>%
  crossing(model_specs) %>%
  mutate(
    model_name = pmap_chr(list(rq, model_kind, followup, cov_choice, type), function(rq, mk, fu, cv, ty) {
      parts <- c(rq, mk, fu, cv, if (!is.na(ty)) ty)
      str_replace_all(paste(parts, collapse = "_"), " ", "_")
    })
  )

pb <- progress_bar$new(total = nrow(phh_plan), format = "  Fitting pharm Cox models [:bar] :percent in :elapsed")

phh_results_list <- pmap(phh_plan, function(rq, cov_choice, model_kind, followup, type, frame_index, ...) {
  on.exit(pb$tick())
  frame <- phh_dataframes[[rq]][[frame_index]]
  tryCatch(
    fit_cox_model(frame, model_kind, cov_choice, covmodel_formula_list, type),
    error = function(e) NULL
  )
}) %>% set_names(phh_plan$model_name)

# 6. export pharmaceutical model results
phh_out_root <- file.path("Output", format(Sys.time(), "%Y%m%d_%H%M"), "CoxResults_pharmaceuticals")

pwalk(
  list(
    res = phh_results_list,
    model_kind = phh_plan$model_kind,
    rq = phh_plan$rq,
    followup = phh_plan$followup,
    cov_choice = phh_plan$cov_choice,
    type = phh_plan$type
  ),
  export_cox_results,
  output_root = phh_out_root
)