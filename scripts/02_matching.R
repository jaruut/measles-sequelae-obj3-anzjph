# ==============================================================================
# 02: Exact 1:10 Case-Control Matching & Demographic Standardisation
# ==============================================================================

# build eligible control pool (NZ residents never hospitalised for measles, non-missing NZDep13)
control_cohort_demo <- mycohort_defn %>%
  filter(!hcu_id %in% measlescases_any_hcu_id$hcu_id, !is.na(dep13)) %>%
  distinct(hcu_id, .keep_all = TRUE) %>%
  rename(dom_cd = domicile_code) %>%
  mutate(
    cc_id = "control",
    nz_res = factor("Y", levels = c("Y", "N")),
    dep13_q = fct_collapse(dep13, "1" = c("1", "2"), "2" = c("3", "4"), "3" = c("5", "6"), "4" = c("7", "8"), "5" = c("9", "10"))
  ) %>%
  select(
    hcu_id, cc_id, dob_monthyear, gender_code, dep13, dep13_q,
    date_of_birth, date_of_death, dhb_2015, tla, nz_res, dom_cd,
    priority_ethnic_code, ethnic_code_1, ethnic_code_2, ethnic_code_3
  )

allcohort_full <- bind_rows(cases_demo, control_cohort_demo) %>%
  arrange(dob_monthyear, gender_code, dep13_q)

# exact matching on birth month/year, sex, and NZDep2013 quintile
matchit_setup <- matchit(
  cc_id ~ dob_monthyear + gender_code + dep13_q,
  data = allcohort_full,
  method = "exact",
  replace = FALSE
)

matched_data <- match.data(matchit_setup)

# sample 10 controls without replacement per case within exact match strata
control_pool <- matched_data %>% filter(cc_id == "control")
case_ids <- unique(matched_data$hcu_id[matched_data$cc_id == "case"])

final_matched_data <- map_dfr(case_ids, function(case_patient) {
  matched_case <- matched_data %>% filter(hcu_id == case_patient)
  
  matched_controls <- control_pool %>%
    filter(
      dob_monthyear == matched_case$dob_monthyear,
      gender_code == matched_case$gender_code,
      dep13_q == matched_case$dep13_q
    )
  
  if (nrow(matched_controls) > 10) {
    matched_controls <- slice_sample(matched_controls, n = 10)
    control_pool <<- control_pool %>% filter(!hcu_id %in% matched_controls$hcu_id)
  }
  
  bind_rows(matched_case, matched_controls)
})

# assign cluster IDs and propagate the index case's admission date (firstevstdate) to matched controls
final_matched_data_10cc <- final_matched_data %>%
  mutate(cluster = as.double(str_extract(subclass, "^\\d+"))) %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  mutate(
    firstevstdate = first(firstevstdate[str_extract(subclass, "_\\d+_") == "_1_"]),
    .after = ethnic_code_3
  ) %>%
  ungroup() %>%
  process_demo()

# verify 1:10 matching ratio across all clusters
cluster_counts <- final_matched_data_10cc %>%
  count(cluster, cc_id) %>%
  pivot_wider(names_from = cc_id, values_from = n, values_fill = 0)

if (any(cluster_counts$control < 10)) {
  warning("At least one case has fewer than 10 matched controls.")
} else {
  message("Verified: Every case has 10 matched controls.")
}

cases_demo <- final_matched_data_10cc %>% filter(cc_id == "case")
ctrls_demo <- final_matched_data_10cc %>% filter(cc_id == "control")

write_csv(final_matched_data_10cc, "Data/Processed/final_matched_data_10cc.csv")

rm(matched_data, matchit_setup, final_matched_data, control_cohort_demo, control_pool, allcohort_full)
gc()