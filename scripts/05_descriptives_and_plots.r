# ==============================================================================
# 05: Cohort Descriptives, Kaplan-Meier Summaries, & Manuscript Figures
# ==============================================================================

# 1. export cohort demographic distribution tables (cases vs controls)
summarise_demo_variable <- function(df, var) {
  df %>%
    count(cc_id, {{ var }}) %>%
    group_by(cc_id) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = cc_id,
      values_from = c(n, prop),
      values_fill = 0
    )
}

write_csv(summarise_demo_variable(final_matched_data_10cc, age_group),            "Output/Demographics/cohort_age_group_proportions.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, gender_code),          "Output/Demographics/cohort_sex_proportions.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, dep13_q),              "Output/Demographics/cohort_dep13q_proportions.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, dhb_2015),             "Output/Demographics/cohort_dhb_proportions.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, dhbregion),            "Output/Demographics/cohort_dhb_region_proportions.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, priority_ethnic_code), "Output/Demographics/cohort_prioritised_ethnicity.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, anyMaori_other),       "Output/Demographics/cohort_maori_vs_nonmaori.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, anyPacific_other),     "Output/Demographics/cohort_pacific_vs_nonpacific.csv")
write_csv(summarise_demo_variable(final_matched_data_10cc, anyMP_other),          "Output/Demographics/cohort_mp_vs_nonmp.csv")

# 2. epidemic curve of index measles hospitalisations (2006–2018)
epicurve_plot <- cases_demo %>%
  mutate(evyear = year(firstevstdate)) %>%
  ggplot(aes(x = evyear)) +
  geom_histogram(binwidth = 1, boundary = 2005.5, colour = "white", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
  labs(x = "Year of index measles hospitalisation", y = "Number of cases") +
  theme_minimal()

ggsave("Output/Figures/index_measles_epicurve.png", epicurve_plot, width = 7, height = 4.5, dpi = 300)

# 3. unadjusted Kaplan-Meier cumulative incidence at 1 year and 5 years
km_summary_table(hospframes_for_surv_analysis[["all pneumonia"]][["ttf 5y"]]) %>%
  write_csv("Output/Demographics/km_1y_5y_all_pneumonia.csv")

km_summary_table(hospframes_for_surv_analysis[["respiratory infection"]][["ttf 5y"]]) %>%
  write_csv("Output/Demographics/km_1y_5y_respiratory_infection.csv")

# 4. combined 5-Year Hospitalisation Forest Plot (Panel A: ttf, Panel B: recurrent PWP-GT)
#    adjusted for sex, age group, NZDep2013 quintile, and Maori/Pacific ethnicity (MPvnMnP)
hosp_outcome_labels <- c(
  "all cause"             = "All-cause",
  "all J"                 = "Respiratory (ICD-10 Chapter J)",
  "respiratory infection" = "Respiratory infection",
  "all pneumonia"         = "All-cause pneumonia",
  "ASH"                   = "Ambulatory sensitive hospitalisations (ASH)",
  "ASH infection"         = "Infectious ASH",
  "inverse ASH"           = "Non-infectious ASH"
)

ttf_5y_hosp_df <- extract_forest_summary(
  results_list   = hosp_results_list,
  frames_list    = hospframes_for_surv_analysis,
  outcome_labels = hosp_outcome_labels,
  model_kind     = "ttf",
  followup       = "5y",
  cov_suffix     = "MPvnMnP"
)

rec_5y_hosp_df <- extract_forest_summary(
  results_list   = hosp_results_list,
  frames_list    = hospframes_for_surv_analysis,
  outcome_labels = hosp_outcome_labels,
  model_kind     = "rec",
  followup       = "5y",
  cov_suffix     = "MPvnMnP",
  rec_type       = "PWPGT"
)

p_ttf_5y_hosp <- build_forest_panel(
  summary_df = ttf_5y_hosp_df,
  domain     = "hosp",
  title      = "Time to first hospitalisation event post-measles",
  xlim       = c(0, 6),
  ticks_at   = c(0.5, 1, 2, 3, 4, 5)
)

p_rec_5y_hosp <- build_forest_panel(
  summary_df = rec_5y_hosp_df,
  domain     = "hosp",
  title      = "Recurrent hospitalisation events post-measles",
  xlim       = c(0, 6),
  ticks_at   = c(0.5, 1, 2, 3, 4, 5, 6)
)

combined_hosp_p <- wrap_elements(p_ttf_5y_hosp) /
  wrap_elements(p_rec_5y_hosp) +
  plot_annotation(tag_levels = "A")

ggsave("Output/Figures/forest_5y_hospitalisations.png", combined_hosp_p, width = 11, height = 10, dpi = 300)

# 5. combined 5-Year Pharmaceutical Dispensing Forest Plot (Panel A: ttf, Panel B: recurrent PWP-GT)
#    note: 04_pharmaceuticals.R stores frames as "ttf_5y" / "rec_5y"; align keys for extraction
phh_frames_aligned <- map(phh_dataframes, ~ set_names(.x, str_replace(names(.x), "_", " ")))

pharm_outcome_labels <- c(
  "anyabsbucket"       = "Any antibacterial",
  "oralrespinfbucket"  = "Oral respiratory/systemic antibacterial",
  "urinarybucket"      = "Urinary tract antibacterial",
  "steroidsrespbucket" = "Respiratory corticosteroids / short-course oral prednisone"
)

ttf_5y_pharms_df <- extract_forest_summary(
  results_list   = phh_results_list,
  frames_list    = phh_frames_aligned,
  outcome_labels = pharm_outcome_labels,
  model_kind     = "ttf",
  followup       = "5y",
  cov_suffix     = "MPvnMnP"
)

rec_5y_pharms_df <- extract_forest_summary(
  results_list   = phh_results_list,
  frames_list    = phh_frames_aligned,
  outcome_labels = pharm_outcome_labels,
  model_kind     = "rec",
  followup       = "5y",
  cov_suffix     = "MPvnMnP",
  rec_type       = "PWPGT"
)

p_ttf_5y_pharms <- build_forest_panel(
  summary_df = ttf_5y_pharms_df,
  domain     = "pharm",
  title      = "Time to first pharmaceutical dispensing event post-measles",
  xlim       = c(0.5, 2.5),
  ticks_at   = c(0.5, 1, 2)
)

p_rec_5y_pharms <- build_forest_panel(
  summary_df = rec_5y_pharms_df,
  domain     = "pharm",
  title      = "Recurrent pharmaceutical dispensing events post-measles",
  xlim       = c(0.5, 2.5),
  ticks_at   = c(0.5, 1, 2)
)

combined_pharms_p <- wrap_elements(p_ttf_5y_pharms) /
  wrap_elements(p_rec_5y_pharms) +
  plot_annotation(tag_levels = "A")

ggsave("Output/Figures/forest_5y_pharmaceuticals.png", combined_pharms_p, width = 11, height = 8, dpi = 300)

