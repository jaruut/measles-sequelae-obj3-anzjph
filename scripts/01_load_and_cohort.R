# ==============================================================================
# 01: Load Reference Tables, Define NZ Resident Cohort, Identify Measles Cases
# ==============================================================================

# 1. load Health NZ coding reference tables
Level2_and_Priority_EthnicCodes <- read_csv(
  "Data/Level2EthnicCodes.csv",
  col_types = cols(EthnicGroupPriority = col_character(), EthnicGroupCode = col_character())
) %>% clean_names()

admission_sourceCodes    <- read_xls("Data/admission-source.xls") %>% clean_names()
admission_typeCodes      <- read_xls("Data/admission-type.xls") %>% clean_names()
MDC_Codes                <- read_excel("Data/AR-DRG-v10.0_MDC_code_tab.xlsx", sheet = "MDC", skip = 2) %>% clean_names()
clin_sysCodes            <- read_xls("Data/clinical-coding-system-july-2019.xls") %>% clean_names()
domCodes                 <- read_excel("Data/domicile_code_table_census_v3_2013_0.xlsx") %>% clean_names()
event_typeCodes          <- read_xls("Data/event-type-feb-2014.xls") %>% clean_names()
event_clincode_typeCodes <- read_excel("Data/event-clinical-code-type.xls") %>% clean_names()
event_end_typeCodes      <- read_xls("Data/copy_of_event-end-type_01.xls") %>% clean_names()
facility_typeCodes       <- read_xls("Data/facility-type_0.xls") %>% clean_names()
healthspecCodes          <- read_excel("Data/health_specialty_code_table_july_2020.xls") %>% clean_names()
sex_typeCodes            <- read_xls("Data/sex-type.xls") %>% clean_names()

dhbCodes <- read_xlsx(
  "Data/dhbcodes-june-2014.xlsx",
  col_types = rep("text", 5),
  skip = 1
) %>% clean_names()

purchaserCodes <- read_excel(
  "Data/purchaser-code-july-2017.xls",
  col_types = c("text", "text", "date", "date"),
  skip = 3
) %>% clean_names()

DepLevels <- as.character(1:10)

# ICD-10-AM codes defining measles hospitalisation
measlesICD <- c(
  "B059",  # Measles without complication
  "B051",  # Measles complicated by meningitis
  "B052",  # Measles complicated by pneumonia
  "B053",  # Measles complicated by otitis media
  "B054",  # Measles with intestinal complications
  "B058",  # Measles with other complications
  "B050",  # Measles complicated by encephalitis
  "B0581", # Measles keratitis and keratoconjunctivitis
  "B0589", # Other measles complications
  "A811"   # Subacute sclerosing panencephalitis (SSPE)
)

yearly_list  <- setNames(1998:2021, as.character(1998:2021))
study_years  <- 2005:2018
active_years <- 2006:2018

# 2. load NHI cohort definition
cohort_defn <- read_delim(
  "Data/text files/mis4907_cohort_defn.txt",
  delim = "|", escape_double = FALSE,
  col_types = cols(
    DATE_OF_BIRTH = col_date(format = "%d/%m/%Y"),
    DATE_OF_DEATH = col_date(format = "%d/%m/%Y"),
    GENDER_CODE = col_factor(ordered = FALSE),
    PRIORITY_ETHNIC_CODE = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
    ETHNIC_CODE_1 = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
    ETHNIC_CODE_2 = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
    ETHNIC_CODE_3 = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
    DOMICILE_CODE = col_factor(levels = domCodes$dom),
    TLA = col_factor(levels = domCodes$tla),
    Dep01 = col_factor(levels = DepLevels, ordered = TRUE),
    Dep06 = col_factor(levels = DepLevels, ordered = TRUE),
    Dep13 = col_factor(levels = DepLevels, ordered = TRUE),
    Dep18 = col_factor(levels = DepLevels, ordered = TRUE)
  ),
  trim_ws = TRUE
) %>% clean_names()

# 3. load NMDS diagnoses and hospitalisation events (1998–2021)
all_pus_diags_txt <- sort(dir(path = "Data/text files/", pattern = "^pus11357_diags"))

all_pus_diags_df <- map(all_pus_diags_txt, function(file) {
  read_delim(
    file = file.path("Data/text files", file),
    delim = "|", escape_double = FALSE,
    col_types = cols(
      DIAG_TYP = col_factor(levels = event_clincode_typeCodes$event_clinical_code_type_code, ordered = TRUE),
      DIAG_SEQ = col_factor(),
      CLIN_SYS = col_factor(levels = clin_sysCodes$clinical_code_system),
      OP_ACDTE = col_date(format = "%d/%m/%Y"),
      OP_ACFLG = col_factor(levels = c("D", "M"))
    ),
    col_select = -c(OP_ACDTE:OP_ACFLG),
    trim_ws = TRUE
  ) %>% clean_names()
}) %>% set_names(yearly_list)

pus_diags_df <- all_pus_diags_df[as.character(study_years)]

all_pus_events_txt <- sort(dir(path = "Data/text files/", pattern = "^pus11357_events"))

all_pus_events_df <- map(all_pus_events_txt, function(file) {
  read_delim(
    file = file.path("Data/text files", file),
    delim = "|", escape_double = FALSE,
    col_types = cols(
      NZ_RES = col_factor(levels = c("Y", "N")),
      EVSTDATE = col_date(format = "%d/%m/%Y"),
      EVENDATE = col_date(format = "%d/%m/%Y"),
      short_stay_ED = col_factor(levels = c("Y", "N")),
      ADM_SRC = col_factor(levels = admission_sourceCodes$admission_source_code),
      ADM_TYPE = col_factor(levels = admission_typeCodes$admission_type_code),
      DOB = col_date(format = "%d/%m/%Y"),
      GENDER = col_factor(levels = sex_typeCodes$sex_type_code),
      ETHNICGP = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
      ETHNICG1 = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
      ETHNICG2 = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
      ETHNICG3 = col_factor(levels = Level2_and_Priority_EthnicCodes$ethnic_group_code),
      DOM_CD = col_factor(levels = domCodes$dom),
      TLA = col_factor(levels = domCodes$tla),
      Dep01 = col_factor(levels = DepLevels, ordered = TRUE),
      Dep06 = col_factor(levels = DepLevels, ordered = TRUE),
      Dep13 = col_factor(levels = DepLevels, ordered = TRUE),
      Dep18 = col_factor(levels = DepLevels, ordered = TRUE),
      EVENT_TYPE = col_factor(levels = event_typeCodes$event_type),
      END_TYPE = col_factor(levels = event_end_typeCodes$event_end_type_code),
      EVNTLVD = col_double(),
      LOS = col_double(),
      FACILITY = col_skip(),
      FACILITY_TYPE = col_factor(levels = facility_typeCodes$facility_type_code),
      HLTHSPEC = col_factor(levels = healthspecCodes$health_specialty_code),
      PURCHASER = col_factor(levels = purchaserCodes$purchaser),
      HOURS_ON_CPAP = col_double(),
      TOTAL_ICU_HOURS = col_double(),
      CCL = col_factor(levels = as.character(0:4), ordered = TRUE),
      ACC_FLAG = col_factor(levels = c("N", "Y")),
      ACC_CLAIM_NUMBER_REPORTED_FLAG = col_factor(levels = c("N", "Y")),
      PCCL = col_factor(levels = as.character(0:4), ordered = TRUE),
      DRG_GROUPER_TYPE = col_factor(levels = sprintf("%02d", 1:9)),
      MDC_CODE = col_factor(levels = MDC_Codes$mdc_code),
      MDC_TYPE = col_factor(levels = LETTERS[1:7])
    ),
    col_select = -c(AGENCY:COST_WEIGHT_CODE),
    trim_ws = TRUE
  ) %>%
    clean_names() %>%
    filter(nz_res == "Y") %>%
    mutate(
      dob_monthyear = floor_date(dob, "month"), .after = dob,
      evyear = year(evstdate), .after = evstdate
    ) %>%
    rename(hcu_id = new_master_encrypted_hcu_id)
}) %>% set_names(yearly_list)

pus_events_df <- all_pus_events_df[as.character(study_years)]

# 4. identify NZ residents via NMDS resident flag or PHO enrolment (2003 Q4 – 2019 Q1)
nzresidents_nmds <- map_dfr(all_pus_events_df, ~ select(.x, hcu_id)) %>% distinct()

all_pho_txt <- sort(dir(path = "Data/text files/", pattern = "^pho0776_2"))[1:62]

pho_hcu_id <- map_dfr(all_pho_txt, function(file) {
  read_delim(
    file = file.path("Data/text files", file),
    delim = "|", escape_double = FALSE,
    col_types = cols_only(NEW_MASTER_ENCRYPTED_HCU_ID = col_character()),
    trim_ws = TRUE
  ) %>%
    clean_names() %>%
    rename(hcu_id = new_master_encrypted_hcu_id) %>%
    distinct()
}) %>% distinct()

nzresidents_all <- bind_rows(nzresidents_nmds, pho_hcu_id) %>% distinct()

mycohort_defn <- cohort_defn %>%
  rename(hcu_id = new_master_encrypted_hcu_id) %>%
  filter(hcu_id %in% nzresidents_all$hcu_id) %>%
  mutate(dob_monthyear = floor_date(date_of_birth, unit = "month"), .after = date_of_birth) %>%
  distinct(hcu_id, .keep_all = TRUE)

rm(cohort_defn, nzresidents_nmds, pho_hcu_id, nzresidents_all, all_pho_txt)
gc()

# 5. identify all historical measles hospitalisations (to exclude from control pool)
#    and index measles hospitalisations (<30 years, 2006–2018)
all_pus_diags_measles_df <- map(all_pus_diags_df, ~ filter(.x, clin_cd %in% measlesICD))

all_pus_events_measles_df <- map2(
  all_pus_events_df,
  all_pus_diags_measles_df,
  ~ filter(.x, event_id %in% .y$event_id)
)

all_measles_events <- bind_rows(all_pus_events_measles_df) %>%
  mutate(age_at_evst = interval(dob, evstdate) %/% years(1))

# any individual with a measles hospitalisation across 1998–2021 (any age) is ineligible as a control
measlescases_any_hcu_id <- all_measles_events %>% distinct(hcu_id)

# study cases: <30 years at admission, first measles hospitalisation between 2006 and 2018
forevst_measlescohort_2006_2018 <- mycohort_defn %>%
  inner_join(
    all_measles_events %>%
      filter(age_at_evst < 30) %>%
      select(hcu_id:adm_type, event_type:age_at_evst),
    by = "hcu_id"
  ) %>%
  rename(dom_cd = domicile_code) %>%
  filter(between(evstdate, as.Date("2006-01-01"), as.Date("2018-12-31"))) %>%
  arrange(evstdate, hcu_id) %>%
  distinct(hcu_id, .keep_all = TRUE) %>%
  rename(firstevstdate = evstdate) %>%
  mutate(
    cc_id = "case", .after = hcu_id,
    dep13_q = fct_collapse(dep13, "1" = c("1", "2"), "2" = c("3", "4"), "3" = c("5", "6"), "4" = c("7", "8"), "5" = c("9", "10"))
  )

cases_demo <- forevst_measlescohort_2006_2018 %>%
  select(
    hcu_id, cc_id, dob_monthyear, gender_code, dep13, dep13_q, dhb_2015,
    date_of_birth, date_of_death, tla, nz_res, dom_cd,
    priority_ethnic_code, ethnic_code_1, ethnic_code_2, ethnic_code_3, firstevstdate
  )

rm(all_pus_diags_df, all_pus_events_df, all_pus_diags_measles_df, all_pus_events_measles_df, all_measles_events)
gc()