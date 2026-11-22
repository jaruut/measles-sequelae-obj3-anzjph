# Long-term Hospitalisation and Pharmaceutical Outcomes Following Measles Infection in Aotearoa New Zealand (2006–2018)
Author: Vix (github: jaruut)
Permanent Archive Figshare: DOI [10.17608/k6.auckland.33315432](https://10.17608/k6.auckland.33315432)

This repository contains the R analysis scripts accompanying our manuscript submitted to the *Australian and New Zealand Journal of Public Health*.

## AI Assistance Disclosure
This code is provided primarily for the purpose of workflow transparency; the R codebase in this repository was refactored and modularised from the original working scripts using Google Gemini (Google LLC; accessed November 2026) to improve readability, remove local file paths and identifying comments, and structure the analytical pipeline for public archiving. 

## Data Availability & Privacy
The analyses use de-identified administrative health records provided by Health New Zealand | Te Whatu Ora (National Minimum Dataset [NMDS], Primary Health Organisation [PHO] Enrolment Collection, Pharmaceutical Collection [PHARMS], and National Health Index [NHI] demographics). Due to ethical and privacy restrictions, raw patient-level data cannot be shared publicly.

## Repository Structure
- `run_all.R`: Master execution script.
- `R/functions.R`: Custom functions for demographic processing, time-to-event data structuring, Cox proportional hazards modelling (`ttf` and `PWP-GT`), and result export.
- `scripts/01_load_and_cohort.R`: Imports Health NZ reference tables, defines the NZ resident cohort, and identifies hospitalised measles cases (<30 years, 2006–2018).
- `scripts/02_matching.R`: Performs exact 1:10 matching on birth month/year, sex, and NZDep2013 quintile.
- `scripts/03_hospitalisations.R`: Estimates baseline hospitalisation rates (GEE Poisson) and 5- and 10-year post-acute hospitalisation hazard ratios (Cox PH & PWP-GT) across ICD-10 diagnostic groupings.
- `scripts/04_pharmaceuticals.R`: Estimates baseline dispensing rates (GEE Poisson) and 5- and 10-year post-acute community pharmaceutical hazard ratios across therapeutic groupings.
- `scripts/05_descriptives_and_plots.R`: Generates cohort demographic tables, Kaplan–Meier cumulative incidence summaries, epidemic curves, and multi-panel hazard ratio forest plots (`forestploter` + `patchwork`).

## Software Requirements
- R (>= 4.2.0)
- Packages: `tidyverse`, `lubridate`, `readxl`, `janitor`, `progress`, `MatchIt`, `survival`, `survminer`, `geepack`, `grid`, `forestploter`, `patchwork`

## LICENSE
This code is released under the **Community Research and Academic Programming License (CRAPL)** (`CRAPL-LICENSE.txt`; see https://matt.might.net/articles/crapl/). It is provided as-is to support transparency and peer review of the accompanying manuscript, with no warranty or obligation of ongoing maintenance.