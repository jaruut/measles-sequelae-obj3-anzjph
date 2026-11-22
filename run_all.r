# Master pipeline script

# open the .Rproj file before running so working directory is set to project root.

library(tidyverse)
library(lubridate)
library(readxl)
library(janitor)
library(progress)
library(MatchIt)
library(survival)
library(survminer)
library(geepack)
library(geepack)
library(forestploter)
library(patchwork)

set.seed(4402373)

# create standard directory structure if not present
dir.create("Data/Processed", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Demographics", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Figures", recursive = TRUE, showWarnings = FALSE)

# load shared functions
source("R/functions.R")

# run pipeline sequentially (or run individual scripts using saved intermediate .rds files)
source("scripts/01_load_and_cohort.R")
source("scripts/02_matching.R")
source("scripts/03_hospitalisations.R")
source("scripts/04_pharmaceuticals.R")
source("scripts/05_descriptives_and_plots.R")