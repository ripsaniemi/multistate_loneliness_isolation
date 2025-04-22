
##############################
#### MULTISTATE ANALYSES #####
##############################

# Code for paper by Komulainen et al. (2025)
# Loneliness and social isolation in transitions to adverse health conditions and mortality: an analysis of data from the UK Biobank study

# Code written by Ripsa Niemi & Mai Gutvilig
# 22.04.2025

# This file includes the needed packages, reading in the data, running the analyses and saving the results in .xlsx-files.

#############################################################

# PACKAGES
library(haven)
library(tidyverse)
library(mstate)
library(ggplot2)
library(purrr)
library(broom)
library(writexl)
library(fastDummies)

# Set working directory
setwd("")

# Read in the variables and functions of the other file
source("ms_functions.R")

#############################################################

# Read in the data
dat <- read_dta("analysis_data.dta")

# About the data
# E stands for Any endocrine, nutritional or metabolic disease
# I stands for Any circulatory disease
# C stands for Any neoplasm

# Variable definitions in the data
# sex = 0/1 (female/male)
# white = 0/1/missing (other ethnicity/white ethnicity/missing)
# depgroup = 0/1/2/missing (area deprivation lest/average/most/missing)
# edugroup = 0/1/2/missing (no secondary education/secondary education/university degree/missing)
# loneliness = 0/1 (not lonely/lonely at baseline)
# isolate_bin = 0/1 (not isolated/isolated at baseline)
# time = time (days) from baseline to death in E & I analyses
# time_c = time (days) from baseline to death in C analyses
# death = 0/1 died in E & I analyses
# death_c = 0/1 died in C analyses

# diagnosis-specific variables (E as an example)
# time_E = time (days) from baseline to E diagnosis 
# E_true = 0/1/. (no diagnosis/has diagnosis/excluded due to diagnosis prior baseline)

#############################################################

# Parts in this file
# 1. Running loneliness main analyses
# 2. Running social isolation main analyses
# 3. Saving the main model results
# 4. Run Markov assumption diagnostic models for loneliness
# 5. Run Markov assumption diagnostic models for social isolation
# 6. Saving the diagnostic models 

#############################################################

# 1. Running loneliness main analyses

#############################################################
#################### LONELINESS #############################
#############################################################

# run the models and save them in a list

# E & I
ms_analyses_lon <- lapply(list_EI2, ms_analysis_lon, time_var = "time", death_var = "death")

# cancers
ms_analyses_lon_C <- lapply(list_C2, ms_analysis_lon, time_var = "time_c", death_var = "death_c")

# sex-specific cancers
ms_analyses_lon_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_lon_sex, 
                             time_var = "time_c", death_var = "death_c", sex_val = 0)
ms_analyses_lon_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_lon_sex, 
                           time_var = "time_c", death_var = "death_c", sex_val = 1)

# combine lists
lon_res_list <- c(ms_analyses_lon, ms_analyses_lon_C, ms_analyses_lon_female, ms_analyses_lon_male)

# order list items
lon_res_list <- lon_res_list[diags_order]

#####################################################

# Get estimates for loneliness

# mod tables of both models
lon_res_tab_m0 <- lapply(lon_res_list, mod_table, model = "m0")
lon_res_tab_m1 <- lapply(lon_res_list, mod_table, model = "m1")

# comb m0 and m1
lon_tabs_m0m1 <- list(m0 = lon_res_tab_m0, m1 = lon_res_tab_m1)

# combine diagnoses to a same table
lon_tabs_m0m1 <- lapply(lon_tabs_m0m1, res_list_to_tab)

#############################################################

# 2. Running social isolation main analyses

#############################################################
################ SOCIAL ISOLATION ###########################
#############################################################

########### run the models and save them in a list ##########

# E & I
ms_analyses_si <- lapply(list_EI2, ms_analysis_si, time_var = "time", death_var = "death")

# cancers
ms_analyses_si_C <- lapply(list_C2, ms_analysis_si, time_var = "time_c", death_var = "death_c")

# sex-specific cancers
ms_analyses_si_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_si_sex, 
                                time_var = "time_c", death_var = "death_c", sex_val = 0)
ms_analyses_si_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_si_sex, 
                              time_var = "time_c", death_var = "death_c", sex_val = 1)

# combine lists
si_res_list <- c(ms_analyses_si, ms_analyses_si_C, ms_analyses_si_female, ms_analyses_si_male)

# order list items
si_res_list <- si_res_list[diags_order]

#####################################################

# Get estimates for socia isolation

# modify tables of both models
si_res_tab_m0 <- lapply(si_res_list, mod_table, model = "m0")
si_res_tab_m1 <- lapply(si_res_list, mod_table, model = "m1")

# comb m0 and m1
si_tabs_m0m1 <- list(m0 = si_res_tab_m0, m1 = si_res_tab_m1)

# combine diagnoses to a same table
si_tabs_m0m1 <- lapply(si_tabs_m0m1, res_list_to_tab)

#############################################################

# 3. Saving the main model results

sheets <- list("M0_loneliness" = lon_tabs_m0m1[["m0"]],
               "M1_loneliness" = lon_tabs_m0m1[["m1"]],
               "M0_social_iso" = si_tabs_m0m1[["m0"]],
               "M1_social_iso" = si_tabs_m0m1[["m1"]])

write_xlsx(sheets, "multistate_cox_results.xlsx")

###################################################

# 4. Run Markov assumption diagnostic models for loneliness

# Markov assumption = transitions are independent of the previous transitions 
# HR for time to diagnosis in transition C) morbidity -> mortality should be close to 1 and/or insignificant

# E & I
markov_lon <- lapply(list_EI2, ms_analysis_lon_markov, time_var = "time", death_var = "death")

# cancers
markov_lon_C <- lapply(list_C2, ms_analysis_lon_markov, time_var = "time_c", death_var = "death_c")

# sex-specific cancers
markov_lon_C_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_lon_markov_sex, 
                              time_var = "time_c", death_var = "death_c", sex_val = 0)
markov_lon_C_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_lon_markov_sex, 
                            time_var = "time_c", death_var = "death_c", sex_val = 1)

# combine all
markov_lon_all <- c(markov_lon, markov_lon_C, markov_lon_C_female, markov_lon_C_male)

###################################################

# 5. Run Markov assumption diagnostic models for social isolation

# E & I
markov_si <- lapply(list_EI2, ms_analysis_si_markov, time_var = "time", death_var = "death")

# cancers
markov_si_C <- lapply(list_C2, ms_analysis_si_markov, time_var = "time_c", death_var = "death_c")

# sex-specific cancers
markov_si_C_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_si_markov_sex, 
                             time_var = "time_c", death_var = "death_c", sex_val = 0)
markov_si_C_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_si_markov_sex, 
                           time_var = "time_c", death_var = "death_c", sex_val = 1)

# combine all
markov_si_all <- c(markov_si, markov_si_C, markov_si_C_female, markov_si_C_male)

###################################################

# 6. Saving the diagnostic models 

# format tables

# loneliness
markov_lon_hr <- lapply(markov_lon_all, get_time_diag_hr)

markov_lon_hr <- imap(markov_lon_hr, ~.x %>%
                        mutate(diagnosis = .y)) |>
  bind_rows()

# social isolation
markov_si_hr <- lapply(markov_si_all, get_time_diag_hr)

markov_si_hr <- imap(markov_si_hr, ~.x %>%
                       mutate(diagnosis = .y)) |>
  bind_rows()


sheets <- list("Loneliness" = markov_lon_hr, 
               "Social isolation" = markov_si_hr)

write_xlsx(sheets, "table_markov_hr.xlsx")

###################################################
################### CODE ENDS #####################
###################################################

