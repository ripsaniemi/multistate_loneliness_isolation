
##############################
#### MULTISTATE ANALYSES #####
##############################

# Code for paper by Komulainen et al.
# Loneliness and social isolation in transitions to adverse health conditions 
# and mortality: an analysis of data from the UK Biobank study

# This file includes The needed packages, reading in the data, running the analyses and saving the results in .xlsx-files.

####################################

# PACKAGES
library(haven)
library(tidyverse)
library(mstate)
library(ggplot2)
library(purrr)
library(broom)
library(writexl)
library(fastDummies)

# set working directory
setwd("")

# read in the variables and functions of the other file
source("ms_functions.R")

# read in the data
dat <- read_dta("analysis_data.dta")

# Parts in this file
# 1. Loneliness
# 2.

#############################################################
#################### 1. LONELINESS ##########################
#############################################################

########### run the models and save them in a list ##########

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
lon_res_list <- c(ms_analyseslon_, ms_analyses_lon_C, ms_analyses_lon_female, ms_analyses_lon_male)

# change the order
lon_res_list <- lon_res_list[diags_order]

#####################################################

########### Save estimates for loneliness ###########

# mod tables of both models
lon_res_tab_m1 <- lapply(lon_res_list, mod_table, model = "c1")
lon_res_tab_m0 <- lapply(lon_res_list, mod_table, model = "c0")

# comb m0 and m1
lon_tabs_m0m1 <- list(m0 = lon_res_tab_m0, m1 = lon_res_tab_m1)

# combine diagnoses to a same table
lon_tabs_m0m1 <- lapply(lon_tabs_m0m1, res_list_to_tab)


#############################################################

################ 2. SOCIAL ISOLATION ########################

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

########### Save estimates for loneliness ###########

# modify tables of both models
si_res_tab_m1 <- lapply(si_res_list, mod_table, model = "c1")
si_res_tab_m0 <- lapply(si_res_list, mod_table, model = "c0")

# comb m0 and m1
si_tabs_m0m1 <- list(m0 = si_res_tab_m0, m1 = si_res_tab_m1)

# combine diagnoses to a same table
si_tabs_m0m1 <- lapply(si_tabs_m0m1, res_list_to_tab)

############## save files ##########

sheets <- list("M0_loneliness" = lon_tabs_m0m1[["m0"]],
               "M1_loneliness" = lon_tabs_m0m1[["m1"]],
               "M0_social_iso" = si_tabs_m0m1[["m0"]],
               "M1_social_iso" = si_tabs_m0m1[["m1"]])

write_xlsx(sheets, "multistate_cox_results.xlsx")

###################################################

# Markov assumption diagnostic models

# Time to diagnosis in transition 3 should not matter, HR should be close to 1 and/or unsignificant

# loneliness
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

####################

# isolation
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

#################################
####### format tables ######

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


sheets <- list("Loneliness" = markov_lon_hr, "Social isolation" = markov_si_hr)

write_xlsx(sheets, "table_markov_hr.xlsx")


####################################################################################

