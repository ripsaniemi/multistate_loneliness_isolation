##################### MULTISTATE ANALYYSI BIOPANKKI-DATALLA ##########################

# aloitettu 12.10.2023 yhdistämään vanhoja tiedostoja yhteen R-fileen
# päivitys tammikuu 2025

# Ripsa ja Mai

# malli 1, jossa vakioitu ikä, sukupuoli, etnisyys, koulutus

library(haven)
library(tidyverse)
library(mstate)
library(ggplot2)
library(purrr)
library(broom)
library(writexl)
library(fastDummies)

# koodi pohjautuu tutoriaaliin: Putter (2021) Tutorial in biostatistics: Competing risks and multi-state models Analyses using the mstate package
# linkki: https://cran.r-project.org/web/packages/mstate/vignettes/Tutorial.pdf

# josta myös artikkeli Putter (2006) Tutorial in biostatistics: Competing risks and multi-state models
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/sim.2712

# kirjasto
#setwd("Z:/Documents/BioBank/Data/multistate/")
setwd("C:/Users/ripsanie/Documents/Biobank/Multistate")

# funktiot
source("ms_functions.R")

# datat, huom. kumpaa ajat (tätä ei automatisoitu)
dat <- read_dta("Z:/Documents/BioBank/Data/multistate/analysis_data.dta")
#dat <- read_dta("Z:/Documents/BioBank/Data/multistate/analysis_data_sensitivity.dta")


#################### 1. YKSINÄISYYS ##########################

#################### mallien tallennus listaan ################
# E & I
ms_analyses <- lapply(list_EI2, ms_analysis_lon, time_var = "time", death_var = "death")

# syövät
ms_analyses_C <- lapply(list_C2, ms_analysis_lon, time_var = "time_c", death_var = "death_c")

# gender syövät
ms_analyses_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_lon_sex, 
                             time_var = "time_c", death_var = "death_c", sex_val = 0)
ms_analyses_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_lon_sex, 
                           time_var = "time_c", death_var = "death_c", sex_val = 1)

lon_res_list <- c(ms_analyses, ms_analyses_C, ms_analyses_female, ms_analyses_male)

lon_res_list <- lon_res_list[diags_order]

####################

########### Mallit 0 ja 1 talteen ###########

# mod tables of both models
lon_res_tab_m1 <- lapply(lon_res_list, mod_table, model = "c1")
lon_res_tab_m0 <- lapply(lon_res_list, mod_table, model = "c0")

# comb m0 and m1
lon_tabs_m0m1 <- list(m0 = lon_res_tab_m0, m1 = lon_res_tab_m1)

# combine diagnoses to same table
lon_tabs_m0m1 <- lapply(lon_tabs_m0m1, res_list_to_tab)


########################################################################

#################### 2. SOSIAALINEN ISOLAATIO ##########################

#################### mallien tallennus listaan ################

# E & I
ms_analyses_si <- lapply(list_EI2, ms_analysis_si, time_var = "time", death_var = "death")

# syövät
ms_analyses_si_C <- lapply(list_C2, ms_analysis_si, time_var = "time_c", death_var = "death_c")

# rinta- ja eturauhassyövät
# gender syövät
ms_analyses_si_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_si_sex, 
                                time_var = "time_c", death_var = "death_c", sex_val = 0)
ms_analyses_si_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_si_sex, 
                              time_var = "time_c", death_var = "death_c", sex_val = 1)


si_res_list <- c(ms_analyses_si, ms_analyses_si_C, ms_analyses_si_female, ms_analyses_si_male)

si_res_list <- si_res_list[diags_order]

########### Mallit 0 ja 1 talteen ###########

si_res_tab_m1 <- lapply(si_res_list, mod_table, model = "c1")
si_res_tab_m0 <- lapply(si_res_list, mod_table, model = "c0")

# comb m0 and m1
si_tabs_m0m1 <- list(m0 = si_res_tab_m0, m1 = si_res_tab_m1)

si_tabs_m0m1 <- lapply(si_tabs_m0m1, res_list_to_tab)

####################### Järjestetään diagnoosit ja tallennetaan datat

m0_results_comb_lon <- lon_tabs_m0m1[["m0"]] |>
  arrange(factor(diagnosis, levels = diags_order), trans)
m1_results_comb_lon <- lon_tabs_m0m1[["m1"]] |>
  arrange(factor(diagnosis, levels = diags_order), trans)

m0_results_comb_si <- si_tabs_m0m1[["m0"]] |>
  arrange(factor(diagnosis, levels = diags_order), trans)
m1_results_comb_si <- si_tabs_m0m1[["m1"]] |>
  arrange(factor(diagnosis, levels = diags_order), trans)

sheets <- list("M0_loneliness" = m0_results_comb_lon,
               "M1_loneliness" = m1_results_comb_lon,
               "M0_social_iso" = m0_results_comb_si,
               "M1_social_iso" = m1_results_comb_si)

write_xlsx(sheets, "results/Multistate_cox_results_m0m1.xlsx")
#write_xlsx(sheets, "results/Multistate_cox_results_m0m1_sens.xlsx")

# huom, kumpaa ajamassa
#write_xlsx(sheets, "~/BioBank/Results/Multistate models/Multistate_cox_results.xlsx")
#write_xlsx(sheets, "~/BioBank/Results/Multistate models/Multistate_cox_results_sensitivity.xlsx")

###################################################

# Markov assumption

# Time to diagnosis in transition 3 should not matter, HR should be close to 1 and/or unsignificant

# ajetaan funktiot

# loneliness
# E & I
markov_lon <- lapply(list_EI2, ms_analysis_lon_markov, time_var = "time", death_var = "death")

# syövät
markov_lon_C <- lapply(list_C2, ms_analysis_lon_markov, time_var = "time_c", death_var = "death_c")

# naisten syövät
markov_lon_C_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_lon_markov_sex, 
                              time_var = "time_c", death_var = "death_c", sex_val = 0)

# miesten syövät
markov_lon_C_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_lon_markov_sex, 
                            time_var = "time_c", death_var = "death_c", sex_val = 1)

markov_lon_all <- c(markov_lon, markov_lon_C, markov_lon_C_female, markov_lon_C_male)

####################

# isolation
# E & I
markov_si <- lapply(list_EI2, ms_analysis_si_markov, time_var = "time", death_var = "death")

# syövät
markov_si_C <- lapply(list_C2, ms_analysis_si_markov, time_var = "time_c", death_var = "death_c")

# naisten syövät
markov_si_C_female <- lapply(list_C2_sex[c("C50", "C5158")], ms_analysis_si_markov_sex, 
                             time_var = "time_c", death_var = "death_c", sex_val = 0)

# miesten syövät
markov_si_C_male <- lapply(list_C2_sex[c("C6063")], ms_analysis_si_markov_sex, 
                           time_var = "time_c", death_var = "death_c", sex_val = 1)

markov_si_all <- c(markov_si, markov_si_C, markov_si_C_female, markov_si_C_male)

# funktio taulukon muotoiluun
get_time_diag_hr <- function(model) {
  
  hr <- model[["c2"]] |> 
    tidy(conf.int = TRUE, exponentiate = TRUE) |> 
    filter(term == "time_diag.3") |>
    mutate(HR_CI = paste0(round(estimate, 3), " (", round(conf.low, 3), "–", round(conf.high, 3), ")")) |>
    select(HR_CI)
}

# loneliness
markov_lon_hr <- lapply(markov_lon_all, get_time_diag_hr)

markov_lon_hr <- imap(markov_lon_hr, ~.x %>%
                        mutate(diagnosis = .y)) |>
  bind_rows() |>
  mutate(Diagnosis = case_when(
    diagnosis == "C" ~ "Any neoplasm",
    diagnosis == "C0014" ~ "Malignant neoplasms of lip, oral cavity and pharynx (C00–C14)",
    diagnosis == "C1526" ~ "Malignant neoplasms of digestive organs (C15–C26)",
    diagnosis == "C3039" ~ "Malignant neoplasms of respiratory and intrathoracic organs (C30–C39)",
    diagnosis == "C4041" ~ "Malignant neoplasms of bone and articular cartilage (C40–C41)",
    diagnosis == "C4344" ~ "Melanoma & other malignant neoplasms of skin (C43–C44)",
    diagnosis == "C4549" ~ "Malignant neoplasms of mesothelial and soft tissue (C45–C49)",
    diagnosis == "C50" ~ "Malignant neoplasms of breast (C50)",
    diagnosis == "C5158" ~ "Malignant neoplasms of female genital organs (C51–C58)",
    diagnosis == "C6063" ~ "Malignant neoplasms of male genital organs (C60–C63)",
    diagnosis == "C6468" ~ "Malignant neoplasms of urinary tract (C64–C68)",
    diagnosis == "C6972" ~ "Malignant neoplasms of eye, brain and other parts of central nervous system (C69–C72)",
    diagnosis == "C7375" ~ "Malignant neoplasms of thyroid and other endocrine glands (C73–C75)",
    diagnosis == "C8196" ~ "Malignant neoplasms of lymphoid, hematopoietic and related tissue (C81–C96)",
    diagnosis == "D0009" ~ "In situ neoplasms (D00–D09)",
    diagnosis == "E" ~ "Any endocrine, nutritional or metabolic disease",
    diagnosis == "E0007" ~ "Disorders of thyroid gland (E00–E07)",
    diagnosis == "E1014" ~ "Diabetes mellitus (E10–E14)",
    diagnosis == "E1516" ~ "Other disorders of glucose regulation and pancreatic internal secretion (E15–E16)",
    diagnosis == "E2035" ~ "Disorders of other endocrine glands (E20–E35)",
    diagnosis == "E4064" ~ "Malnutrition and other nutritional deficiencies (E40–E64)",
    diagnosis == "E6568" ~ "Obesity & other hyperalimentation (E65–E68)",
    diagnosis == "E7088" ~ "Metabolic disorders (E78–E88)",
    diagnosis == "I" ~ "Any circulatory disease",
    diagnosis == "I0509" ~ "Chronic rheumatic heart diseases (I05–I09)",
    diagnosis == "I1015" ~ "Hypertensive diseases (I10–I15)",
    diagnosis == "I2025" ~ "Ischaemic heart diseases (I20–I25)",
    diagnosis == "I2628" ~ "Pulmonary heart disease and diseases of pulmonary circulation (I26–I28)",
    diagnosis == "I3052" ~ "Other forms of heart disease (I30–I52)",
    diagnosis == "I6069" ~ "Cerebrovascular diseases (I60–I69)",
    diagnosis == "I7079" ~ "Diseases of arteries, arterioles, and capillaries (I70–I79)",
    diagnosis == "I8089" ~ "Other diseases of veins, lymphatic vessels, and lymph nodes (I80–I89)",
    diagnosis == "I9599" ~ "Other & unspecified disorders of the circulatory system (I95–I99)"),
    Diagnosis = fct_relevel(Diagnosis, diags_names_order)) |>
  arrange(Diagnosis) |>
  relocate(Diagnosis, .before = HR_CI) |>
  select(-diagnosis)

# social isolation
markov_si_hr <- lapply(markov_si_all, get_time_diag_hr)

markov_si_hr <- imap(markov_si_hr, ~.x %>%
                       mutate(diagnosis = .y)) |>
  bind_rows() |>
  mutate(Diagnosis = case_when(
    diagnosis == "C" ~ "Any neoplasm",
    diagnosis == "C0014" ~ "Malignant neoplasms of lip, oral cavity and pharynx (C00–C14)",
    diagnosis == "C1526" ~ "Malignant neoplasms of digestive organs (C15–C26)",
    diagnosis == "C3039" ~ "Malignant neoplasms of respiratory and intrathoracic organs (C30–C39)",
    diagnosis == "C4041" ~ "Malignant neoplasms of bone and articular cartilage (C40–C41)",
    diagnosis == "C4344" ~ "Melanoma & other malignant neoplasms of skin (C43–C44)",
    diagnosis == "C4549" ~ "Malignant neoplasms of mesothelial and soft tissue (C45–C49)",
    diagnosis == "C50" ~ "Malignant neoplasms of breast (C50)",
    diagnosis == "C5158" ~ "Malignant neoplasms of female genital organs (C51–C58)",
    diagnosis == "C6063" ~ "Malignant neoplasms of male genital organs (C60–C63)",
    diagnosis == "C6468" ~ "Malignant neoplasms of urinary tract (C64–C68)",
    diagnosis == "C6972" ~ "Malignant neoplasms of eye, brain and other parts of central nervous system (C69–C72)",
    diagnosis == "C7375" ~ "Malignant neoplasms of thyroid and other endocrine glands (C73–C75)",
    diagnosis == "C8196" ~ "Malignant neoplasms of lymphoid, hematopoietic and related tissue (C81–C96)",
    diagnosis == "D0009" ~ "In situ neoplasms (D00–D09)",
    diagnosis == "E" ~ "Any endocrine, nutritional or metabolic disease",
    diagnosis == "E0007" ~ "Disorders of thyroid gland (E00–E07)",
    diagnosis == "E1014" ~ "Diabetes mellitus (E10–E14)",
    diagnosis == "E1516" ~ "Other disorders of glucose regulation and pancreatic internal secretion (E15–E16)",
    diagnosis == "E2035" ~ "Disorders of other endocrine glands (E20–E35)",
    diagnosis == "E4064" ~ "Malnutrition and other nutritional deficiencies (E40–E64)",
    diagnosis == "E6568" ~ "Obesity & other hyperalimentation (E65–E68)",
    diagnosis == "E7088" ~ "Metabolic disorders (E78–E88)",
    diagnosis == "I" ~ "Any circulatory disease",
    diagnosis == "I0509" ~ "Chronic rheumatic heart diseases (I05–I09)",
    diagnosis == "I1015" ~ "Hypertensive diseases (I10–I15)",
    diagnosis == "I2025" ~ "Ischaemic heart diseases (I20–I25)",
    diagnosis == "I2628" ~ "Pulmonary heart disease and diseases of pulmonary circulation (I26–I28)",
    diagnosis == "I3052" ~ "Other forms of heart disease (I30–I52)",
    diagnosis == "I6069" ~ "Cerebrovascular diseases (I60–I69)",
    diagnosis == "I7079" ~ "Diseases of arteries, arterioles, and capillaries (I70–I79)",
    diagnosis == "I8089" ~ "Other diseases of veins, lymphatic vessels, and lymph nodes (I80–I89)",
    diagnosis == "I9599" ~ "Other & unspecified disorders of the circulatory system (I95–I99)"),
    Diagnosis = fct_relevel(Diagnosis, diags_names_order)) |>
  arrange(Diagnosis) |>
  relocate(Diagnosis, .before = HR_CI) |>
  select(-diagnosis)

markov_si_hr

sheets <- list("Loneliness" = markov_lon_hr, "Social isolation" = markov_si_hr)

write_xlsx(sheets, "table_markov_hr.xlsx")


####################################################################################

