##############################
#### MULTISTATE FUNCTIONS ####
##############################

# Code for paper by Komulainen et al. (2025)
# Loneliness and social isolation in transitions to adverse health conditions and mortality: an analysis of data from the UK Biobank study

# Code written by Ripsa Niemi & Mai Gutvilig
# 22.04.2025

# This file includes all functions and variable definitions needed to run the analyses.

#########################################

# State-matrix 
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, 3] <- 3
dimnames(tmat) <- list(from = c("Healthy", "Morbidity", "Mortality"), 
                       to = c("Healthy", "Morbidity", "Mortality"))

################ Endocrine and circulatory diseases  (E & I) #######

############## E & I diagnoses to list ###########

dgs <-  c( "E", "E0007", "E1014", "E1516", "E2035", "E4064", "E6568", "E7088", 
           "I", "I0509", "I1015", "I2025", "I2628", "I3052", "I6069", "I7079", "I8089", "I9599")

diags <- paste0(dgs, "_true")
time_diags <- paste0("time_", dgs)

list_EI <- list(diags, time_diags)
list_EI2 <- list()

for(i in 1:length(dgs)) {
  list_EI2[[dgs[i]]] <- c(list_EI[[1]][i], list_EI[[2]][i])
}

################ cancers C to list ###################

cans <- c( "C", "C0014", "C1526", "C3039", "C4041", "C4344", "C4549", "C6468",
          "C6972", "C7375", "C8196", "D0009")

cancers <- paste0(cans, "_true")
time_cancers <- paste0("time_", cans)

list_C <- list(cancers, time_cancers)
list_C2 <- list()

for(i in 1:length(cans)) {
  list_C2[[cans[i]]] <- c(list_C[[1]][i], list_C[[2]][i])
}


################ sex-specific cancers #################

cans_sex <- c("C50", "C5158", "C6063")

cancers_sex <- paste0(cans_sex, "_true")
time_cancers_sex <- paste0("time_", cans_sex)

list_C_sex <- list(cancers_sex, time_cancers_sex)
list_C2_sex <- list()

for(i in 1:length(cans_sex)) {
  list_C2_sex[[cans_sex[i]]] <- c(list_C_sex[[1]][i], list_C_sex[[2]][i])
}

#####################################################
# long names of diagnoses in order
diags_names_order <- c("Any neoplasm", 
                       "Malignant neoplasms of lip, oral cavity and pharynx (C00–C14)", 
                       "Malignant neoplasms of digestive organs (C15–C26)", 
                       "Malignant neoplasms of respiratory and intrathoracic organs (C30–C39)", 
                       "Malignant neoplasms of bone and articular cartilage (C40–C41)", 
                       "Melanoma & other malignant neoplasms of skin (C43–C44)", 
                       "Malignant neoplasms of mesothelial and soft tissue (C45–C49)", 
                       "Malignant neoplasms of breast (C50)", 
                       "Malignant neoplasms of female genital organs (C51–C58)", 
                       "Malignant neoplasms of male genital organs (C60–C63)", 
                       "Malignant neoplasms of urinary tract (C64–C68)", 
                       "Malignant neoplasms of eye, brain and other parts of central nervous system (C69–C72)", 
                       "Malignant neoplasms of thyroid and other endocrine glands (C73–C75)", 
                       "Malignant neoplasms of lymphoid, hematopoietic and related tissue (C81–C96)", 
                       "In situ neoplasms (D00–D09)", 
                       "Any endocrine, nutritional or metabolic disease", 
                       "Disorders of thyroid gland (E00–E07)", 
                       "Diabetes mellitus (E10–E14)", 
                       "Other disorders of glucose regulation and pancreatic internal secretion (E15–E16)", 
                       "Disorders of other endocrine glands (E20–E35)", 
                       "Malnutrition and other nutritional deficiencies (E40–E64)", 
                       "Obesity & other hyperalimentation (E65–E68)", 
                       "Metabolic disorders (E78–E88)", 
                       "Any circulatory disease", 
                       "Chronic rheumatic heart diseases (I05–I09)", 
                       "Hypertensive diseases (I10–I15)", 
                       "Ischaemic heart diseases (I20–I25)",
                       "Pulmonary heart disease and diseases of pulmonary circulation (I26–I28)", 
                       "Other forms of heart disease (I30–I52)", 
                       "Cerebrovascular diseases (I60–I69)", 
                       "Diseases of arteries, arterioles, and capillaries (I70–I79)", 
                       "Other diseases of veins, lymphatic vessels, and lymph nodes (I80–I89)", 
                       "Other & unspecified disorders of the circulatory system (I95–I99)")

# short names of diagnoses in order
diags_order <- c("C", "C0014", "C1526", "C3039", "C4041","C4344", "C4549","C50","C5158",
                 "C6063", "C6468", "C6972", "C7375", "C8196", "D0009",
                 "E" , "E0007", "E1014", "E1516", "E2035", "E4064", "E6568", "E7088",
                 "I", "I0509", "I1015", "I2025", "I2628", "I3052", "I6069", "I7079",
                 "I8089", "I9599")

dg_keys <- as.list(diags_names_order) |> setNames(diags_order)

########################

# Model covariates

covs_lon <- c("loneliness",
              "sex", "white_1", "white_999", "age", 
              "depgroup_1", "depgroup_2", "depgroup_999",
              "edugroup_1", "edugroup_2", "edugroup_999")

covs_si <- covs_lon
covs_si[1] <- "isolate_bin"

# dataframe for estimating baseline hazards
newd_lon <- data.frame(loneliness = rep(0, 3), age = rep(0, 3), sex = rep(0,3),
                       white_1= rep(0, 3), white_999= rep(0, 3),
                       edugroup_1 = rep(0,3), edugroup_2 = rep(0,3), edugroup_999 = rep(0,3),
                       depgroup_1 = rep(0,3), depgroup_2 = rep(0,3), depgroup_999 = rep(0,3),
                       death= rep(0, 3), trans = 1:3, 
                       time_C = rep(0, 3), time_C0014 = rep(0, 3), time_C1526 = rep(0, 3), 
                       time_C3039 = rep(0, 3), time_C4041 = rep(0, 3), time_C4344 = rep(0, 3),
                       time_C4549 = rep(0, 3), time_C50 = rep(0, 3), time_C5158 = rep(0, 3),
                       time_C6063 = rep(0, 3), time_C6468 = rep(0, 3), time_C6972 = rep(0, 3), 
                       time_C7375 = rep(0, 3),  time_C8196 = rep(0, 3), time_D0009 = rep(0, 3),
                       time_E = rep(0, 3), time_E0007 = rep(0, 3), time_E1014 = rep(0, 3),
                       time_E1516 = rep(0, 3), time_E2035 = rep(0, 3), time_E4064 = rep(0, 3),
                       time_E6568 = rep(0, 3), time_E7088 = rep(0, 3),
                       time_I = rep(0, 3), time_I0509 = rep(0, 3), time_I1015 = rep(0, 3), 
                       time_I2025 = rep(0, 3), time_I2628 = rep(0, 3), time_I3052 = rep(0, 3),
                       time_I6069 = rep(0, 3), time_I7079 = rep(0, 3), time_I8089 = rep(0, 3), 
                       time_I9599 = rep(0, 3))

newd_si <- newd_lon |> rename(isolate_bin = loneliness)

#########################################
######### ANALYSIS FUNCTIONS ############

# Data prepping function
data_prep_func <- function(diag, time_var, death_var, data) {
  
  # Data according to current diagnosis
  dat_diag <- data %>% 
    select(diag[1], diag[2], loneliness, isolate_bin, death, time, death_c, time_c, sex, age, white, edugroup, depgroup) %>% 
    mutate(death2 = get(death_var),
           time2 = get(time_var)) %>%
    select(-time, -death, -death_c, -time_c) %>%
    rename(death = death2,
           time = time2) %>%
    drop_na()  %>%
    dummy_cols(select_columns = "white") %>%
    dummy_cols(select_columns = "edugroup") %>%
    dummy_cols(select_columns = "depgroup")
  
  # which diagnosis is running
  print(diag[1])
  
  # covariates
  covs <- c(diag[2], "sex", 
            "white_1", "white_999", 
            "age", "loneliness", "isolate_bin",
            "edugroup_1", "edugroup_2", "edugroup_999",
            "depgroup_1", "depgroup_2", "depgroup_999")
  
  # multistate prep function
  dat_prep <- msprep(time = c(NA, diag[2], "time"), 
                     status = c(NA, diag[1], "death"), 
                     data = dat_diag, 
                     trans = tmat, 
                     keep = covs)
  
  # expand dataset
  dat_ms <- expand.covs(dat_prep, covs=covs, append = TRUE)
  
  # diagnosis 0/1 variable to third transition
  dat_ms$diagnosis <- 0
  dat_ms$diagnosis[dat_ms$trans==3] <- 1
  
  # return dataset
  dat_ms
  
}

#################################
########## MAIN MODELS ##########
#################################

# Loneliness
ms_analysis_lon <- function(diag, time_var, death_var) {
  
  # run data prep function
  dat_ms <- data_prep_func(diag, time_var, death_var, dat)

  # MODEL 0: stratified baseline hazards
  m0 <- coxph(Surv(Tstart, Tstop, status) ~  loneliness.1 + loneliness.2 + loneliness.3 +
                sex.1 + sex.2 + sex.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                strata(trans), data=dat_ms)
  
  # MODEL 1: proportional baseline hazards
  m1 <- coxph(Surv(Tstart, Tstop, status) ~  loneliness.1 + loneliness.2 + loneliness.3 +
                sex.1 + sex.2 + sex.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + strata(to), data=dat_ms)
  
  models <- list(m0 = m0, m1 = m1)
  models
  
}

# Sex-specific loneliness
ms_analysis_lon_sex <- function(diag, time_var, death_var, sex_val) {
  
  # keep only wanted sex
  data_sex <- dat %>% 
    filter(sex == sex_val)
  
  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, data_sex)
  
  #########
  # MODEL 0: stratified baseline hazards
  m0 <- coxph(Surv(Tstart, Tstop, status) ~  loneliness.1 + loneliness.2 + loneliness.3 + 
                white_1.1 + white_1.2 + white_1.3 + 
                white_999.1 + white_999.2 + white_999.3 + 
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 + 
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 + 
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 + 
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                strata(trans), data=dat_ms)
  
  # MODEL 1: proportional baseline hazards
  m1 <- coxph(Surv(Tstart, Tstop, status) ~  loneliness.1 + loneliness.2 + loneliness.3 + 
                white_1.1 + white_1.2 + white_1.3 + 
                white_999.1 + white_999.2 + white_999.3 + 
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 + 
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 + 
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 + 
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + strata(to), data=dat_ms)

  models <- list(m0 = m0, m1 = m1)
  models
}

# Social isolation
ms_analysis_si <- function(diag, time_var, death_var) {

  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, dat)
  
  # MODEL 0: stratified baseline hazards
  m0 <- coxph(Surv(Tstart, Tstop, status) ~  isolate_bin.1 + isolate_bin.2 + isolate_bin.3 +
                sex.1 + sex.2 + sex.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                strata(trans), data=dat_ms)
  
  # MODEL 1: proportional baseline hazards
  m1 <- coxph(Surv(Tstart, Tstop, status) ~  isolate_bin.1 + isolate_bin.2 + isolate_bin.3 +
                sex.1 + sex.2 + sex.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + strata(to), data=dat_ms)

  models <- list(m0 = m0, m1 = m1)
  models
}

# Sex-specific social isolation
ms_analysis_si_sex <- function(diag, time_var, death_var, sex_val) {
  
  # keep only wanted sex
  data_sex <- dat %>% 
    filter(sex == sex_val)
  
  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, data_sex)
  
  # MODEL 0: stratified baseline hazards
  m0 <- coxph(Surv(Tstart, Tstop, status) ~  isolate_bin.1 + isolate_bin.2 + isolate_bin.3 + 
                white_1.1 + white_1.2 + white_1.3 + 
                white_999.1 + white_999.2 + white_999.3 + 
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 + 
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 + 
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 + 
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                strata(trans), data=dat_ms)
  
  # MODEL 1: proportional baseline hazards
  m1 <- coxph(Surv(Tstart, Tstop, status) ~  isolate_bin.1 + isolate_bin.2 + isolate_bin.3 + 
                white_1.1 + white_1.2 + white_1.3 + 
                white_999.1 + white_999.2 + white_999.3 + 
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 + 
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 + 
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 + 
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + strata(to), data=dat_ms)

  models <- list(m0 = m0, m1 = m1)
  models
}



#############################################
#### MARKOV ASSUMPTION DIAGNOSTIC MODELS ####
#############################################

# Loneliness
ms_analysis_lon_markov <- function(diag, time_var, death_var) {
  
  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, dat)
  
  # time to diagnosis variable
  time_dg <- paste0("time_", str_remove(diag[1], "_true"), ".3")
  dat_ms <- dat_ms %>% mutate(time_diag.3 = get(time_dg))
  
  # MODEL 2: time to diagnosis (for trans 3) included in the model
  m2 <- coxph(Surv(Tstart, Tstop, status) ~  loneliness.1 + loneliness.2 + loneliness.3 +
                sex.1 + sex.2 + sex.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + time_diag.3 + strata(to), data=dat_ms)

  models <- list(m2 = m2)
  models
  
}

# Sex-specific loneliness
ms_analysis_lon_markov_sex <- function(diag, time_var, death_var, sex_val) {
  
  # keep only wanted sex
  data_sex <- dat %>% 
    filter(sex == sex_val)
  
  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, data_sex)
  
  # time to diagnosis variable
  time_dg <- paste0("time_", str_remove(diag[1], "_true"), ".3")
  dat_ms <- dat_ms %>% mutate(time_diag.3 = get(time_dg))
  
  # MODEL 2: time to diagnosis (for trans 3) included in the model
  m2 <- coxph(Surv(Tstart, Tstop, status) ~  loneliness.1 + loneliness.2 + loneliness.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + time_diag.3 + strata(to), data=dat_ms)
  
  models <- list(m2 = m2)
  models
  
}

# Social isolation
ms_analysis_si_markov <- function(diag, time_var, death_var) {
  
  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, dat)

  # time to diagnosis variable
  time_dg <- paste0("time_", str_remove(diag[1], "_true"), ".3")
  dat_ms <- dat_ms %>% mutate(time_diag.3 = get(time_dg))
  
  # MODEL 2: time to diagnosis (for trans 3) included in the model
  m2 <- coxph(Surv(Tstart, Tstop, status) ~  isolate_bin.1 + isolate_bin.2 + isolate_bin.3 +
                sex.1 + sex.2 + sex.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + time_diag.3 + strata(to), data=dat_ms)
  
  models <- list(m2 = m2)
  models
}

# Sex-specific social isolation 
ms_analysis_si_markov_sex <- function(diag, time_var, death_var, sex_val) {
  
  # keep only wanted sex
  data_sex <- dat %>% 
    filter(sex == sex_val)
  
  # prep the data
  dat_ms <- data_prep_func(diag, time_var, death_var, data_sex)

  # time to diagnosis variable
  time_dg <- paste0("time_", str_remove(diag[1], "_true"), ".3")
  dat_ms <- dat_ms %>% mutate(time_diag.3 = get(time_dg))
  
  # MODEL 2: time to diagnosis (for trans 3) included in the model
  m2 <- coxph(Surv(Tstart, Tstop, status) ~  isolate_bin.1 + isolate_bin.2 + isolate_bin.3 +
                white_1.1 + white_1.2 + white_1.3 +
                white_999.1 + white_999.2 + white_999.3 +
                age.1 + age.2 + age.3 +
                edugroup_1.1 + edugroup_1.2 + edugroup_1.3 +
                edugroup_2.1 + edugroup_2.2 + edugroup_2.3 +
                edugroup_999.1 + edugroup_999.2 + edugroup_999.3 +
                depgroup_1.1 + depgroup_1.2 + depgroup_1.3 +
                depgroup_2.1 + depgroup_2.2 + depgroup_2.3 +
                depgroup_999.1 + depgroup_999.2 + depgroup_999.3 +
                diagnosis + time_diag.3 + strata(to), data=dat_ms)
  
  models <- list(m2 = m2)
  models
}

########################

# Result table modification functions

# take loneliness/isolation estimates
mod_table <- function(tab, model) {
  
  m1_diag <- tab[[model]] |> 
    tidy(conf.int = TRUE, exponentiate = TRUE) |>
    select(term, estimate, std.error, starts_with("conf"), p.value) |>
    slice(1:3) %>%
    mutate(trans=c(1:3))
}

# combine diagnoses to same table
res_list_to_tab <- function(tab_list) {
  imap(tab_list, ~.x %>%
         mutate(diagnosis = .y)) |>
    bind_rows() |>
    relocate(diagnosis, .before = term) |>
    relocate(trans, .before=term) %>% 
    mutate(p.value = round(p.value, 5)) %>%
    select(-term) 
}


# function to get markov assumption diagnostic model estimate
get_time_diag_hr <- function(model) {
  
  hr <- model[["m2"]] |> 
    tidy(conf.int = TRUE, exponentiate = TRUE) |> 
    filter(term == "time_diag.3") |>
    mutate(HR_CI = paste0(round(estimate, 3), " (", round(conf.low, 3), "–", round(conf.high, 3), ")")) |>
    select(HR_CI)
}

