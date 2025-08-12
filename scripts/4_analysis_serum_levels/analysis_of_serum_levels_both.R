
# Load packages -------
library(lm.beta)
library(dplyr)
library(broom)
library(purrr)

# Key function to iterate the linear models --------
source("./scripts/recursive_lm_function.R")

# Import and merge datasets 
pubmep_t1 = readRDS("./data/processed/vasn_serum_levels_pubmep_t1.RDS")
pubmep_t2 = readRDS("./data/processed/vasn_serum_levels_pubmep_t2.RDS")

pubmep_t1 = pubmep_t1 %>% 
  dplyr::select(-c(Code_new_T2)) %>% 
  dplyr::mutate(Cohort = 1)

pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(-c(Code_new_T2)) %>% 
  dplyr::mutate(Cohort = 1)

iberomics = readRDS("./data/processed/vasn_serum_levels_iberomics.RDS")

iberomics_pre = iberomics %>% 
  dplyr::filter(Estadio_tanner == 1) %>% 
  dplyr::mutate(Cohort = 2)

iberomics_pub = iberomics %>% 
  dplyr::filter(Estadio_tanner != 1) %>% 
  dplyr::mutate(Cohort = 2)

rm(iberomics)

prepuber = rbind(pubmep_t1, iberomics_pre) %>% 
  dplyr::mutate(Cohort = as.factor(Cohort))

puber = rbind(pubmep_t2, iberomics_pub) %>% 
  dplyr::mutate(Cohort = as.factor(Cohort))

rm(pubmep_t1, pubmep_t2, iberomics_pre, iberomics_pub)


# Select the outcomes and predictors -----
phenotypes = prepuber %>% 
  dplyr::select(BMI_zscore, WC_zscore,  
                DBP_zscore_Stavnsbo:SBP_zscore_Stavnsbo, 
                DBP_zscore_NHBPEP:BP_zscore_Stavnsbo, 
                HDLc_zscore, HOMA_IR_zscore, 
                FMI
                ) %>% 
  colnames() 

# 1) Prepubertal -------
assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = prepuber, 
                          formula = "+ Age + Sex + Cohort")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = prepuber, 
                          formula = "+ Age + Sex + Cohort + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  
  assoc = rbind(assoc, assoc2) 
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "prepubertal", cohort = "both")

# write.csv(assoc, "./results/Serum_vasorin_levels/both/2024_12_11_prepubertal_both_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)

rm(assoc, assoc2)

# 2) Pubertal -------
assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = puber, 
                          formula = "+ Age + Sex + Cohort")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = puber, 
                          formula = "+ Age + Sex + Cohort + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "pubertal", cohort = "both")


# write.csv(assoc, "./results/Serum_vasorin_levels/both/2024_12_11_pubertal_both_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)

rm(assoc, assoc2)

# 3) Prepubertal Male -------
prepuber_boys = prepuber %>% 
  dplyr::filter(Sex == 0)

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = prepuber_boys, 
                          formula = "+ Age + Cohort")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = prepuber_boys, 
                          formula = "+ Age + Cohort + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  
  assoc = rbind(assoc, assoc2) 
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "prepubertal", cohort = "both", sex = "male")

# write.csv(assoc, "./results/Serum_vasorin_levels/both/by_sex/2024_12_11_prepubertal_males_both_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)

rm(assoc, assoc2)

# 4) Pubertal Male -------
puber_boys = puber %>% 
  dplyr::filter(Sex == 0)

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())



## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = puber_boys, 
                          formula = "+ Age + Cohort")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = puber_boys, 
                          formula = "+ Age + Cohort + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "pubertal", cohort = "both", sex = "male")


# write.csv(assoc, "./results/Serum_vasorin_levels/both/by_sex/2024_12_11_pubertal_males_both_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)

rm(assoc, assoc2)


# 5) Prepubertal Female -------
prepuber_girls = prepuber %>% 
  dplyr::filter(Sex == 1)

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = prepuber_girls, 
                          formula = "+ Age + Cohort")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = prepuber_girls, 
                          formula = "+ Age + Cohort + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  
  assoc = rbind(assoc, assoc2) 
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "prepubertal", cohort = "both", sex = "female")

# write.csv(assoc, "./results/Serum_vasorin_levels/both/by_sex/2024_12_11_prepubertal_females_both_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)

rm(assoc, assoc2)

# 6) Pubertal Female -------
puber_girls = puber %>% 
  dplyr::filter(Sex == 1)

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = puber_girls, 
                          formula = "+ Age + Cohort")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = puber_girls, 
                          formula = "+ Age + Cohort + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "pubertal", cohort = "both", sex = "female")


# write.csv(assoc, "./results/Serum_vasorin_levels/both/by_sex/2024_12_11_pubertal_females_both_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)

rm(assoc, assoc2)
