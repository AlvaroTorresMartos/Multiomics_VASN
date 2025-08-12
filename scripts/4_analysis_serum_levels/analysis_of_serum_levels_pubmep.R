
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
  dplyr::select(-c(Code_new_T2))

pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(-c(Code_new_T2))



# Select the outcomes and predictors -----
phenotypes = pubmep_t1 %>% 
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
                          data = pubmep_t1, 
                          formula = "+ Age + Sex")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = pubmep_t1, 
                          formula = "+ Age + Sex + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "prepubertal", 
                cohort = "pubmep")

# write.csv(assoc, "./results/Serum_vasorin_levels/pubmep/2024_12_11_prepubertal_pubmep_association_between_serum_vasorin_levels_and_outcomes.csv", 
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
                          data = pubmep_t2,
                          formula = "+ Age + Sex")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = pubmep_t2, 
                          formula = "+ Age + Sex + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows) 
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "pubertal", 
                cohort = "pubmep")

# write.csv(assoc, "./results/Serum_vasorin_levels/pubmep/2024_12_11_pubertal_pubmep_association_between_serum_vasorin_levels_and_outcomes.csv", 
#           row.names = FALSE)



