
# Load packages -------
library(lm.beta)
library(dplyr)
library(broom)
library(purrr)

# Key function to iterate the linear models --------
source("./scripts/recursive_lm_function.R")

# Import and merge datasets 
iberomics = readRDS("./data/processed/vasn_serum_levels_iberomics.RDS")

iberomics_pre = iberomics %>% 
  dplyr::filter(Estadio_tanner == 1)

iberomics_pub = iberomics %>% 
  dplyr::filter(Estadio_tanner != 1)

rm(iberomics)


# Select the outcomes and predictors -----
phenotypes = iberomics_pre %>% 
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
                          data = iberomics_pre, 
                          formula = "+ Age + Sex")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = iberomics_pre, 
                          formula = "+ Age + Sex + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "prepubertal", cohort = "iberomics")

# write.csv(assoc, "./results/Serum_vasorin_levels/iberomics/2024_12_11_prepubertal_iberomics_association_between_serum_vasorin_levels_and_outcomes.csv",
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
                          data = iberomics_pub, 
                          formula = "+ Age + Sex")
  }else{
    assoc2 = recursive_lm(outcome = phenotypes[i], predictor = "Vasn_pg_mL", 
                          data = iberomics_pub, 
                          formula = "+ Age + Sex + Origen")
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)
  assoc = rbind(assoc, assoc2)
}

assoc = assoc %>% 
  dplyr::arrange(p.value) %>% # FDR
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
  dplyr::relocate(fdr, .after = p.value) %>% 
  dplyr::mutate(time = "pubertal", cohort = "iberomics")


# write.csv(assoc, "./results/Serum_vasorin_levels/iberomics/2024_12_11_pubertal_iberomics_association_between_serum_vasorin_levels_and_outcomes.csv",
#           row.names = FALSE)


