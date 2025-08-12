
# Load packages -------
library(lm.beta)
library(dplyr)
library(broom)
library(readxl)
library(purrr)
library(tidyr)
library(tibble)
library(gtsummary)


# Key function to iterate the linear models --------
source("./scripts/recursive_lm_function.R")

# Import and merge datasets ----
pubmept1 = readRDS("./data/processed/EWAS_pubmep_t1.RDS")
pubmept2 = readRDS("./data/processed/EWAS_pubmep_t2.RDS")

iberomics = readRDS("./data/processed/EWAS_iberomics.RDS")

pubmept1 = pubmept1 %>% 
  dplyr::select(Code, TINAGL1_cg20378147:EEFSEC_cg18073874, Sex:BMI_zscore, 
                WC_zscore, DBP_zscore_Stavnsbo:SBP_zscore_Stavnsbo, 
                DBP_zscore_NHBPEP:HOMA_IR_zscore, 
                QUICKI, FMI, CD8T:Neu) %>% 
  dplyr::mutate(Cohort = 1)

pubmept2 = pubmept2 %>% 
  dplyr::select(Code, TINAGL1_cg20378147:EEFSEC_cg18073874, Sex:BMI_zscore, 
                WC_zscore, DBP_zscore_Stavnsbo:SBP_zscore_Stavnsbo, 
                DBP_zscore_NHBPEP:HOMA_IR_zscore, 
                QUICKI, FMI, CD8T:Neu) %>% 
  dplyr::mutate(Cohort = 1)

iberomics = iberomics %>% 
  dplyr::select(Code, TINAGL1_cg20378147:EEFSEC_cg18073874, Sex:BMI_zscore, 
                WC_zscore, DBP_zscore_Stavnsbo:SBP_zscore_Stavnsbo, 
                DBP_zscore_NHBPEP:HOMA_IR_zscore, 
                QUICKI, FMI, CD8T:Neu) %>% 
  dplyr::mutate(Cohort = 2)

iberomics_pre = iberomics %>% dplyr::filter(Estadio_tanner == 1)
iberomics_post = iberomics %>% dplyr::filter(Estadio_tanner != 1)
rm(iberomics)

colnames(pubmept1)[183] = "MIR1-1_cg23792592"

prepuber = rbind(pubmept1, iberomics_pre) %>% 
  dplyr::mutate(Cohort = factor(Cohort))

puber = rbind(pubmept2, iberomics_post) %>% 
  dplyr::mutate(Cohort = factor(Cohort))

rm(pubmept1, pubmept2, iberomics_pre, iberomics_post)

# Select the outcomes -----
phenotypes = prepuber %>% 
  dplyr::select(BMI_zscore, WC_zscore:SBP_zscore_Stavnsbo, DBP_zscore_NHBPEP:HOMA_IR_zscore, QUICKI, FMI) %>% 
  colnames() 

# 1) Prepubertal -------
col = colnames(prepuber)
col[2:194] = base::paste(col[2:194],"prepubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(prepuber) = col
# CpGs predictors
predictors = colnames(prepuber)[2:194] 

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = prepuber, formula = "+ Age + Sex + Cohort + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = prepuber, formula = "+ Age + Sex + Cohort + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}


# 2) Pubertal -------
col = colnames(puber)
col[2:194] = base::paste(col[2:194],"pubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(puber) = col
predictors = colnames(puber)[2:194] 

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = puber, formula = "+ Age + Sex + Cohort + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = puber, formula = "+ Age + Sex + Cohort + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}

assoc

# assoc_def = assoc %>% dplyr::filter(fdr < 0.05)

# write.csv2(assoc, "./results/Validation_EWAS/both/2024_12_04_DNAm_both_cohorts.csv", 
#            row.names = FALSE)



# 3) Prepubertal Male -------
prepuber_boys = prepuber %>% 
  dplyr::filter(Sex == 0)

col = colnames(prepuber_boys)
col[2:194] = base::paste(col[2:194],"prepubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(prepuber_boys) = col
# CpGs predictors
predictors = colnames(prepuber_boys)[2:194] 

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = prepuber_boys, formula = "+ Age + Cohort + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = prepuber_boys, formula = "+ Age + Cohort + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}

# write.csv2(assoc, "./results/Validation_EWAS/both/by_sex/2024_12_13_DNAm_both_cohorts_prepubertal_males.csv",
#            row.names = FALSE)

rm(assoc, assoc2)

# 4) Pubertal Male -------
puber_boys = puber %>% 
  dplyr::filter(Sex == 0)

col = colnames(puber_boys)
col[2:194] = base::paste(col[2:194],"pubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(puber_boys) = col
predictors = colnames(puber_boys)[2:194] 

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = puber_boys, formula = "+ Age + Cohort + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = puber_boys, formula = "+ Age + Cohort + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}



# write.csv2(assoc, "./results/Validation_EWAS/both/by_sex/2024_12_13_DNAm_both_cohorts_pubertal_males.csv",
#            row.names = FALSE)

rm(assoc, assoc2)


# 5) Prepubertal Female -------
prepuber_girls = prepuber %>% 
  dplyr::filter(Sex == 1)

col = colnames(prepuber_girls)
col[2:194] = base::paste(col[2:194],"prepubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(prepuber_girls) = col
# CpGs predictors
predictors = colnames(prepuber_girls)[2:194] 

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = prepuber_girls, formula = "+ Age + Cohort + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = prepuber_girls, formula = "+ Age + Cohort + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}

# write.csv2(assoc, "./results/Validation_EWAS/both/by_sex/2024_12_13_DNAm_both_cohorts_prepubertal_females.csv",
#            row.names = FALSE)

rm(assoc, assoc2)

# 6) Pubertal Female -------
puber_girls = puber %>% 
  dplyr::filter(Sex == 1)

col = colnames(puber_girls)
col[2:194] = base::paste(col[2:194],"pubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(puber_girls) = col
predictors = colnames(puber_girls)[2:194] 

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())

## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = puber_girls, formula = "+ Age + Cohort + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = puber_girls, formula = "+ Age + Cohort + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}


# write.csv2(assoc, "./results/Validation_EWAS/both/by_sex/2024_12_13_DNAm_both_cohorts_pubertal_females.csv",
#            row.names = FALSE)

rm(assoc, assoc2)


