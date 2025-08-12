
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

# Import datasets  ----
data = readRDS("./data/processed/EWAS_iberomics.RDS")

# Select the outcomes -----
phenotypes = data %>% 
  dplyr::select(BMI_zscore:WC_zscore, QUICKI, DBP_zscore_Stavnsbo:SBP_zscore_Stavnsbo, 
                DBP_zscore_NHBPEP:FMI) %>% 
  colnames() 

# Check the posibility to adjust by counfounders -----
# data %>%
#   dplyr::filter(Estadio_tanner == 1) %>% 
#   dplyr::select(Origen, all_of(phenotypes)) %>%
#   tbl_summary(by = Origen)
# data %>%
#   dplyr::filter(Estadio_tanner != 1) %>% 
#   dplyr::select(Origen, all_of(phenotypes)) %>%
#   tbl_summary(by = Origen)
# When we evaluate FMI: 
## The prepubertal children from Zaragoza haven't got FMI (DONT ADJUST by ORIGEN)
## The pubertal children from Zaragoza haven't got FMI (DONT ADJUST by ORIGEN)



## 1) Prepubertal -------
data2 = data %>% dplyr::filter(Estadio_tanner == 1) 
col = colnames(data2)
# dna methylation prepubertal 
col[2:194] = base::paste(col[2:194],"prepubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(data2) = col
# CpGs predictors
predictors = colnames(data2)[2:194] 

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
                    data = data2, formula = "+ Age + Sex + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = data2, formula = "+ Age + Sex + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}




## 2) Pubertal -------
data2 = data %>% dplyr::filter(Estadio_tanner != 1) 
col = colnames(data2)
col[2:194] = base::paste(col[2:194],"pubertal", sep="_")
col[2:194] = sapply(col[2:194], function(x){
  gsub("^_cg","cg", x)
})
col[2:194] = sapply(col[2:194], function(x){
  gsub("1-1","1_1", x)
})
colnames(data2) = col
predictors = colnames(data2)[2:194] 



## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = data2, formula = "+ Age + Sex + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = data2, formula = "+ Age + Sex + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}

## Export results -----
write.csv2(assoc, "./results/Validation_EWAS/iberomics/2024_12_04_iberomics_associations_between_DNA_methylation_193_CpGs_and_outcomes.csv",
           row.names = FALSE)


