
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

# Import datasets ----
data_t1 = readRDS("./data/processed/EWAS_pubmep_t1.RDS")
data_t2 = readRDS("./data/processed/EWAS_pubmep_t2.RDS")

# Select the outcomes -----
phenotypes = data_t1 %>% 
  dplyr::select(BMI_zscore, WC_zscore:SBP_zscore_Stavnsbo, DBP_zscore_NHBPEP:HOMA_IR_zscore, QUICKI, FMI) %>% 
  colnames() 

# Check the posibility to adjust by counfounders -----
# data_t1 %>%
#   dplyr::select(Origen, all_of(phenotypes)) %>%
#   tbl_summary(by = Origen)
# data_t2 %>%
#   dplyr::select(Origen, all_of(phenotypes)) %>%
#   tbl_summary(by = Origen)
# When we evaluate FMI: 
## In time point 1, only the children from Cordoba are regarded 
## In time point 2, the children from Cordoba and Zaragoza are regarded (DONT ADJUST BY ORIGEN)

# 1) Prepubertal -------
col = colnames(data_t1)
# dna methylation prepubertal 
col[3:195] = base::paste(col[3:195],"prepubertal", sep="_")
col[3:195] = sapply(col[3:195], function(x){
  gsub("^_cg","cg", x)
})
col[3:195] = sapply(col[3:195], function(x){
  gsub("1-1","1_1", x)
})
colnames(data_t1) = col
# CpGs predictors
predictors = colnames(data_t1)[3:195] 

assoc = data.frame(Outcome = as.character(), 
                   Predictor = as.character(), 
                   estimate = as.numeric(), 
                   std_estimate = as.numeric(), 
                   p.value = as.numeric(), 
                   all = as.character())




## Iterative linear models  ----
for (i in 1:length(phenotypes)){
  assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = data_t1, formula = "+ Age + Sex + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}




# 2) Pubertal ------

col = colnames(data_t2)
# dna methylation prepubertal 
col[3:195] = base::paste(col[3:195],"pubertal", sep="_")
col[3:195] = sapply(col[3:195], function(x){
  gsub("^_cg","cg", x)
})
col[3:195] = sapply(col[3:195], function(x){
  gsub("1-1","1_1", x)
})
colnames(data_t2) = col
# CpGs predictors
predictors = colnames(data_t2)[3:195] 


## Iterative linear models ----
for (i in 1:length(phenotypes)){
  if(i == length(phenotypes)){
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = data_t2, formula = "+ Age + Sex + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }else{
    assoc2 = lapply(predictors, recursive_lm, outcome = phenotypes[i], 
                    data = data_t2, formula = "+ Age + Sex + Origen + CD8T + CD4T + NK + Bcell + Mono + Neu" )
  }
  assoc2 = purrr::map_dfr(assoc2, ~.x, bind_rows)  %>% 
    dplyr::arrange(p.value) %>% # FDR
    dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr")) %>% 
    dplyr::relocate(fdr, .after = p.value)
  assoc = rbind(assoc, assoc2)
}


## Export results -----
# write.csv2(assoc, "./results/Validation_EWAS/pubmep/2024_12_04_pubmep_associations_between_DNA_methylation_193_CpGs_and_outcomes.csv",
#            row.names = FALSE)

# assoc = read.csv2("./results/Validation_EWAS/pubmep/2024_12_04_pubmep_associations_between_DNA_methylation_193_CpGs_and_outcomes.csv")
# 
# assoc2 = assoc %>% dplyr::filter(fdr < 0.05)
# # assoc2 = assoc %>% filter(p.value < 0.005)
# 
# outcomes = names(table(assoc2$outcome))
# predictors = names(table(assoc2$predictor))
# 
# assoc3 = assoc[assoc$outcome %in% outcomes & assoc$predictor %in% predictors, ]
# 
# 
# 
# wide_beta = assoc3 %>% select(outcome, predictor, std_coefficient)
# wide_pvalue = assoc3 %>% select(outcome, predictor, p.value)
# 
# wide_beta = tidyr::pivot_wider(wide_beta, names_from = outcome,
#                                values_from = std_coefficient) %>%
#   tibble::column_to_rownames(var = "predictor")
# 
# wide_pvalue = tidyr::pivot_wider(wide_pvalue, names_from = outcome,
#                                  values_from = p.value) %>%
#   tibble::column_to_rownames(var = "predictor")
# 
# corrplot::corrplot(data.matrix(wide_beta), method="color", sig.level =0.05,tl.cex=0.8,tl.col = "black",
#                    cl.ratio=0.4, cl.length=7, p.mat = data.matrix(wide_pvalue),
#                    insig = "blank")
# Filtro de 0.05
# text(col(data.matrix(wide_beta)), row(data.matrix(wide_beta))[36:1,], round(data.matrix(wide_pvalue), 2), cex=0.54)
# Filtro de 0.005
# text(col(data.matrix(wide_beta)), row(data.matrix(wide_beta))[15:1,], round(data.matrix(wide_pvalue), 2), cex=0.54)

# 0.05:  15 20
# Device size: 30 x 30
# Orientaiton : landscape 
# Option Use cairo marked 


