
# Load packages -------
library(dplyr)
library(readxl)
library(Hmisc)

# Import databases, merge and preprocess ------
pubmep_long = read.csv2("./data/raw/1_discovery_population_pubmep/BASE_PUBMEP_LONGITUDINAL_LONGformat_NIÑAS_Y_NIÑOS_213inds_31_01_2022_EPIC.csv",
                        header=TRUE, sep=";", stringsAsFactors=F, dec=",", 
                        na.strings=c(""," ","NaN","NA"))

## Split prepubertal and pubertal children -----
pubmep_t1 = pubmep_long %>% dplyr::filter(time == 1)
pubmep_t2 = pubmep_long %>% dplyr::filter(time == 2)
# vars = data.frame(colnames(pubmep_long))

## Import DXA information -----
pubmep_dexa = readxl::read_xls("./data/raw/1_discovery_population_pubmep/BasePUBMEPSantiago_DEXA_PUBMEP.xls", 
                               sheet = 2) 
pubmep_dexa =  pubmep_dexa %>%
  dplyr::select(Code_new, Mgrasatg) %>% 
  dplyr::rename(Code_new_T2 = Code_new) 

genobox_dexa = readRDS("./data/raw/1_discovery_population_pubmep/2023_03_01_GENOBOX_BASE_COMPLETA_INTEGRADA.RDS")
genobox_dexa = genobox_dexa %>% 
  dplyr::select(4, Mgrasatg) %>% 
  dplyr::filter(Code %nin% "R15604") 


pubmep_t1 = dplyr::left_join(pubmep_t1, genobox_dexa, by = "Code")
pubmep_t2 = dplyr::left_join(pubmep_t2, pubmep_dexa, by = "Code_new_T2")

pubmep_t1 = pubmep_t1 %>% 
  dplyr::select(Code, Code_new_T2, Sex_Time, Age_Time, Origen_T1, Tanner_Time, 
                BMIz__Orbegozo__Time, WC_Time, Mgrasatg, Talla_m_Time, 
                Glucose__mg_dl__Time, Insulin__mU_l__Time, 
                HOMA_Time, QUICKI_Time, TAG__mg_dl__Time, LDLc__mg_dl__Time, 
                HDLc__mg_dl__Time, SBP_Time, DBP_Time, VASN_T1
  )

pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(Code, Code_new_T2, Sex_Time, Age_Time, Origen_T1, Tanner_Time, 
                BMIz__Orbegozo__Time, WC_Time, Mgrasatg, Talla_m_Time, 
                Glucose__mg_dl__Time, Insulin__mU_l__Time, 
                HOMA_Time, QUICKI_Time, TAG__mg_dl__Time, LDLc__mg_dl__Time, 
                HDLc__mg_dl__Time, SBP_Time, DBP_Time, VASN_T2
  )

colnames(pubmep_t1) = c("Code", "Code_new_T2", "Sex", "Age", "Origen", "Estadio_tanner", 
                        "BMI_zscore","Waist_circumference", "Mgrasatg", "Talla_m",  
                        "Glucose_mg_dL","Insulin_uU_ml","HOMA_IR",
                        "QUICKI","Triglycerides_mg_dL","LDLc_mg_dL","HDLc_mg_dL",
                        "SBP_mmHg", "DBP_mmHg", "Vasn_pg_mL"
)

colnames(pubmep_t2) = c("Code", "Code_new_T2", "Sex", "Age", "Origen", "Estadio_tanner", 
                        "BMI_zscore","Waist_circumference", "Mgrasatg", "Talla_m",
                        "Glucose_mg_dL","Insulin_uU_ml","HOMA_IR",
                        "QUICKI","Triglycerides_mg_dL","LDLc_mg_dL","HDLc_mg_dL",
                        "SBP_mmHg", "DBP_mmHg", "Vasn_pg_mL"
)

## Recoding -----
data_t1 = pubmep_t1 %>% dplyr::mutate(Sex = ifelse(Sex == "Varon_T1", 0, 1), 
                                    Sex = as.factor(Sex), 
                                    Origen = as.factor(Origen),
                                    Estadio_tanner = as.factor(Estadio_tanner), 
                                    log_vasn_pg_mL = log(Vasn_pg_mL)
)

data_t2 = pubmep_t2 %>%  dplyr::mutate(Sex = ifelse(Sex == "Varon_T2", 0, 1), 
                                     Sex = as.factor(Sex), 
                                     Origen = as.factor(Origen), 
                                     Estadio_tanner = as.factor(Estadio_tanner), 
                                     log_vasn_pg_mL = log(Vasn_pg_mL)
)


## Calculate FMI -----
data_t1 = data_t1 %>% 
  dplyr::mutate(FMI = (Mgrasatg/1000)/((Talla_m)**2)) %>% 
  dplyr::select(-c(Mgrasatg))
data_t2 = data_t2 %>% 
  dplyr::mutate(FMI = (Mgrasatg/1000)/((Talla_m)**2)) %>% 
  dplyr::select(-c(Mgrasatg))

## Import functions and reference values from ObMetrics -----
source("./scripts/functions_referencevalues.R")

## Calculate ObMetrics outputs -----
data_t1 = data_t1 %>% 
  calculate_wc_zscore(ref_wc_f = wc_ferrandez_score_girls, 
                      ref_wc_m = wc_ferrandez_score_boys) %>% 
  dbp_zscore_stravnsbo() %>% sbp_zscore_stravnsbo() %>% 
  calculate_height_zscore %>% calculate_bp_zscore() %>%  
  bp_zscore() %>% 
  tag_zscore_stravnsbo() %>% hdl_zscore_stravnsbo() %>%
  ldl_zscore_stravnsbo() %>% glu_zscore_stravnsbo() %>% 
  insulin_zscore_stravnsbo() %>% homa_zscore_stravnsbo() 

data_t2 = data_t2 %>% 
  calculate_wc_zscore(ref_wc_f = wc_ferrandez_score_girls, 
                      ref_wc_m = wc_ferrandez_score_boys) %>% 
  dbp_zscore_stravnsbo() %>% sbp_zscore_stravnsbo() %>% 
  calculate_height_zscore %>% calculate_bp_zscore() %>%  
  bp_zscore() %>% 
  tag_zscore_stravnsbo() %>% hdl_zscore_stravnsbo() %>%
  ldl_zscore_stravnsbo() %>% glu_zscore_stravnsbo() %>% 
  insulin_zscore_stravnsbo() %>% homa_zscore_stravnsbo() 

## Export final datasets -----
# saveRDS(data_t1, "./data/processed/vasn_serum_levels_pubmep_t1.RDS")
# saveRDS(data_t2, "./data/processed/vasn_serum_levels_pubmep_t2.RDS")
