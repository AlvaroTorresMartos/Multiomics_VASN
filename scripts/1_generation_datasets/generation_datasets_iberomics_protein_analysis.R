
# Load packages -------
library(dplyr)

# Import databases, merge and preprocess ------
iberomics = readRDS("/home/usuario/Escritorio/CIBM/BASES_IBEROMICS/OMICS/OMICS_BD/2023_12_19_Base_integrada_OMICS_def.RDS")
iberomics2 = readRDS("/home/usuario/Escritorio/CIBM/BASES_IBEROMICS/OMICS/OMICS_BD/2024_02_14_Base_integrada_OMICS_def.RDS")

data = iberomics %>% 
  dplyr::select(Code, Sex, Age, Origen, Estadio_tanner, 
                BMI_zscore_Orbegozo_longi, Perimetro_de_cintura,
                WC_zscore_Ferrandez, Height,
                Glucose__mg_dl_, Insulin__mU_ml_, HOMA_IR, 
                QUICKI, TAG__mg_dl_, LDLc__mg_dl_, HDLc__mg_dl_, TS, TD, 
                VASN_pg_mL
  )

colnames(data) = c("Code", "Sex", "Age", "Origen", "Estadio_tanner", 
                   "BMI_zscore", "Waist_circumference", 
                   "WC_zscore", "Talla_cm",
                   "Glucose_mg_dL","Insulin_uU_ml","HOMA_IR",
                   "QUICKI","Triglycerides_mg_dL","LDLc_mg_dL","HDLc_mg_dL",
                   "SBP_mmHg", "DBP_mmHg", "Vasn_pg_mL"
)

data = data %>% 
  dplyr::mutate(log_vasn_pg_mL = log(Vasn_pg_mL))

dxa = iberomics2 %>% dplyr::select(Code, FMI_dexa) %>% 
  dplyr::rename(FMI = FMI_dexa)

data = data %>% 
  merge(dxa, by = "Code") %>% 
  dplyr::mutate(Sex = ifelse(Sex == 1, 0, 1),
                Sex = as.factor(Sex), 
                Origen = as.factor(Origen), 
                Estadio_tanner = as.factor(Estadio_tanner)) %>% 
  dplyr::mutate(Talla_m = Talla_cm/100, .keep = "unused") %>% 
  dplyr::relocate(Talla_m, .after = "Waist_circumference") %>% 
  dplyr::relocate(WC_zscore, .after = "FMI")

## Import functions and reference values from ObMetrics -----
source("./scripts/functions_referencevalues.R")

## Calculate ObMetrics outputs -----
data = data %>% 
  # calculate_wc_zscore(ref_wc_f = wc_ferrandez_score_girls, 
  #                     ref_wc_m = wc_ferrandez_score_boys) %>% 
  dbp_zscore_stravnsbo() %>% sbp_zscore_stravnsbo() %>% 
  calculate_height_zscore %>% calculate_bp_zscore() %>%  
  bp_zscore() %>% 
  tag_zscore_stravnsbo() %>% hdl_zscore_stravnsbo() %>%
  ldl_zscore_stravnsbo() %>% glu_zscore_stravnsbo() %>% 
  insulin_zscore_stravnsbo() %>% homa_zscore_stravnsbo()

  

## Export final datasets -----
# saveRDS(data, "./data/processed/vasn_serum_levels_iberomics.RDS")
