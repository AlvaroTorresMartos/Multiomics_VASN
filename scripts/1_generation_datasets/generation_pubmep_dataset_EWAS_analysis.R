
# Load packages -------
library(dplyr)
library(readxl)
library(Hmisc)
library(tidyr)
library(tibble)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
library(stringr)

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

bmi_t1 = pubmep_t1 %>% dplyr::select(Code, Code_new_T2, BMI_Time)
bmi_t2 = pubmep_t2 %>% dplyr::select(Code, Code_new_T2, BMI_Time)

# write.csv(bmi_t1, "./data/processed/bmi_pubmep_t1.csv")
# write.csv(bmi_t2, "./data/processed/bmi_pubmep_t2.csv")

pubmep_t1 = pubmep_t1 %>% 
  dplyr::select(Code, Code_new_T2, Sex_Time, Age_Time, Origen_T1, Tanner_Time, 
                BMIz__Orbegozo__Time, WC_Time, Mgrasatg, Talla_m_Time, 
                Glucose__mg_dl__Time, Insulin__mU_l__Time, 
                HOMA_Time, QUICKI_Time, TAG__mg_dl__Time, LDLc__mg_dl__Time, 
                HDLc__mg_dl__Time, SBP_Time, DBP_Time
                )



pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(Code, Code_new_T2, Sex_Time, Age_Time, Origen_T1, Tanner_Time, 
                BMIz__Orbegozo__Time, WC_Time, Mgrasatg, Talla_m_Time, 
                Glucose__mg_dl__Time, Insulin__mU_l__Time, 
                HOMA_Time, QUICKI_Time, TAG__mg_dl__Time, LDLc__mg_dl__Time, 
                HDLc__mg_dl__Time, SBP_Time, DBP_Time
                )


colnames(pubmep_t1) = c("Code", "Code_new_T2", "Sex", "Age", "Origen", "Estadio_tanner", 
                        "BMI_zscore","Waist_circumference", "Mgrasatg", "Talla_m",  
                        "Glucose_mg_dL","Insulin_uU_ml","HOMA_IR",
                        "QUICKI","Triglycerides_mg_dL","LDLc_mg_dL","HDLc_mg_dL",
                        "SBP_mmHg", "DBP_mmHg"
)

colnames(pubmep_t2) = c("Code", "Code_new_T2", "Sex", "Age", "Origen", "Estadio_tanner", 
                        "BMI_zscore","Waist_circumference", "Mgrasatg", "Talla_m",
                        "Glucose_mg_dL","Insulin_uU_ml","HOMA_IR",
                        "QUICKI","Triglycerides_mg_dL","LDLc_mg_dL","HDLc_mg_dL",
                        "SBP_mmHg", "DBP_mmHg"
)


## Import EWAS data -----
mvalue_cpgs = read.csv2("./data/raw/1_discovery_population_pubmep/PUBMEP_mVals_BMIQ_NEW_VAL_LIST.csv", 
                        row.names = 1) %>% tibble::rownames_to_column("CpG")

lista_def2 = readxl::read_xlsx("./data/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.xlsx")
mvalue_cpgs = mvalue_cpgs %>% dplyr::filter(CpG %in% lista_def2$Name)


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[mvalue_cpgs$CpG, ]
# ann850k = ann850k[mvalue_cpgs$CpG,]
head(ann850k)
dim(ann850k)
ann850k = as.data.frame(ann850k)

mvalue_cpgs2 = mvalue_cpgs %>% dplyr::filter(CpG %in% ann850k$Name) %>% 
  tibble::column_to_rownames("CpG") %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname) %>% 
  tidyr::pivot_wider(names_from=rowname, values_from=value)

mvalue_cpgs2 = mvalue_cpgs2 %>% dplyr::rename(Code = name) 
mvalue_cpgs2$Code = stringr::str_replace(mvalue_cpgs2$Code, "\\.", "-")

## Merge -----
data_t1 = pubmep_t1 %>% dplyr::left_join(mvalue_cpgs2, by = "Code") %>% 
  dplyr::relocate(all_of(colnames(mvalue_cpgs2)[2:ncol(mvalue_cpgs2)]),
                  .after = Code_new_T2)

mvalue_cpgs2 = mvalue_cpgs2 %>% dplyr::rename(Code_new_T2 = Code) 

data_t2 = pubmep_t2 %>% dplyr::left_join(mvalue_cpgs2, by = "Code_new_T2") %>% 
  dplyr::relocate(all_of(colnames(mvalue_cpgs2)[2:ncol(mvalue_cpgs2)]),
                  .after = Code_new_T2)

## Recoding -----
data_t1 = data_t1 %>% dplyr::mutate(Sex = ifelse(Sex == "Varon_T1", 0, 1), 
                                    Sex = as.factor(Sex), 
                                    Origen = as.factor(Origen),
                                    Estadio_tanner = as.factor(Estadio_tanner)
)

data_t2 = data_t2 %>%  dplyr::mutate(Sex = ifelse(Sex == "Varon_T2", 0, 1), 
                                     Sex = as.factor(Sex), 
                                     Origen = as.factor(Origen), 
                                     Estadio_tanner = as.factor(Estadio_tanner)
)

## Rename the names of the CpGs -----
predictors = colnames(data_t1)[3:195]
names = data.frame(Name = predictors) %>% 
  dplyr::left_join(ann850k[, c("Name","UCSC_RefGene_Name")], by = "Name")
names$UCSC_RefGene_Name = stringr::str_replace(names$UCSC_RefGene_Name,";.*","")
names$def = base::paste(names$UCSC_RefGene_Name, names$Name, sep="_")
names$def = stringr::str_replace(names$def, "-", "_")
predictors = names$def
colnames(data_t1)[3:195] = predictors
rm(names)
predictors = colnames(data_t2)[3:195]
names = data.frame(Name = predictors) %>% 
  dplyr::left_join(ann850k[, c("Name","UCSC_RefGene_Name")], by = "Name")
names$UCSC_RefGene_Name = stringr::str_replace(names$UCSC_RefGene_Name,";.*","")
names$def = base::paste(names$UCSC_RefGene_Name, names$Name, sep="_")
predictors = names$def
colnames(data_t2)[3:195] = predictors
rm(names)

col = colnames(data_t2)

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

## Import and merge white cells populations -----
estimacion_celulas_blancas = read.csv2("./data/raw/1_discovery_population_pubmep/ESTIMACION_CELULAS_BLANCAS_PUBMEP.csv")
colnames(estimacion_celulas_blancas) = gsub("counts.", "", colnames(estimacion_celulas_blancas))
data_t1 = data_t1 %>% merge(estimacion_celulas_blancas, by.x = "Code", by.y = "X")
data_t2 = data_t2 %>% merge(estimacion_celulas_blancas, by.x = "Code_new_T2", by.y = "X")



## Export final datasets -----
saveRDS(data_t1, "./data/processed/EWAS_pubmep_t1.RDS")
saveRDS(data_t2, "./data/processed/EWAS_pubmep_t2.RDS")
