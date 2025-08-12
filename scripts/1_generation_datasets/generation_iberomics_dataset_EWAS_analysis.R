
# Load packages -------
library(dplyr)
library(readxl)
library(Hmisc)
library(tidyr)
library(tibble)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
library(stringr)

# Import databases, merge and preprocess ------
iberomics = readRDS("/home/usuario/Escritorio/CIBM/BASES_IBEROMICS/OMICS/OMICS_BD/2023_12_19_Base_integrada_OMICS_def.RDS")
# iberomics2 = readRDS("/home/usuario/Escritorio/CIBM/BASES_IBEROMICS/OMICS/OMICS_BD/2024_02_14_Base_integrada_OMICS_def.RDS")

# bmi_ibe = iberomics %>% 
#   dplyr::select(Code, BMI)
# 
# write.csv(bmi_ibe, "./data/processed/bmi_iberomics.csv")

data = iberomics %>% 
  dplyr::select(Code, Sex, Age, Origen, Estadio_tanner, 
                BMI_zscore_Orbegozo_longi, WC_zscore_Ferrandez, Height,
                Glucose__mg_dl_, Insulin__mU_ml_, HOMA_IR, 
                QUICKI, TAG__mg_dl_, LDLc__mg_dl_, HDLc__mg_dl_, TS, TD
  )

colnames(data) = c("Code", "Sex", "Age", "Origen", "Estadio_tanner", 
                   "BMI_zscore","WC_zscore", "Talla_cm",
                   "Glucose_mg_dL","Insulin_uU_ml","HOMA_IR",
                   "QUICKI","Triglycerides_mg_dL","LDLc_mg_dL","HDLc_mg_dL",
                   "SBP_mmHg", "DBP_mmHg"
)


## Import EWAS data -----
lista_def2 = readxl::read_xlsx("./data/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.xlsx")
mvalue_cpgs = read.csv2("./data/raw/2_validation_population_iberomics/mVals_BMIQ_NEW_VAL_LIST.csv", 
                        row.names = 1) %>% tibble::rownames_to_column("CpG") %>% 
  dplyr::filter(CpG %in% lista_def2$Name)

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)
ann850k = as.data.frame(ann850k)
ann850k = ann850k[mvalue_cpgs$CpG,]
head(ann850k)
dim(ann850k)

mvalue_cpgs2 = mvalue_cpgs %>% dplyr::filter(CpG %in% ann850k$Name) %>% 
  tibble::column_to_rownames("CpG") %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname) %>% 
  tidyr::pivot_wider(names_from=rowname, values_from=value)

## Merge -----
data = mvalue_cpgs2 %>% dplyr::rename(Code = name) %>%
  dplyr::left_join(data, by = "Code")

data = data %>% dplyr::mutate(Sex = ifelse(Sex == 1, 0, 1),
                              Sex = as.factor(Sex), 
                              Origen = as.factor(Origen), 
                              Estadio_tanner = as.factor(Estadio_tanner)) %>% 
  dplyr::mutate(Talla_m = Talla_cm/100, .keep = "unused")



## Rename the names of the CpGs -----
predictors = colnames(data)[2:194]
names = data.frame(Name = predictors) %>% 
  dplyr::left_join(ann850k[, c("Name","UCSC_RefGene_Name")], by = "Name")
names$UCSC_RefGene_Name = stringr::str_replace(names$UCSC_RefGene_Name,";.*","")
names$def = base::paste(names$UCSC_RefGene_Name, names$Name, sep="_")
predictors = names$def
colnames(data)[2:194] = predictors
rm(names)

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

## Import and merge white cells populations and DXA -----
estimacion_celulas_blancas = read.csv2("./data/raw/2_validation_population_iberomics/ESTIMACION_CELULAS_BLANCAS_IBEROMICS.csv") %>% 
  dplyr::rename(Code = X)
dxa = iberomics2 %>% dplyr::select(Code, FMI_dexa) %>% 
  dplyr::rename(FMI = FMI_dexa)
data = data %>% merge(dxa, by = "Code") %>% 
  merge(estimacion_celulas_blancas, by = "Code")

## Export final datasets -----
saveRDS(data, "./data/processed/EWAS_iberomics.RDS")

