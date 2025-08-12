
# Load packages -------
library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)



# Import data -----
info = read_xlsx("./data/raw/1_discovery_population_pubmep/BD_METABOLOMICS.xlsx") 
info = info %>% 
  dplyr::select(Code, RNAseqT2) %>% 
  dplyr::mutate(RNAseqT2 = as.numeric(RNAseqT2)) %>% 
  dplyr::filter(RNAseqT2 == 1)

pubmep_t2 = readRDS("./data/processed/EWAS_pubmep_t2.RDS")

pubmep_t2prot = readRDS("./data/processed/vasn_serum_levels_pubmep_t2.RDS") %>% 
  tidyr::drop_na(Vasn_pg_mL) %>% 
  dplyr::select(Code, Vasn_pg_mL)

pubmep_t2 = pubmep_t2 %>% 
  dplyr::left_join(pubmep_t2prot, by = "Code") %>% 
  dplyr::filter(Code %in% info$Code)

# Use the script pubmep t1 t2  instead of this
pubmep_t2 = pubmep_t2 %>% 
  # dplyr::left_join(pubmep_t2prot, by = "Code") %>% 
  dplyr::filter(Code %in% info$Code) %>% 
  dplyr::select(-c(Code))

rm(pubmep_t2prot, info)

pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(Sex:BMI_zscore, WC_zscore, FMI, Height_zscore, 
                Glucose_zscore, Insulin_zscore, HOMA_IR_zscore, 
                DBP_zscore_NHBPEP, SBP_zscore_NHBPEP, BP_zscore_NHBPEP, 
                TAG_zscore, HDLc_zscore, LDLc_zscore, Vasn_pg_mL)


pubmep_t2 = pubmep_t2 %>% 
  dplyr::rename(Recluiment_center = Origen, 
                Tanner_scale = Estadio_tanner, 
                DBP_zscore = DBP_zscore_NHBPEP, 
                SBP_zscore = SBP_zscore_NHBPEP, 
                BP_zscore = BP_zscore_NHBPEP, 
                Vasorin_pg_mL = Vasn_pg_mL
                ) %>% 
  dplyr::mutate(Sex = ifelse(Sex == 0, "Male", "Female"), 
                Sex = forcats::fct_relevel(Sex, "Male", "Female"),
                Recluiment_center = ifelse(Recluiment_center == 0, "Córdoba",
                                    ifelse(Recluiment_center == 1, "Santiago de Compostela", "Zaragoza")), 
                Tanner_scale = ifelse(Tanner_scale == 1, 2, Tanner_scale))

# Descriptive table  -----
pubmep_t2 %>% 
  tbl_summary(
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})"),
    digits = list(
      all_continuous() ~ 2,
      all_categorical() ~ c(0,2)
    ), 
    missing_text = "Missing values"
    ) %>%
  modify_table_styling(
    columns = "label",
    rows = row_type == "label",   
    text_format = "bold"
  )
