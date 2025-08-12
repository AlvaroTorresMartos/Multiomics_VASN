
# Load packages -------
library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)



# Import data -----
pubmep_t1 = readRDS("./data/processed/EWAS_pubmep_t1.RDS")

pubmep_t1prot = readRDS("./data/processed/vasn_serum_levels_pubmep_t1.RDS") %>% 
  tidyr::drop_na(Vasn_pg_mL) %>% 
  dplyr::select(Code, Vasn_pg_mL)

pubmep_t1 = pubmep_t1 %>% 
  dplyr::left_join(pubmep_t1prot, by = "Code")

rm(pubmep_t1prot)

pubmep_t1 = pubmep_t1 %>% 
  dplyr::select(Sex:BMI_zscore, WC_zscore, FMI, Height_zscore, 
                Glucose_zscore, Insulin_zscore, HOMA_IR_zscore, 
                DBP_zscore_NHBPEP, SBP_zscore_NHBPEP, BP_zscore_NHBPEP, 
                TAG_zscore, HDLc_zscore, LDLc_zscore, Vasn_pg_mL)


pubmep_t1 = pubmep_t1 %>% 
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
                Tanner_scale = 1)

# Descriptive table  -----
pubmep_t1 %>% 
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
