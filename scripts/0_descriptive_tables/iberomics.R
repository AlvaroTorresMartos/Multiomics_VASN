
# Load packages -------
library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)



# Import data -----
iberomics = readRDS("./data/processed/EWAS_iberomics.RDS")

iberomics_prot = readRDS("./data/processed/vasn_serum_levels_iberomics.RDS") %>% 
  tidyr::drop_na(Vasn_pg_mL) %>% 
  dplyr::select(Code, Vasn_pg_mL)

iberomics = iberomics %>% 
  dplyr::left_join(iberomics_prot, by = "Code")

rm(iberomics_prot)

iberomics = iberomics %>% 
  dplyr::select(Sex:BMI_zscore, WC_zscore, FMI, Height_zscore, 
                Glucose_zscore, Insulin_zscore, HOMA_IR_zscore, 
                DBP_zscore_NHBPEP, SBP_zscore_NHBPEP, BP_zscore_NHBPEP, 
                TAG_zscore, HDLc_zscore, LDLc_zscore, Vasn_pg_mL)


iberomics = iberomics %>% 
  dplyr::rename(Recluiment_center = Origen, 
                Tanner_scale = Estadio_tanner, 
                DBP_zscore = DBP_zscore_NHBPEP, 
                SBP_zscore = SBP_zscore_NHBPEP, 
                BP_zscore = BP_zscore_NHBPEP, 
                Vasorin_pg_mL = Vasn_pg_mL
                ) %>% 
  dplyr::mutate(Sex = ifelse(Sex == 0, "Male", "Female"), 
                Sex = forcats::fct_relevel(Sex, "Male", "Female"),
                Recluiment_center = ifelse(Recluiment_center == 1, "Santiago de Compostela", "Zaragoza"))

# Descriptive table  -----
iberomics %>% 
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
