

library(janitor)
library(dplyr)

vasn = readxl::read_xls("./data/raw/3_additional_cohort/vasn_vat.xls", sheet = 1) %>% 
  janitor::clean_names()

# Cole, T.J.; Bellizi, M.C.; FLegal, K.M.; Dietz, W.H. Establishing a standard definition for child overweight and obesity worldwide: International survey. BMJ 2000, 320, 1240–1243.
metadata = readxl::read_xls("./data/raw/3_additional_cohort/metadata_ines_rtpcr.xls", sheet = 1) %>% 
  janitor::clean_names() %>% 
  dplyr::rename(id = no_rna) %>% 
  dplyr::select(id, sexo, edad, imc, obesidad, adiponectin_ug_ml)


data = merge(vasn, metadata, by = "id") %>% 
  dplyr::arrange(obesidad, vat_expression) %>% 
  dplyr::mutate(sexo = case_when(sexo == "MUJER" ~ 1, TRUE ~ 0)) %>% 
  dplyr::rename(BMI = imc, sex = sexo, decimal_age = edad)
  
ref_val = load("/home/usuario/Escritorio/Shiny/ObMetrics/obmetrics_app/www/inputs/reference_tables.RData")


source("/home/usuario/Escritorio/Shiny/ObMetrics/obmetrics_app/outputs_functions/1_basic_functions.R")
source("/home/usuario/Escritorio/Shiny/ObMetrics/obmetrics_app/outputs_functions/3_mets_components_zscores.R")

data2 = data %>% 
     cole_zscore_iotf(cole_bmi_f = cole_ref_f, cole_bmi_m = cole_ref_m) %>% 
       dplyr::rename(bmi_zscore = Excess_adiposity)
     

library(gtsummary)
data2 %>% 
  dplyr::filter(id != 18) %>% 
  dplyr::select(sex, decimal_age, bmi_zscore, obesidad, adiponectin_ug_ml, vat_expression) %>% 
  dplyr::mutate(sex = as.factor(sex), 
                decimal_age = as.numeric(decimal_age), 
                adiponectin_ug_ml = as.numeric(adiponectin_ug_ml)) %>% 
  tbl_summary(
    by = obesidad, 
    type = list(
      decimal_age ~ "continuous2",
      all_continuous() ~ "continuous2"), 
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
  ) %>% 
  add_p(pvalue_fun = label_style_pvalue(digits = 2))




