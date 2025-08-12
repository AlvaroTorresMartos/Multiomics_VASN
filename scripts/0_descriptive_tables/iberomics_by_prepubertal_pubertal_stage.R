
# Load packages -------
library(dplyr)
library(tidyr)
library(forcats)
library(stats)
library(gtsummary)



# Import data -----
iberomics = readRDS("./data/processed/EWAS_iberomics.RDS")

iberomics_prot = readRDS("./data/processed/vasn_serum_levels_iberomics.RDS") %>% 
  tidyr::drop_na(Vasn_pg_mL) %>% 
  dplyr::select(Code, Vasn_pg_mL)

bmi_ibe = read.csv("./data/processed/bmi_iberomics.csv", row.names = 1)

iberomics = iberomics %>% 
  dplyr::left_join(iberomics_prot, by = "Code") %>% 
  dplyr::left_join(bmi_ibe, by = "Code")

rm(iberomics_prot, bmi_ibe)

iberomics = iberomics %>% 
  dplyr::select(Sex:Estadio_tanner, BMI, BMI_zscore, WC_zscore, FMI, Height_zscore, 
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
                Recluiment_center = ifelse(Recluiment_center == 1, "Santiago de Compostela", "Zaragoza"), 
                Experimental_group = ifelse(Tanner_scale == 1, "Prepubertal", "Pubertal"), 
                Experimental_group = forcats::fct_relevel(Experimental_group, "Prepubertal", "Pubertal")
                )

# Load BMI z-score (Cole et al 2012)
cole_zscore_iotf = function(data, cole_bmi_f, cole_bmi_m){
  data$Obesity_status = rep(NA)
  for(i in 1:nrow(data)){
    if(is.na(data$Age[i]) | is.na(data$Sex[i])){
    }else{
      if(data$Age[i] >= 2 & data$Age[i] <= 18 & data$Sex[i] == "Male"){
        pos = which.min(abs(cole_bmi_m$Age - data$Age[i]))
        data$Obesity_status[i] = (((data$BMI[i]/cole_bmi_m$M[pos])**cole_bmi_m$L[pos]) - 1)/(cole_bmi_m$L[pos] * cole_bmi_m$S[pos])
        print("Male")
      }
      if(data$Age[i] >= 2 & data$Age[i] <= 18 & data$Sex[i] == "Female"){
        pos = which.min(abs(cole_bmi_f$Age - data$Age[i]))
        data$Obesity_status[i] = (((data$BMI[i]/cole_bmi_f$M[pos])**cole_bmi_f$L[pos]) - 1)/(cole_bmi_f$L[pos] * cole_bmi_f$S[pos])
        print("Female")
      }
    }
  }
  data = data %>% 
    dplyr::mutate(Obesity_status2 = case_when(Age < 5 & Obesity_status > 2 &  Obesity_status < 3 ~ "Overweight", 
                                              Age > 5 & Obesity_status > 1 &  Obesity_status < 2 ~ "Overweight", 
                                              Age < 5 & Obesity_status >= 3 ~ "Obesity", 
                                              Age > 5 & Obesity_status >= 2 ~ "Obesity",
                                              TRUE ~ "Normal-weight"), 
                  Obesity_status2 = as.factor(Obesity_status2), 
                  Obesity_status2 = forcats::fct_relevel(Obesity_status2, "Normal-weight", "Overweight", "Obesity")
    ) %>% 
    dplyr::select(-c(Obesity_status)) %>% dplyr::rename(Obesity_status = Obesity_status2)
  return(data)
}


cole_ref_m = readxl::read_xlsx("/home/usuario/Escritorio/Shiny/ObMetrics/obmetrics_app/www/data/bmi_zscores/2012_cole_bmi.xlsx", sheet = 1, col_types = "numeric")
cole_ref_f = readxl::read_xlsx("/home/usuario/Escritorio/Shiny/ObMetrics/obmetrics_app/www/data/bmi_zscores/2012_cole_bmi.xlsx", sheet = 2, col_types = "numeric")

iberomics = iberomics %>% 
  cole_zscore_iotf(cole_bmi_f = cole_ref_f, cole_bmi_m = cole_ref_m)

# Statistical testing ----

vars_continuas = iberomics %>% 
  dplyr::select(Age, BMI:Vasorin_pg_mL) %>% 
  colnames()


tests_normalidad = sapply(vars_continuas, function(v) {
  x = iberomics[[v]]
  x1 = iberomics[[v]][iberomics$Experimental_group == "Prepubertal"] %>% 
    na.omit()
  x2 = iberomics[[v]][iberomics$Experimental_group == "Pubertal"] %>% 
    na.omit()
  n1 = length(x1)
  n2 = length(x2)
  print(v)
  if(n1 < 50 & n2 < 50){
    print("Shapiro-Wilk test")
    p1 = shapiro.test(x1)$p.value
    p2 = shapiro.test(x2)$p.value
  } else{
    print("Anderson-Darling test") # Because Kolmogorov-Smirnov Tests gave us problems due to ties
    p1 = nortest::ad.test(x1)$p.value
    p2 = nortest::ad.test(x2)$p.value
  }
  print(p1)
  print(p2)
  # h0: normality
  # h0: no normality
  # TRUE no-parametric
  # FALSE parametric 
  any(c(p1, p2) < 0.05)
})


parametrics = vars_continuas[tests_normalidad == FALSE]
non_parametrics = vars_continuas[tests_normalidad == TRUE]



# Descriptive table  -----
iberomics %>% 
  tbl_summary(by = Experimental_group, 
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
  ) %>% 
  add_p(
    test = list(all_of(parametrics) ~ "t.test", all_of(non_parametrics) ~ "wilcox.test"), 
    test.args = all_continuous() ~ list(var.equal = FALSE),
    pvalue_fun = label_style_pvalue(digits = 2)
    )

# Information about the statistical tests
# help("add_p")
# help("add_p.tbl_summary")

# For Two Groups
# Continuous Data
# 
# - t-test: For normal data. Compares the means of two groups.
# - Wilcoxon rank-sum test: For non-normal or heterogeneous data. Compares distributions between two groups.
# - Paired t-test: For paired (before/after) data that follows a normal distribution.
# - Paired Wilcoxon test: For paired (before/after) data that does not follow a normal distribution.
# 
# Categorical Data
# 
# - Chi-square test of independence: Compares proportions between two groups.
# - Chi-square test without correction: Same as above but without correction, suitable for large samples.
# - Fisher's exact test: For small contingency tables (frequencies < 5).
# - McNemar's test: For paired or dependent data (e.g., before and after a treatment).
# 
# For More Than Two Groups
# Continuous Data
# 
# - One-way ANOVA (aov): For normal data with homogeneous variances. Compares means across multiple groups.
# - One-way ANOVA (oneway.test): For normal data with heterogeneous variances.
# - Kruskal-Wallis test: For non-normal or heterogeneous data. Compares distributions across multiple groups.
# - Mood's test: Evaluates differences in dispersion across multiple groups.
# 
# Categorical Data
# 
# - Chi-square test of independence: Compares proportions across multiple groups.
# - Fisher's exact test: For small contingency tables.
# - Proportions test (prop.test): Compares proportions across multiple groups.
# 
# Data Adjusted for Covariates
# 
# - ANCOVA: For continuous data, adjusts comparisons for covariates.
# - Random intercept logistic regression (lme4): For categorical data with hierarchical or nested structures.