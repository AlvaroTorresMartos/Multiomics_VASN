
# Load packages -------
library(dplyr)
library(ggplot2)
library(forcats)

# Import results -----

both_pre = read.csv("./results/Serum_vasorin_levels/both/2024_12_11_prepubertal_both_association_between_serum_vasorin_levels_and_outcomes.csv")
both_pub = read.csv("./results/Serum_vasorin_levels/both/2024_12_11_pubertal_both_association_between_serum_vasorin_levels_and_outcomes.csv")

pubmed_pre = read.csv("./results/Serum_vasorin_levels/pubmep/2024_12_11_prepubertal_pubmep_association_between_serum_vasorin_levels_and_outcomes.csv")
pubmep_pub = read.csv("./results/Serum_vasorin_levels/pubmep/2024_12_11_pubertal_pubmep_association_between_serum_vasorin_levels_and_outcomes.csv")

iberomics_pre = read.csv("./results/Serum_vasorin_levels/iberomics/2024_12_11_prepubertal_iberomics_association_between_serum_vasorin_levels_and_outcomes.csv")
iberomics_pub = read.csv("./results/Serum_vasorin_levels/iberomics/2024_12_11_pubertal_iberomics_association_between_serum_vasorin_levels_and_outcomes.csv")



results = rbind(pubmed_pre, pubmep_pub, iberomics_pre, iberomics_pub, 
               both_pre, both_pub)

rm(pubmed_pre, pubmep_pub, iberomics_pre, iberomics_pub, both_pre, both_pub)

results = results %>% 
  dplyr::mutate(sex = "both")
  
male_pre = read.csv("./results/Serum_vasorin_levels/both/by_sex/2024_12_11_prepubertal_males_both_association_between_serum_vasorin_levels_and_outcomes.csv")
fema_pre = read.csv("./results/Serum_vasorin_levels/both/by_sex/2024_12_11_prepubertal_females_both_association_between_serum_vasorin_levels_and_outcomes.csv")

male_pub = read.csv("./results/Serum_vasorin_levels/both/by_sex/2024_12_11_pubertal_males_both_association_between_serum_vasorin_levels_and_outcomes.csv")
fema_pub = read.csv("./results/Serum_vasorin_levels/both/by_sex/2024_12_11_pubertal_females_both_association_between_serum_vasorin_levels_and_outcomes.csv")

results = rbind(results, male_pre, fema_pre, male_pub, fema_pub)

rm(male_pre, fema_pre, male_pub, fema_pub)

results = results %>% 
  dplyr::mutate(CI_low_std = (CI_low*std_coefficient)/beta_coefficient, 
                CI_high_std = (CI_high*std_coefficient)/beta_coefficient) %>% 
  dplyr::filter(outcome %in% c("BMI_zscore","WC_zscore", "FMI",  
                               "HOMA_IR_zscore", "BP_zscore_NHBPEP",
                               "HDLc_zscore")) %>% 
  dplyr::mutate(outcome = ifelse(outcome %in% "BP_zscore_NHBPEP", 
                                 "BP_zscore", outcome))


results %>%
  dplyr::mutate(
    cohort = ifelse(cohort == "pubmep", "Discovery Population", 
                    ifelse(cohort == "iberomics", "Validation Population", "Pool")
                    ), 
    cohort = forcats::fct_relevel(cohort, "Discovery Population", "Validation Population",  "Pool"),
    time = ifelse(time == "prepubertal", "Prepubertal", "Pubertal"), 
    outcome = factor(outcome), 
    outcome = forcats::fct_relevel(outcome, "HDLc_zscore", "BP_zscore",
                                   "HOMA_IR_zscore", "FMI",
                                   "WC_zscore", "BMI_zscore"), 
    significant = p.value < 0.05,
    outcome_sex = ifelse(sex %in% c("female", "male"),
                         paste0(outcome, "_", sex),
                         as.character(outcome)),
    outcome_sex = factor(outcome_sex),
    outcome_sex = forcats::fct_relevel(outcome_sex,
                                       c("HDLc_zscore_female", "HDLc_zscore_male", "HDLc_zscore",
                                         "BP_zscore_female", "BP_zscore_male", "BP_zscore",
                                         "HOMA_IR_zscore_female", "HOMA_IR_zscore_male", "HOMA_IR_zscore",
                                         "FMI_female", "FMI_male", "FMI",
                                         "WC_zscore_female", "WC_zscore_male", "WC_zscore",
                                         "BMI_zscore_female", "BMI_zscore_male", "BMI_zscore"))
  ) %>% 
  ggplot(aes(x = std_coefficient, y = outcome_sex, shape = sex)) +
  facet_grid(cohort ~ time, scales = "free", space = "free_y") + 
  geom_point(aes(color = significant), size = 4) +
  geom_errorbar(aes(xmin = CI_low_std, xmax = CI_high_std, color = significant), 
                width=0, size = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "Standardized Coefficient",
     y = "Outcomes") +
  scale_shape_manual(
    name = "Sex",
    values = c("both" = 16, "male" = 15, "female" = 17),
    labels = c("both" = "Both", "male" = "Male", "female" = "Female")
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "#FF6961"), 
    labels = c("FALSE" = "p-value > 0.05", "TRUE" = "p-value < 0.05"),
    name = "Statistically Significant"
                     ) + 
  scale_y_discrete(
    labels = function(y) {
      ifelse(grepl("_male|_female", y), "", gsub("_both", "", y))
    }
  ) +
  scale_x_continuous(breaks=seq(-0.5, 0.3, by=0.2)) + 
  theme_grey(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9), 
    legend.title = element_text(size = 9, face = "bold"),
    axis.ticks.y = element_blank(), 
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black") 
  )

# geom_text(
#   aes(label = ifelse(p.value < 0.05, "*", "")), 
#   nudge_x = -0.2, # Ajusta el desplazamiento horizontal
#   size = 10,       # Ajusta el tamaño del asterisco
#   show.legend = FALSE
# ) +
  


  
  
