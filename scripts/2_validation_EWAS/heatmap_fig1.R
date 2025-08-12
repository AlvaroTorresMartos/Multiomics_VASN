
# https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html
# https://divingintogeneticsandgenomics.com/post/how-to-make-a-triangle-correlation-heatmap-with-p-values-labeled/

# Load packages -------
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(Hmisc)
library(ComplexHeatmap)
library(grid)
library(circlize)

# By cohort -----
assoc_pub = read.csv2("./results/Validation_EWAS/pubmep/2024_12_04_pubmep_associations_between_DNA_methylation_193_CpGs_and_outcomes.csv") 
assoc_ibe = read.csv2("./results/Validation_EWAS/iberomics/2024_12_04_iberomics_associations_between_DNA_methylation_193_CpGs_and_outcomes.csv") 

# Estos son los que pasan p.value y que vamos a mostrar (quitamos QUICKI, LDLc, Insulin porque no forman parte del MetS)
outcomes = c("BMI_zscore", "WC_zscore", "FMI",  
             "HOMA_IR_zscore", 
             # "Glucose_zscore",  
             "BP_zscore_NHBPEP", "BP_zscore_Stavnsbo", "DBP_zscore_NHBPEP", 
             "DBP_zscore_Stavnsbo", "SBP_zscore_NHBPEP", "SBP_zscore_Stavnsbo", 
             # "TAG_zscore", 
             "HDLc_zscore")

assoc_pub2 = assoc_pub %>% 
  dplyr::filter(outcome %in% outcomes) %>% 
  dplyr::filter(p.value < 0.05)

outcomes1 = names(table(assoc_pub2$outcome))
predictor1 = names(table(assoc_pub2$predictor))

assoc_ibe2 = assoc_ibe %>% 
  dplyr::filter(outcome %in% outcomes) %>% 
  dplyr::filter(p.value < 0.05)

outcomes2 = names(table(assoc_ibe2$outcome))
predictor2 = names(table(assoc_ibe2$predictor))

outcomes = outcomes1[outcomes1 %in% outcomes2]
predictors = predictor1[predictor1 %in% predictor2]

rm(outcomes1, outcomes2, predictor1, predictor2)

assoc_pub3 = assoc_pub2[assoc_pub2$outcome %in% outcomes & assoc_pub2$predictor %in% predictors, ] 
assoc_ibe3 = assoc_ibe2[assoc_ibe2$outcome %in% outcomes & assoc_ibe2$predictor %in% predictors, ]

assoc_pub3 = assoc_pub3 %>%
  dplyr::mutate(outcome = ifelse(stringr::str_detect(outcome, "BP_zscore"), "BP_zscore_(any)", outcome)) %>%
  dplyr::group_by(outcome, predictor) %>%
  dplyr::arrange(fdr, desc(abs(std_coefficient))) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

assoc_ibe3 = assoc_ibe3 %>%
  dplyr::mutate(outcome = ifelse(stringr::str_detect(outcome, "BP_zscore"), "BP_zscore_(any)", outcome)) %>%
  dplyr::group_by(outcome, predictor) %>%
  dplyr::arrange(fdr, desc(abs(std_coefficient))) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

rm(assoc_ibe, assoc_ibe2, assoc_pub, assoc_pub2)

out_pred = dplyr::inner_join(assoc_pub3, assoc_ibe3, by = c("outcome", "predictor")) %>% 
  dplyr::select(outcome, predictor)

assoc_pub3 = left_join(out_pred, assoc_pub3, by = c("outcome", "predictor")) %>% 
  dplyr::select(outcome, predictor, std_coefficient, fdr, p.value) %>% 
  dplyr::mutate(outcome = ifelse(stringr::str_detect(outcome, "\\(PUBMEP\\)"), 
                                 outcome, paste0(outcome, "_(PUBMEP)")))

assoc_ibe3 = left_join(out_pred, assoc_ibe3, by = c("outcome", "predictor")) %>% 
  dplyr::select(outcome, predictor, std_coefficient, fdr, p.value) %>% 
  dplyr::mutate(outcome = ifelse(stringr::str_detect(outcome, "\\(IBEROMICS\\)"), 
                                 outcome, paste0(outcome, "_(IBEROMICS)")))

assoc3 = rbind(assoc_pub3, assoc_ibe3)

rm(assoc_pub3, assoc_ibe3, out_pred, predictors, outcomes)

wide_beta = assoc3 %>% dplyr::select(outcome, predictor, std_coefficient)
wide_pvalue = assoc3 %>% dplyr::select(outcome, predictor, p.value)
wide_fdr = assoc3 %>% dplyr::select(outcome, predictor, fdr)


# wide beta 
wide_beta = tidyr::pivot_wider(wide_beta, names_from = outcome,
                               values_from = std_coefficient) %>%
  tibble::column_to_rownames(var = "predictor")

outcomes = c("BMI_zscore_(PUBMEP)", "BMI_zscore_(IBEROMICS)",
             "WC_zscore_(PUBMEP)" , "WC_zscore_(IBEROMICS)", 
             "FMI_(PUBMEP)", "FMI_(IBEROMICS)", 
             "HOMA_IR_zscore_(PUBMEP)",  "HOMA_IR_zscore_(IBEROMICS)", 
             # "Glucose_zscore_(PUBMEP)", "Glucose_zscore_(IBEROMICS)", 
             "BP_zscore_(any)_(PUBMEP)", "BP_zscore_(any)_(IBEROMICS)", 
             # "TAG_zscore_(PUBMEP)", "TAG_zscore_(IBEROMICS)", 
             "HDLc_zscore_(PUBMEP)", "HDLc_zscore_(IBEROMICS)"
)

wide_beta = wide_beta[, c(outcomes)]

# wide_beta = wide_beta %>% 
#   dplyr::mutate(across(everything(), ~replace_na(.x, 0))) 

wide_beta = wide_beta %>% 
  tibble::rownames_to_column(var = "var") %>% 
  dplyr::mutate(pub = ifelse(grepl("prepubertal", var), 0, 1)) %>% 
  dplyr::arrange(pub, var) %>% 
  tibble::column_to_rownames(var = "var") %>% 
  dplyr::select(-c(pub))

# wide_fdr
wide_fdr = tidyr::pivot_wider(wide_fdr, names_from = outcome,
                              values_from = fdr) %>%
  tibble::column_to_rownames(var = "predictor")

wide_fdr = wide_fdr[, c(outcomes)]

# wide_fdr = wide_fdr %>%
#   dplyr::mutate(across(everything(), ~replace_na(.x, 1)))

wide_fdr = wide_fdr %>% 
  tibble::rownames_to_column(var = "var") %>% 
  dplyr::mutate(pub = ifelse(grepl("prepubertal", var), 0, 1)) %>% 
  dplyr::arrange(pub, var) %>% 
  tibble::column_to_rownames(var = "var") %>% 
  dplyr::select(-c(pub))

wide_fdr = wide_fdr %>%
  dplyr::mutate(across(everything(), ~replace_na(.x, 1)))

# wide_pvalue
wide_pvalue = tidyr::pivot_wider(wide_pvalue, names_from = outcome,
                                 values_from = p.value) %>%
  tibble::column_to_rownames(var = "predictor")

wide_pvalue = wide_pvalue[, c(outcomes)]

# wide_pvalue = wide_pvalue %>%
#   dplyr::mutate(across(everything(), ~replace_na(.x, 1)))

wide_pvalue = wide_pvalue %>% 
  tibble::rownames_to_column(var = "var") %>% 
  dplyr::mutate(pub = ifelse(grepl("prepubertal", var), 0, 1)) %>% 
  dplyr::arrange(pub, var) %>% 
  tibble::column_to_rownames(var = "var") %>% 
  dplyr::select(-c(pub))

wide_pvalue = wide_pvalue %>%
  dplyr::mutate(across(everything(), ~replace_na(.x, 1)))

# quitamos un gen con comportamiento extraño 

wide_beta = wide_beta[rownames(wide_beta) %nin% c("FGD4_cg03418231_pubertal", "TNXB_cg16416464_pubertal"), ]
wide_fdr = wide_fdr[rownames(wide_fdr) %nin% c("FGD4_cg03418231_pubertal", "TNXB_cg16416464_pubertal"), ]
wide_pvalue = wide_pvalue[rownames(wide_pvalue) %nin% c("FGD4_cg03418231_pubertal", "TNXB_cg16416464_pubertal"), ]


cell_fun_full <- function(j, i, x, y, w, h, fill){
  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white", lwd = 1.5))
  
  if (wide_fdr[i, j]  < 0.001) {
    grid.text(paste0(sprintf("%.2f", wide_beta[i, j]), "***"), 
              x, y, gp = gpar(fontsize = 13))
  } else if (wide_fdr[i, j]  <= 0.01) {
    grid.text(paste0(sprintf("%.2f", wide_beta[i, j]), "**"), 
              x, y, gp = gpar(fontsize = 13))
  } else if (wide_fdr[i, j]  <= 0.05) {
    grid.text(paste0(sprintf("%.2f", wide_beta[i, j]), "*"), 
              x, y, gp = gpar(fontsize = 13))
  } else if(wide_pvalue[i, j] <= 0.05) {
    grid.text(sprintf("%.2f", wide_beta[i, j]),
              x, y, gp = gpar(fontsize = 11))
  }
}

final_names = colnames(wide_beta)

final_names = final_names %>% 
  gsub(pattern = "PUBMEP", replacement = "Discovery_Population")

final_names = final_names %>%
  gsub(pattern = "IBEROMICS", replacement = "Validation_Population")

colnames(wide_beta) = final_names

hp_full = ComplexHeatmap::Heatmap(as.matrix(wide_beta), 
                                  name = "Standarized Coefficient",
                                  col = circlize::colorRamp2(
                                    c(-0.7, 0, 0.7), 
                                    c("#4575b4", "#FAF3E0", "#f46d43")),
                                  na_col = "grey",
                                  cluster_columns = FALSE, 
                                  cluster_rows = FALSE, 
                                  row_names_side = "left",
                                  column_names_side = "bottom", 
                                  column_names_rot = 45,
                                  row_names_gp = gpar(fontsize = 11.5),  
                                  column_names_gp = gpar(fontsize = 12),
                                  cell_fun = cell_fun_full,
                                  rect_gp = gpar(col = "white", lwd = 0.5), 
                                  heatmap_legend_param = list(
                                    legend_direction = "horizontal", 
                                    legend_width = unit(3, "cm"))
)



# Dibujar el heatmap
draw(hp_full,  ht_gap = unit(1, "cm"), heatmap_legend_side = "top")



# cpgs_validadas_fig1 = wide_beta %>% row.names()
# write.csv2(cpgs_validadas_fig1, "./results/Validation_EWAS/CpGs_validadas.csv", 
#            row.names = FALSE)

