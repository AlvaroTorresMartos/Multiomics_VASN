
# https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html
# https://divingintogeneticsandgenomics.com/post/how-to-make-a-triangle-correlation-heatmap-with-p-values-labeled/

# Load packages -------
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(grid)
library(circlize)

# BOTH -----
assoc = read.csv2("./results/Validation_EWAS/pubmep/2024_12_04_pubmep_associations_between_DNA_methylation_193_CpGs_and_outcomes.csv") 
# Estos son los que pasan FDR y que vamos a mostrar
# (quitamos QUICKI e Insulin porque no forman parte del MetS)

outcomes = c("BMI_zscore", "WC_zscore", "FMI",  "HOMA_IR_zscore",
             "BP_zscore_NHBPEP", "BP_zscore_Stavnsbo", 
             "DBP_zscore_NHBPEP", "DBP_zscore_Stavnsbo", 
             "SBP_zscore_NHBPEP", "SBP_zscore_Stavnsbo", 
             "HDLc_zscore"
             )



assoc2 = assoc %>% 
  dplyr::filter(outcome %in% outcomes) %>% 
  dplyr::filter(fdr < 0.05)

outcomes = names(table(assoc2$outcome))
predictors = names(table(assoc2$predictor))

assoc3 = assoc[assoc$outcome %in% outcomes & assoc$predictor %in% predictors, ] 

assoc3 = assoc3 %>%
  dplyr::mutate(outcome = ifelse(stringr::str_detect(outcome, "BP_zscore"), "BP_zscore_(any)", outcome)) %>%
  dplyr::group_by(outcome, predictor) %>%
  dplyr::arrange(fdr, desc(abs(std_coefficient))) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

wide_beta = assoc3 %>% dplyr::select(outcome, predictor, std_coefficient)
wide_pvalue = assoc3 %>% dplyr::select(outcome, predictor, p.value)
wide_fdr = assoc3 %>% dplyr::select(outcome, predictor, fdr)

# wide beta 
wide_beta = tidyr::pivot_wider(wide_beta, names_from = outcome,
                               values_from = std_coefficient) %>%
  tibble::column_to_rownames(var = "predictor")

outcomes = c("BMI_zscore", "WC_zscore", "FMI",  "HOMA_IR_zscore",
             "BP_zscore_(any)", "HDLc_zscore"
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


for(i in 1:nrow(wide_beta)){
  for(j in 1:ncol(wide_beta)){
    if(wide_pvalue[i, j] >= 0.05){
      wide_beta[i, j] = NA
    }
  }
}


cell_fun_full <- function(j, i, x, y, w, h, fill){
  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white", lwd = 1.5))
  
  if (wide_fdr[i, j]  < 0.001) {
    grid.text(paste0(sprintf("%.2f", wide_beta[i, j]), "***"), 
              x, y, gp = gpar(fontsize = 11))
  } else if (wide_fdr[i, j]  <= 0.01) {
    grid.text(paste0(sprintf("%.2f", wide_beta[i, j]), "**"), 
              x, y, gp = gpar(fontsize = 11))
  } else if (wide_fdr[i, j]  <= 0.05) {
    grid.text(paste0(sprintf("%.2f", wide_beta[i, j]), "*"), 
              x, y, gp = gpar(fontsize = 11))
  } else if(wide_pvalue[i, j] < 0.05) {
    grid.text(sprintf("%.2f", wide_beta[i, j]),
              x, y, gp = gpar(fontsize = 9))
  }
}

hp_full = ComplexHeatmap::Heatmap(as.matrix(wide_beta), 
                                  name = "Standarized Coefficient",
                                  col = circlize::colorRamp2(
                                    c(-0.8, 0, 0.8), 
                                    c("#4575b4", "#FAF3E0", "#f46d43")),
                                  na_col = "grey",
                                  cluster_columns = FALSE, 
                                  cluster_rows = FALSE, 
                                  row_names_side = "left",
                                  column_names_side = "bottom", 
                                  column_names_rot = 45,
                                  row_names_gp = gpar(fontsize = 9.5),  
                                  column_names_gp = gpar(fontsize = 12),
                                  cell_fun = cell_fun_full,
                                  rect_gp = gpar(col = "white", lwd = 0.5), 
                                  heatmap_legend_param = list(
                                    legend_direction = "horizontal", 
                                    legend_width = unit(3, "cm"))
)
                                   


# Dibujar el heatmap
draw(hp_full,  ht_gap = unit(1, "cm"), heatmap_legend_side = "top")

# CpG_fig1b = rownames(wide_beta)
