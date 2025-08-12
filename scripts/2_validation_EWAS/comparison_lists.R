

# g5vsg3_between = read.csv2("./data/processed/comparisons_between_groups_lists/between_OBIRNOCHANGE VS OBIRtononIR_RAW0001andFDR05_AGE.csv")
# 
# g3vsg4_between = read.csv2("./data/processed/comparisons_between_groups_lists/between_OBOWIRtononIR VS OBOWnonIRtoIR_RAW0001_AGE.csv")
# 
# g2vsg4_between = read.csv2("./data/processed/comparisons_between_groups_lists/between_OBOWnonIRnochange VS OBOWnonIRtoIR_RAW0001_AGE.csv")
# 
# g3_within = read.csv2("./data/processed/comparisons_between_groups_lists/DMPs_significants_adjusted_OB_OW_IR_to_non_IR_WITHIN_raw0001_age.csv")
# 
# g4_within = read.csv2("./data/processed/comparisons_between_groups_lists/DMPs_significants_adjusted_OB_OW_non_IR_to_IR_WITHIN_raw0001_age.csv")

library(dplyr)
library(janitor)
library(dplyr)




lista_def = readxl::read_xlsx("./data/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.xlsx")

lista_def2 = lista_def %>% 
  dplyr::mutate(gene = paste0(Name, "_", UCSC_RefGene_Name), 
                gene = gsub(pattern = "sep-01", replacement = "SEPT1", x = gene),
                gene = gsub(pattern = "_NA", replacement = "", x = gene)
                ) %>% 
  janitor::clean_names() %>% 
  dplyr::select(chr, gene, pos, p_value, approach_group_comparison)

g5vsg3_between = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G5 vs. G3") %>% 
  dplyr::select(chr, gene, pos, p_value) 

g3vsg4_between = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G3 vs. G4") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)

g2vsg4_between = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G2 vs. G4") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)

g3_within = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G3 (within)") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)

g4_within = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G4 (within)") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)

nw_owob_cross = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Pubertal non-IR NW vs. IR Obese & Overweight") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)

owob_owb_cross = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Pubertal non-IR Obese & Overweight vs. IR Obese & Overweight") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)

noir_ir_cross = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Pubertal IR vs. non-IR") %>% 
  dplyr::select(chr, gene, pos, p_value) %>% 
  rbind(df)












