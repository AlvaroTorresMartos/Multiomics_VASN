
# Load packages -------
library(dplyr)
library(ggplot2)
library(forcats)

# Import and merge datasets 
pubmep_t1 = readRDS("./data/processed/vasn_serum_levels_pubmep_t1.RDS")
pubmep_t2 = readRDS("./data/processed/vasn_serum_levels_pubmep_t2.RDS")

pubmep_t1 = pubmep_t1 %>% 
  dplyr::select(-c(Code_new_T2)) %>% 
  dplyr::mutate(Cohort = 1)

pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(-c(Code_new_T2)) %>% 
  dplyr::mutate(Cohort = 1)

iberomics = readRDS("./data/processed/vasn_serum_levels_iberomics.RDS")

iberomics_pre = iberomics %>% 
  dplyr::filter(Estadio_tanner == 1) %>% 
  dplyr::mutate(Cohort = 2)

iberomics_pub = iberomics %>% 
  dplyr::filter(Estadio_tanner != 1) %>% 
  dplyr::mutate(Cohort = 2)

rm(iberomics)

prepuber = rbind(pubmep_t1, iberomics_pre) %>% 
  dplyr::mutate(Cohort = as.factor(Cohort), 
                Pubertal_stage = 1)

puber = rbind(pubmep_t2, iberomics_pub) %>% 
  dplyr::mutate(Cohort = as.factor(Cohort), 
                Pubertal_stage = 2)


data = rbind(prepuber, puber)

rm(pubmep_t1, pubmep_t2, iberomics_pre, iberomics_pub, prepuber, puber)

data %>% 
  tidyr::drop_na(Sex, log_vasn_pg_mL) %>% 
  dplyr::mutate(Cohort = ifelse(Cohort == 1, "Discovery Population", "Validation Population"), 
                Cohort = forcats::fct_relevel(Cohort, "Discovery Population", "Validation Population"),
                Pubertal_stage = ifelse(Pubertal_stage == 1, "Prepubertal", "Pubertal"), 
                Sex = ifelse(Sex == 0, "Female", "Male")
                ) %>% 
  ggplot(aes(x = Sex, y = log_vasn_pg_mL, fill = Sex)) +
  geom_boxplot() + # Crear el boxplot sin outliers visibles outlier.shape = NA
  geom_jitter(width = 0.2, alpha = 0.5) + # Agregar puntos para ver distribución individual
  facet_grid(Pubertal_stage ~ Cohort) + # Crear facetas con Pubertal_stage y Cohort
  labs(
    x = "Sex",
    y = "log(vasorin (pg/mL))", # "Vasorin (pg/mL)"
    title = "Serum Vasorin Protein Levels across Sex, Study Population, and Pubertal Stage"
  ) +
  scale_fill_manual(values = c("#FBB4AE", "#B3CDE3")) + # Colores personalizados
  theme(
    legend.position = "none", # Posicionar la leyenda arriba
    # panel.grid.major = element_line(color = "gray80"), # Color de las líneas del fondo
    panel.grid.minor = element_blank(),  # Eliminar líneas menores
    axis.title.x = element_text(size = 12, face = "bold"), # Tamaño del título del eje x
    axis.title.y = element_text(size = 12, face = "bold"), # Tamaño del título del eje y
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 11),  # Etiquetas del eje x
    axis.text.y = element_text(size = 11)
    )
 





 
# Crear el scatterplot agrupado
# data %>%
#   tidyr::drop_na(log_vasn_pg_mL) %>%
#   dplyr::mutate(Cohort = ifelse(Cohort == 1, "PUBMEP", "IBEROMICS"),
#                 Cohort = forcats::fct_relevel(Cohort, "PUBMEP", "IBEROMICS"),
#                 Pubertal_stage = ifelse(Pubertal_stage == 1, "Prepubertal", "Pubertal")
#   )  %>%
#   grouped_ggscatterstats(
#   ## Datos y variables
#   x = Vasn_pg_mL,
#   y = HOMA_IR_zscore,
#   grouping.var = interaction(Pubertal_stage, Cohort), # Agrupación por Pubertal_stage y Cohort
#   type = "p", # Correlación de Pearson
#   label.var = NULL, # Si quieres etiquetar puntos, puedes usar una columna
#   label.expression = NULL, # Puedes agregar una expresión lógica si quieres etiquetar valores específicos
# 
#   ## Argumentos para personalización
#   ggtheme = ggplot2::theme_minimal(), # Tema minimalista
#   # annotation.args = list(
#   #   title = "Relationship between Serum Vasorin Levels and FMI"#,
#   #   # caption = "Data grouped by pubertal stage"
#   # ),
#   plotgrid.args = list(nrow = 2, ncol = 2), # Diseño de la cuadrícula
#   messages = FALSE # Desactiva mensajes extra de ggstatsplot
# )
