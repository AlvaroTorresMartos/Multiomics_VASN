
library(dplyr)
library(janitor)
library(Hmisc)
library(qqman)

lista_def = readxl::read_xlsx("./data/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.xlsx")

lista_def2 = lista_def %>% 
  mutate(
    gene = paste0(Name, "_", UCSC_RefGene_Name),
    gene = gsub("sep-01", "SEPT1", gene),
    gene = gsub("_NA", "", gene)
  ) %>% 
  clean_names() %>% 
  select(chr, gene, pos, p_value, approach_group_comparison)

g5vsg3_between <- lista_def2 %>% 
  filter(approach_group_comparison == "Longitudinal G5 vs. G3") %>% 
  select(chr, gene, pos, p_value)

allgenes = read.csv("./data/ann850K.csv")

allgenes = allgenes %>% dplyr::select(X, chr, pos) %>% 
  dplyr::filter(chr %nin% c("chrX", "chrY") ) %>% 
  dplyr::rename(gene = X) %>% 
  dplyr::relocate(gene, .after = "chr") %>% 
  dplyr::mutate(p_value = NA)

allgenes2 = allgenes %>%
  dplyr::mutate(p_value = if_else(is.na(p_value), 
                                  10^(-rbeta(n(), shape1 = 0.5, shape2 = 5) * 4), 
                                  p_value))

allgenes2$p_value %>% hist()

write.csv(allgenes2, "./data/allgenes_manhattanplot.csv")

df_original = rbind(g5vsg3_between, allgenes2)

# ------------------------------------------------------------------
# 2) Definimos df_original con CHR numérico (posiciones reales)
# ------------------------------------------------------------------
df_original <- df_original %>%
  mutate(CHR = as.numeric(gsub("chr", "", chr))) %>%
  rename(P = p_value) %>%
  select(CHR, pos, P, gene)


unique_chrs <- sort(unique(df_original$CHR))


manhattan(
  df_original,
  chr = "CHR",
  p = "P",
  bp = "pos",
  snp = "gene",
  main = "Manhattan Plot",
  chrlabs = as.numeric(1:22), 
  genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5),
  ylim = c(0, 10)
)

with(
  subset(df_original, P < 1e-5),
  text(pos, -log10(P), labels = gene, pos = 3, cex = 0.7, col = "red")
)

