# Code
# library used
library(readr)
library(readxl)
library(dplyr)

##   : 'dplyr'

## The following objects are masked from 'package:stats':
##
## filter, lag

## The following objects are masked from 'package:base':
##
## intersect, setdiff, setequal, union

library(tidyverse)
## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
## v ggplot2 3.3.5 v purrr 0.3.4
## v tibble 3.1.4 v stringr 1.4.0
## v tidyr 1.1.4 v forcats 0.5.1
## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag() masks stats::lag()

library(ggplot2)
library(ggridges)
library(GGally)

## Registered S3 method overwritten by 'GGally':
## method from
## +.gg ggplot2

library(cowplot)

# basic data
## loading data
dir <- "C:/Users/USER/Desktop/sue/portfolio_tables"
list.files(path = dir)
## [1] "~$tab_1.xlsx" "~$tab_2.xlsx" "~$tab_5.xlsx" "~$tab_7.xlsx"
## [5] "tab_1.xlsx" "tab_2.xlsx" "tabe_3.xlsx" "tab_4.xlsx"
## [9] "tab_5.xlsx" "tab_6.xlsx" "tab_7.xlsx"
fullpath_1 <- file.path(dir, "tab_1.xlsx")
## had to earn these since the list of patients in two data were different
## patients with adc
adc_patients <- as.data.frame(read_excel(fullpath_1, sheet = 2)) %>%
  filter((`Histology Type` == 'ADC')) %>%
  pull(ID)
##snv column name
snv_columns <- read_excel(fullpath_1, sheet = 4, skip = 1, n_max = 0) %>%
  names()
data_1_female <- as.data.frame(read_excel(fullpath_1, sheet = 2)) %>%
  filter(`Histology Type` == 'ADC') %>%
  filter(ID %in% snv_columns & Gender == 'Female') %>%
  pull(ID)

#snv data of adc patients
data_1_snv <- as.data.frame(read_excel(fullpath_1, sheet = 4, skip = 1)) %>%
  select(Gene, which(snv_columns %in% adc_patients))
rbm10_snv_patient <- as.data.frame(read_excel(fullpath_1, sheet = 4, skip = 1)) %>%
  select(Gene, which(snv_columns %in% adc_patients)) %>%
  filter(Gene == "RBM10") %>%
  pivot_longer(cols = -Gene, names_to = "patient", values_to = "RBM10 snv") %>%
  mutate(have_rbm10_snv = ifelse(`RBM10 snv` %in% NA, "N", "Y")) %>%
  filter(`have_rbm10_snv` == "Y") %>% pull(patient)

## transcriptsome data (log2T/N) of ADC patients
transcriptsome_columns <- read_excel(fullpath_1, sheet = 5, n_max = 0) %>% names()
data_1_transcriptome <- as.data.frame(read_excel(fullpath_1, sheet = 5)) %>% 
  select(gene, which(transcriptsome_columns %in% adc_patients & transcriptsome_columns %in% snv_columns))

## proteome data (log2T/N) of ADC patients
proteome_columns <- read_excel(fullpath_1, sheet = 6, na="NA", n_max = 0) %>%
  names()
data_1_proteome <- as.data.frame(read_excel(fullpath_1, sheet = 6, na = "NA")) %>%
  select(Gene, which(proteome_columns %in% adc_patients & proteome_columns %in% snv_columns))
## egfr
### egfr gene
rbm10_cancer_egfr <- c("RBM10", "EGFR")
### transcriptome
egfr_transcriptome <- data_1_transcriptome %>% filter(gene %in% rbm10_cancer_egfr) %>%
  pivot_longer(cols = -gene, names_to = "patient", values_to = "transcriptome")%>%
  mutate(with_rbm10_snv = ifelse(patient %in% rbm10_snv_patient, "Y", "N")) %>%
  mutate(gender = ifelse(patient %in% data_1_female, "F", "M")) %>%
  mutate(groups = case_when(gender == "F"& with_rbm10_snv == "Y" ~ "Female with RBM10 SNV",
                            gender == "F"& with_rbm10_snv == "N" ~ "Female without RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "Y" ~ "Male with RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "N" ~ "Male without RBM10 SNV")) %>%
  select(-c(with_rbm10_snv, gender)) %>%
  pivot_wider( names_from = "gene", values_from = "transcriptome")
mylabels <- c(`1` = "Female with RBM10 SNV",
              `2` = "Female without RBM10 SNV",
              `3` = "Male with RBM10 SNV",
              `4` = "Male without RBM10 SNV")
plot_egfr_trans <- egfr_transcriptome %>% ggparcoord(
  columns = 3:4,
  groupColumn = NULL,
  showPoints = TRUE,
  alphaLines = 1,
  boxplot = TRUE) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(limits = rbm10_cancer_egfr) +
  theme_bw() +
  labs(y = 'Log2T/N',
       x = NULL,
       title = "RBM10 and EGFR RNA Expression")+
  facet_wrap(. ~ groups, ncol = 4, labeller = labeller(groups = mylabels))

### proteome
egfr_proteome <- data_1_proteome %>% filter(Gene %in% rbm10_cancer_egfr) %>%
  pivot_longer(cols = -Gene, names_to = "patient", values_to = "proteome") %>%
  mutate(with_rbm10_snv = ifelse(patient %in% rbm10_snv_patient, "Y", "N")) %>%
  mutate(gender = ifelse(patient %in% data_1_female, "F", "M")) %>%
  mutate(groups = case_when(gender == "F"& with_rbm10_snv == "Y" ~ "Female with RBM10 SNV",
                            gender == "F"& with_rbm10_snv == "N" ~ "Female without RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "Y" ~ "Male with RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "N" ~ "Male without RBM10 SNV")) %>%
  select(-c(with_rbm10_snv, gender)) %>%
  mutate_if(is.numeric, round, 3) %>%
  pivot_wider( names_from = "Gene", values_from = "proteome")                         
                                         
plot_egfr_prot <- egfr_proteome %>% ggparcoord(
  columns = 3:4,
  groupColumn = NULL,
  showPoints = TRUE,
  alphaLines = 1,
  boxplot = TRUE) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(limits = rbm10_cancer_egfr) +
  theme_bw()+
  labs(y = 'Log2T/N',
       x = NULL,
       title = "RBM10 and EGFR Protein Expression")+
  facet_wrap(. ~ groups, ncol = 4,labeller = labeller(groups = mylabels))

## mapk
### mapk pathway genes
rbm10_cancer_mapk <- c("RBM10", "AKT2", "AKT3", "BAD", "FOXO3", "MAP2K2", "MAPK1", "MAPK3", "PLCG1", "PLCG2")

### filtering basic data and changing into tidy form
mapk_transcriptome <- data_1_transcriptome %>% filter(gene %in% rbm10_cancer_mapk) %>%

pivot_longer(cols = -gene, names_to = "patient", values_to = "transcriptome")
mapk_proteome <- data_1_proteome %>% filter(Gene %in% rbm10_cancer_mapk) %>%
pivot_longer(cols = -Gene, names_to = "patient", values_to = "proteome")

#### unifying column name
names(mapk_proteome)[names(mapk_proteome) == 'Gene'] <- 'gene'

### merging data and mutating into four groups
mapk <- merge(mapk_transcriptome, mapk_proteome, by= c("gene","patient"), all = TRUE) %>%
  pivot_longer(cols = -c(gene, patient), names_to = "data_type", values_to = "value") %>%
  mutate(with_rbm10_snv = ifelse(patient %in% rbm10_snv_patient, "Y", "N")) %>%
  mutate(gender = ifelse(patient %in% data_1_female, "F", "M")) %>%
  mutate(groups = case_when(gender == "F"& with_rbm10_snv == "Y" ~ "Female with RBM10 SNV",
                            gender == "F"& with_rbm10_snv == "N" ~ "Female without RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "Y" ~ "Male with RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "N" ~ "Male without RBM10 SNV"))

### plot
plot_mapk <- mapk%>% ggplot(aes(value, gene, fill = groups)) +
  geom_density_ridges(alpha = 0.5) +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_ridges(grid = FALSE) +
  scale_y_discrete(limits = rbm10_cancer_mapk) +
  scale_color_brewer(palette = "Set1") +
  labs(x = 'Log2T/N',
       y = 'MAPK pathway genes',
       title = "RBM10 and MAPK Pathway Gene Expression") +
  theme_bw() +
  facet_grid(data_type ~ .) +
  theme(plot.title = element_text(size=11),
        legend.position = "none")   

## PI3K
rbm10_cancer_pi3k <- c("RBM10", "AKT", "MTOR", "PTEN", "FRAP", "FRAP1", "FRAP2", "RAFT1", "RAPT1")

### filtering basic data and changing into tidy form
pi3k_transcriptome <- data_1_transcriptome %>% filter(gene %in% rbm10_cancer_pi3k) %>%
  pivot_longer(cols = -gene, names_to = "patient", values_to = "transcriptome")
pi3k_proteome <- data_1_proteome %>% filter(Gene %in% rbm10_cancer_pi3k) %>%
  pivot_longer(cols = -Gene, names_to = "patient", values_to = "proteome")

#### unifying column name
names(pi3k_proteome)[names(pi3k_proteome) == 'Gene'] <- 'gene'

### merging data and mutating into four groups
pi3k <- merge(pi3k_transcriptome, pi3k_proteome, by= c("gene","patient"), all = TRUE) %>%
  pivot_longer(cols = -c(gene, patient), names_to = "data_type", values_to = "value") %>%
  mutate(with_rbm10_snv = ifelse(patient %in% rbm10_snv_patient, "Y", "N")) %>%
  mutate(gender = ifelse(patient %in% data_1_female, "F", "M")) %>%
  mutate(groups = case_when(gender == "F"& with_rbm10_snv == "Y" ~ "Female with RBM10 SNV",
                            gender == "F"& with_rbm10_snv == "N" ~ "Female without RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "Y" ~ "Male with RBM10 SNV",
                            gender == "M"& with_rbm10_snv == "N" ~ "Male without RBM10 SNV"))

### plot
plot_pi3k <- pi3k%>% ggplot(aes(value, gene, fill = groups)) +
  geom_density_ridges(alpha = 0.5) +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_ridges(grid = FALSE) +
  scale_y_discrete(limits = rbm10_cancer_pi3k) +
  scale_color_brewer(palette = "Set1") +
  labs(x = 'Log2T/N',
       y = 'PI3K pathway genes',
       title = "RBM10 and PI3K Pathway Gene Expression") +
  theme_bw() +
  facet_grid(data_type ~ .) +
  theme(plot.title = element_text(size=11),
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.title = element_text(size=9),
        legend.text = element_text(size=9),
        legend.position = c(.5,.9),)

# all of four plots
plots_top <- plot_grid(
  plot_mapk, plot_pi3k,
  labels = "AUTO", ncol = 2)
## Picking joint bandwidth of 0.112
## Picking joint bandwidth of 0.24
## Warning: Removed 26 rows containing non-finite values (stat_density_ridges).
## Picking joint bandwidth of 0.0968
## Picking joint bandwidth of 0.264
## Warning: Removed 3 rows containing non-finite values (stat_density_ridges).
plots_bottom <- plot_grid(
  plot_egfr_trans, plot_egfr_prot, labels = c("C", "D"), ncol = 1)
title <- ggdraw() +
  draw_label(
    "RBM10 and Expression of Possible Cancer Produing Pathways",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 15
  ) +
  theme(plot.margin = margin(5, 5, 5, 7)
  )

figure_1 <- plot_grid(
  title, plots_top, plots_bottom,
  ncol = 1,
  rel_heights = c(0.1, 6, 3 )
)
ggsave("figure_1.png", figure_1, width = 8, height = 14)


# The Figure is attached as in separate files.

