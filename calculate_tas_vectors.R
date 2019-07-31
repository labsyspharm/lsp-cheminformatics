# Script to calculate TAS vectors from all combined data sources

library(tidyverse)

dir_cheminformatics_files = "~/repo/tas_vectors/chembl24_1"
dir_hmsl_data = "~/repo/tas_vectors"

setwd(dir_cheminformatics_files)
complete_table_Q1 <- read_csv("Q1_table_alldata_20190708.csv.gz", col_types = "cid")

setwd(dir_hmsl_data)
hmsl_kinomescan <- read_csv("hmsl_kinomescan_q1_20190619.csv.gz", col_types = "cidd")

# Checking if all target/compound combinations are unique
complete_table_Q1 %>%
  count(chembl_id_cmpd, gene_id) %>%
  count(n)
# # A tibble: 1 x 2
# n     nn
# <int>  <int>
#   1     1 952128

hmsl_kinomescan %>%
  count(chembl_id, gene_id, cmpd_conc_nM) %>%
  count(n)
# # A tibble: 1 x 2
# n    nn
# <int> <int>
#   1     1 64784


complete_table_tas <- complete_table_Q1 %>%
  mutate(
    tas = case_when(
      Q1 < 100 ~ 1L,
      Q1 < 1000 ~ 2L,
      Q1 < 10000 ~ 3L,
      Q1 >= 10000 ~ 10L,
      TRUE ~ NA_integer_
    )
  )

hmsl_kinomescan_tas <- hmsl_kinomescan %>%
  mutate(
    tas = case_when(
      cmpd_conc_nM == 10000 & percent_control_Q1 >= 75 ~ 10L,
      cmpd_conc_nM == 10000 & percent_control_Q1 < 0.1 ~ 2L,
      cmpd_conc_nM == 1000 & percent_control_Q1 >= 90 ~ 10L,
      cmpd_conc_nM == 1000 & percent_control_Q1 < 1 ~ 2L,
      cmpd_conc_nM == 100 & percent_control_Q1 >= 75 ~ 10L,
      cmpd_conc_nM == 100 & percent_control_Q1 < 25 ~ 2L,
      TRUE ~ NA_integer_
    )
  )

# Check for contradicting TAS at different compound concentrations
hmsl_kinomescan_tas %>%
  arrange(chembl_id, gene_id) %>%
  group_by(chembl_id, gene_id) %>%
  filter(
    length(unique(na.omit(tas))) > 1
  ) %>%
  n_groups()
# [1] 52
# Thre are quite a few cases where there is contradicting info.
# In this case, take the minimum TAS. According to Nienke false positives
# are less likely than false negatives so this seems like the prudent approach.


hmsl_kinomescan_tas_agg <- hmsl_kinomescan_tas %>%
  group_by(chembl_id, gene_id) %>%
  summarize(tas = if(all(is.na(tas))) NA_integer_ else min(tas, na.rm = TRUE)) %>%
  ungroup() %>%
  drop_na()

combined_q1 <- full_join(
  complete_table_tas %>%
    select(chembl_id_cmpd, gene_id, tas_affinity = tas),
    # group_nest(chembl_id_cmpd, gene_id, .key = "kd_data"),
  hmsl_kinomescan_tas_agg %>%
    select(chembl_id_cmpd = chembl_id, gene_id, tas_inhibition = tas),
  by = c("chembl_id_cmpd", "gene_id")
) %>%
  arrange(chembl_id_cmpd, gene_id)

# Again checking for cases where we get contradictory results
combined_q1 %>%
  drop_na() %>%
  filter(tas_affinity != tas_inhibition)

# Prefer results from full affinity measuremennts over the percent inhibition
combined_q1_agg <- combined_q1 %>%
  mutate(tas_combined = ifelse(is.na(tas_affinity), tas_inhibition, tas_affinity))

# When incorporating Verena's manual annotatations, prefer affinity data.
# If affinity data not present, take minimum of of Verena + percent control


setwd(dir_cheminformatics_files)

gene_map <- read_csv(
  "map_targets_chemblID_genID_uniprot_tid.csv.gz",
  col_types = "ccccciciiccccc"
)

compound_map <- read_csv(
  "all_compounds_chembl24_1_parent_mapped.csv.gz",
  col_types = "icciciiccic"
)

hmsl_map <- read_csv(
  "hms_lincs_chembl_mapping.csv.gz",
  col_types = "ccccc"
)

tas_vector <- combined_q1_agg %>%
  select(chembl_id_cmpd, gene_id, tas_combined) %>%
  # left_join(
  #   gene_map %>%
  #     select(uniprot_id, gene_name = pref_name, gene_symbol = symbol, gene_id),
  #   by = "gene_id"
  # )
  mutate(gene_id = as.character(gene_id)) %>%
  genebabel::join_hgnc("gene_id", "entrez_id", c("symbol", "name")) %>%
  left_join(
    compound_map %>%
      distinct(chembl_id_cmpd = chembl_id, compound_name = pref_name),
    by = "chembl_id_cmpd"
  ) %>%
  rename()

tas_vector_lincs <- tas_vector %>%
  left_join(
    hmsl_map %>%
      distinct(hms_id, chembl_id_cmpd = chembl_id),
    by = "chembl_id_cmpd"
  )

write_csv(tas_vector_lincs, "tas_vector_20190720.csv.gz")
