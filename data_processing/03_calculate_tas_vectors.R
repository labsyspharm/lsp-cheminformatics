# Script to calculate TAS vectors from all combined data sources

library(tidyverse)
library(data.table)
library(here)
library(vroom)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# set directories & import files -----------------------------------------------
###############################################################################T

complete_table_Q1 <- syn("syn20693968") %>%
  read_csv(col_types = "iidi")

hmsl_kinomescan <- syn("syn20693970") %>%
  read_csv(col_types = "iidd")

# Checking if all target/compound combinations are unique
complete_table_Q1 %>%
  count(eq_class, entrez_gene_id) %>%
  count(n)
# # A tibble: 1 x 2
# n     nn
# <int>  <int>
#   1     1 877496

hmsl_kinomescan %>%
  count(eq_class, entrez_gene_id, cmpd_conc_nM) %>%
  count(n)
# # A tibble: 1 x 2
# n    nn
# <int> <int>
#   1     1 65398

cmpd_eq_classes <- syn("syn20693885") %>%
  read_csv(col_types = "cci")

literature_annotations <- syn("syn20694521") %>%
  read_csv(col_types = "ii___________________") %>%
  mutate(hms_id = paste0("HMSL", hms_id)) %>%
  left_join(
    cmpd_eq_classes %>%
      select(eq_class, id),
    by = c("hms_id" = "id")
  ) %>%
  select(eq_class, entrez_gene_id = gene_id)

# calculate TAS ----------------------------------------------------------------
###############################################################################T

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
  as.data.table() %>%
  .[
    ,
    .(n = length(unique(na.omit(tas))) > 1),
    keyby = .(eq_class, entrez_gene_id)
  ] %>%
  .[
    ,
    sum(n)
  ]
# [1] 42
# Thre are a few cases where there is contradicting info.
# In this case, take the minimum TAS. According to Nienke false positives
# are less likely than false negatives so this seems like the prudent approach.


hmsl_kinomescan_tas_agg <- hmsl_kinomescan_tas %>%
  as.data.table() %>%
  .[
    ,
    .(tas = if(all(is.na(tas))) NA_integer_ else min(tas, na.rm = TRUE)),
    keyby = .(eq_class, entrez_gene_id)
  ] %>%
  as_tibble() %>%
  # Drop any NAs here, sometimes no TAS assertion can be made
  # if percent_control value is in a certain range
  drop_na(tas)

combined_q1 <- full_join(
  complete_table_tas %>%
    select(eq_class, entrez_gene_id, tas_affinity = tas),
  hmsl_kinomescan_tas_agg %>%
    select(eq_class, entrez_gene_id, tas_inhibition = tas),
  by = c("eq_class", "entrez_gene_id")
) %>%
  full_join(
    literature_annotations %>%
      mutate(tas_literature = 2L),
    by = c("eq_class", "entrez_gene_id")
  ) %>%
  arrange(eq_class, entrez_gene_id)

# Again checking for cases where we get contradictory results
combined_q1 %>%
  drop_na() %>%
  filter(tas_affinity != tas_inhibition)

# Prefer results from full affinity measuremennts over the percent inhibition
# When incorporating Verena's manual annotatations, prefer affinity data.
# If affinity data not present, take minimum of of Verena + percent control
combined_q1_agg <- combined_q1 %>%
  mutate(
    tas_combined = if_else(
      !is.na(tas_affinity),
      tas_affinity,
      pmin(tas_inhibition, tas_literature, na.rm = TRUE)
    )
  )

gene_map <- syn("syn20693721") %>%
  read_csv(col_types = "ciciccciccccc")

compound_map <- syn("syn20692551") %>%
  read_rds()

tas_vector <- combined_q1_agg %>%
  select(eq_class, entrez_gene_id, tas_combined)

# Using vroom here instead of loading the entire csv because it is downright massive
# and vroom is much faster
gene_info <- vroom(
  file.path(dir_release, "gene_info_20190829.gz"),
  delim = "\t",
  col_names = c(
    "tax_id", "entrez_gene_id", "entrez_symbol", "locus_tag", "entrez_synonyms", "db_xrefs", "chromosome",
    "map_location", "entrez_description", "entrez_type_of_gene", "symbol", "entrez_name", "nomenclature_status",
    "other_designations", "modification_date", "feature_type"
  ),
  col_types = "iic_c___cc_c____",
  skip = 1
)

tas_vector_annotated <- tas_vector %>%
  left_join(
    compound_map %>%
      select(eq_class = lspci_id, chembl_id, hms_id, pref_name),
    by = "eq_class"
  ) %>%
  left_join(
    gene_info,
    by = "entrez_gene_id"
  )

tas_vector_annotated_long <- tas_vector %>%
  left_join(
    cmpd_eq_classes %>%
      select(compound_id = id, compound_id_source = source, eq_class),
    by = "eq_class"
  ) %>%
  left_join(
    gene_info %>%
      select(entrez_gene_id, entrez_symbol),
    by = "entrez_gene_id"
  ) %>%
  select(
    eq_class,
    compound_id,
    compound_id_source,
    entrez_gene_id,
    entrez_symbol,
    tas = tas_combined
  )

write_csv(
  tas_vector,
  file.path(dir_release, "tas_vector.csv.gz")
)

write_csv(
  tas_vector_annotated,
  file.path(dir_release, "tas_vector_annotated.csv.gz")
)

write_csv(
  tas_vector_annotated_long,
  file.path(dir_release, "tas_vector_annotated_long.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

tas_activity <- Activity(
  name = "Calculate target affinity spectrum (TAS) vectors",
  used = c(
    "syn20693968",
    "syn20693970",
    "syn20693885",
    "syn20694521",
    "syn20693721",
    "syn20692551"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/03_calculate_tas_vectors.R"
)

c(
  file.path(dir_release, "tas_vector.csv.gz"),
  file.path(dir_release, "tas_vector_annotated.csv.gz"),
  file.path(dir_release, "tas_vector_annotated_long.csv.gz")
) %>%
  synStoreMany(parentId = syn_release, activity = tas_activity)
