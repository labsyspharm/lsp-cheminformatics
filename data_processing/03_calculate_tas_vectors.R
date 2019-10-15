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

complete_dose_response_Q1 <- syn("syn20830834") %>%
  read_rds()

complete_single_dose_Q1 <- syn("syn20830836") %>%
  read_rds()

# Checking if all target/compound combinations are unique
complete_dose_response_Q1 %>%
  unnest(data) %>%
  count(fp_name, fp_type, eq_class, entrez_gene_id) %>%
  count(n)
# # A tibble: 1 x 2
# n      nn
# <int>   <int>
#   1     1 1776340

complete_single_dose_Q1 %>%
  unnest(data) %>%
  count(fp_name, fp_type, eq_class, entrez_gene_id, cmpd_conc_nM) %>%
  count(n)
# # A tibble: 1 x 2
# n     nn
# <int>  <int>
#   1     1 130796

cmpd_eq_classes <- syn("syn20830516") %>%
  read_rds()

literature_annotations_raw <- syn("syn20694521") %>%
  read_csv(col_types = "ii___________________") %>%
  mutate(hms_id = paste0("HMSL", hms_id))

literature_annotations <- cmpd_eq_classes %>%
  mutate(
    data = map(
      data,
      ~literature_annotations_raw %>%
        left_join(
          .x %>%
            select(eq_class, id),
          by = c("hms_id" = "id")
        ) %>%
        select(eq_class, entrez_gene_id = gene_id) %>%
        mutate(tas_literature = 2L) %>%
        # There are some duplicates in the literature annotation file
        # plus, in one case, two HMSL ids map to the same eq_class
        distinct()
    )
  )



# calculate TAS ----------------------------------------------------------------
###############################################################################T

complete_table_tas <- complete_dose_response_Q1 %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
          tas = case_when(
            Q1 < 100 ~ 1L,
            Q1 < 1000 ~ 2L,
            Q1 < 10000 ~ 3L,
            Q1 >= 10000 ~ 10L,
            TRUE ~ NA_integer_
          )
        )
    )
  )


hmsl_kinomescan_tas <- complete_single_dose_Q1 %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
        tas = case_when(
          cmpd_conc_nM == 10000 & percent_control_Q1 >= 75 ~ 10L,
          cmpd_conc_nM == 10000 & percent_control_Q1 < 0.1 ~ 2L,
          cmpd_conc_nM == 1000 & percent_control_Q1 >= 90 ~ 10L,
          cmpd_conc_nM == 1000 & percent_control_Q1 < 1 ~ 2L,
          cmpd_conc_nM == 100 & percent_control_Q1 >= 75 ~ 10L,
          cmpd_conc_nM == 100 & percent_control_Q1 < 25 ~ 2L,
          TRUE ~ NA_integer_
        ),
        entrez_gene_id = as.integer(entrez_gene_id)
      )
    )
  )


# Check for contradicting TAS at different compound concentrations
hmsl_kinomescan_tas %>%
  unnest(data) %>%
  as.data.table() %>%
  .[
    ,
    .(n = length(unique(na.omit(tas))) > 1),
    keyby = .(fp_name, fp_type, eq_class, entrez_gene_id)
  ] %>%
  .[
    ,
    sum(n),
    keyby = .(fp_name, fp_type)
  ]
# fp_name     fp_type V1
# 1:      morgan_chiral      morgan 42
# 2:      morgan_normal      morgan 42
# 3: topological_normal topological 42
# Thre are a few cases where there is contradicting info.
# In this case, take the minimum TAS. According to Nienke false positives
# are less likely than false negatives so this seems like the prudent approach.


hmsl_kinomescan_tas_agg <- hmsl_kinomescan_tas %>%
  mutate(
    data = map(
      data,
      ~.x %>%
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
    )
  )


combined_q1 <- complete_table_tas %>%
  rename(dose_response = data) %>%
  full_join(
    hmsl_kinomescan_tas_agg %>%
      rename(single_dose = data)
  ) %>%
  full_join(
    literature_annotations %>%
      rename(literature = data)
  ) %>%
  mutate(
    data = pmap(
      list(dose_response, single_dose, literature),
      function(dose_response, single_dose, literature) {
        full_join(
          dose_response %>%
            select(eq_class, entrez_gene_id, tas_affinity = tas),
          single_dose %>%
            select(eq_class, entrez_gene_id, tas_single_dose = tas),
          by = c("eq_class", "entrez_gene_id")
        ) %>%
          full_join(
            literature,
            by = c("eq_class", "entrez_gene_id")
          ) %>%
          arrange(eq_class, entrez_gene_id)
      }
    )
  ) %>%
  select(fp_name, fp_type, data)




# Again checking for cases where we get contradictory results
# combined_q1 %>%
#   filter(tas_affinity != tas_single_dose)
#
# combined_q1 %>%
#   filter(tas_affinity > tas_literature) %>%
#   View()

# Prefer results from full affinity measuremennts over the percent inhibition
# When incorporating Verena's manual annotatations, prefer affinity data.
# If affinity data not present, take minimum of of Verena + percent control
combined_q1_agg <- combined_q1 %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
        tas = if_else(
          !is.na(tas_affinity),
          tas_affinity,
          pmin(tas_single_dose, tas_literature, na.rm = TRUE)
        )
      )
    )
  )


gene_map <- syn("syn20693721") %>%
  read_csv(col_types = "ciciccciccccc")

compound_map <- syn("syn20835543") %>%
  read_rds()

tas_vector <- combined_q1_agg %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        select(lspci_id = eq_class, entrez_gene_id, starts_with("tas"))
    )
  )


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
      rename(compound_map = data)
  ) %>%
  mutate(
    data = map2(
      data, compound_map,
      ~.x %>%
        left_join(
          .y %>%
            select(lspci_id, chembl_id, hms_id, pref_name),
          by = "lspci_id"
        ) %>%
        left_join(
          gene_info,
          by = "entrez_gene_id"
        )
    )
  ) %>%
  select(fp_name, fp_type, data)


tas_vector_annotated_long <- tas_vector %>%
  left_join(
    cmpd_eq_classes %>%
      rename(cmpd_eq_classes = data)
  ) %>%
  mutate(
    data = map2(
      data, cmpd_eq_classes,
      ~.x %>%
        left_join(
          .y %>%
            select(compound_id = id, compound_id_source = source, lspci_id = eq_class),
          by = "lspci_id"
        )%>%
        left_join(
          gene_info %>%
            select(entrez_gene_id, entrez_symbol),
          by = "entrez_gene_id"
        ) %>%
        mutate(
          entrez_symbol = if_else(!is.na(entrez_symbol), entrez_symbol, as.character(entrez_gene_id))
        ) %>%
        select(
          lspci_id,
          compound_id,
          compound_id_source,
          entrez_gene_id,
          entrez_symbol,
          starts_with("tas")
        )
    )
  ) %>%
  select(fp_name, fp_type, data)

write_rds(
  tas_vector,
  file.path(dir_release, "tas_vector.rds"),
  compress = "gz"
)

write_rds(
  tas_vector_annotated,
  file.path(dir_release, "tas_vector_annotated.rds"),
  compress = "gz"
)

write_rds(
  tas_vector_annotated_long,
  file.path(dir_release, "tas_vector_annotated_long.rds"),
  compress = "gz"
)

write_csv(
  tas_vector_annotated_long %>%
    unnest(data),
  file.path(dir_release, "tas_vector_annotated_long.csv.gz"),
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

tas_activity <- Activity(
  name = "Calculate target affinity spectrum (TAS) vectors",
  used = c(
    "syn20830834",
    "syn20830836",
    "syn20830516",
    "syn20694521",
    "syn20693721",
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/03_calculate_tas_vectors.R"
)

syn_tas_folder <- Folder("tas", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "tas_vector.rds"),
  file.path(dir_release, "tas_vector_annotated.rds"),
  file.path(dir_release, "tas_vector_annotated_long.rds"),
  file.path(dir_release, "tas_vector_annotated_long.csv.gz")
) %>%
  synStoreMany(parentId = syn_tas_folder, activity = tas_activity)
