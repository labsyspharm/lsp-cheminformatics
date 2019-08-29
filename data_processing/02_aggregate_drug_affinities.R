library(tidyverse)
library(httr)
library(data.table)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# set directories & import files -----------------------------------------------
###############################################################################T

all_cmpds_eq_classes <- syn("syn20693885") %>%
  read_csv(col_types = "cci")

biochem_alldata <- syn("syn20693825") %>%
  read_csv(col_types = "iiiicccdcicccccccicicc")

biochem_neat <- biochem_alldata %>%
  rename(pref_name_target = pref_name) %>%
  filter(tax_id == 9606, !is.na(gene_id)) %>%
  left_join(
    all_cmpds_eq_classes %>%
      select(id, eq_class),
    by = c("chembl_id_compound" = "id")
  )

doseresponse_inhouse <- syn("syn20692433") %>%
  read_csv(col_types = "idcccccicccc")

doseresponse_inhouse_neat <- doseresponse_inhouse %>%
  mutate(hms_id = paste0("HMSL", hms_id)) %>%
  left_join(
    all_cmpds_eq_classes %>%
      select(id, eq_class),
    by = c("hms_id" = "id")
  )

# bind files and make Q1 file --------------------------------------------------
###############################################################################T

biochem_rowbind <- biochem_neat %>%
  mutate(
    reference_id = chembl_id_doc,
    reference_type = "chembl_id",
    file_url = paste0("https://www.ebi.ac.uk/chembl/document_report_card/",chembl_id_doc),
    value = standard_value,
    value_unit = standard_units
  ) %>%
  select(eq_class, value, value_unit, gene_id, reference_id, reference_type, file_url)

doseresponse_inhouse_rowbind <- doseresponse_inhouse_neat %>%
  mutate(
    reference_id = synapse_id,
    reference_type = "synapse_id"
  ) %>%
  select(eq_class, value, value_unit, gene_id, reference_id, reference_type, file_url)

complete_table <- bind_rows(
  biochem_rowbind,
  doseresponse_inhouse_rowbind
) %>%
  # Call to distinct important, since some assays can be recorded multiple times
  # for the same eq_class now, when multiple forms of the same drug where mapped
  # to the same eq_class and an assay was stored in the db for all forms
  distinct()

write_csv(
  complete_table,
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl.csv.gz")
)

#test: unit
complete_table$value_unit%>%unique()

# Using data.table here for speed
complete_table_Q1 <- complete_table %>%
  data.table::as.data.table() %>%
  .[
    ,
    .(Q1 = round(quantile(value, 0.25, names = FALSE), 2)),
    by = .(eq_class, gene_id)
  ] %>%
  as_tibble() %>%
  mutate(binding = if_else(Q1 < 10000, 1L, 0L))

write_csv(
  complete_table_Q1,
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl_Q1.csv.gz")
)

# Also calculate Q1 values for kinomescan data from HMS LINCS for which no
# complete dose response curve is available
hmsl_kinomescan <- syn("syn20692432") %>%
  read_csv(col_types = "cccdcdcc") %>%
  # Some of the HMSL IDs don't start with "HMSL", fix that
  mutate(hms_id = if_else(str_starts(hms_id, "HMSL"), hms_id, paste0("HMSL", hms_id)))

# check how common the situation is where a single gene was tested in different
# variants (mutants, post-translational modification, etc.)
hmsl_kinomescan %>%
  count(hms_id, gene_symbol, cmpd_conc_nM) %>%
  count(n)
# # A tibble: 21 x 2
# n    nn
# <int> <int>
#   1     1 60448
# 2     2  4467
# 3     3   283
# 4     4   250
# 5     5    32
# 6     6    57
# 7     7    56
# 8     8    88
# 9     9   122
# 10    10   161
# # … with 11 more rows
# It's very common, so we have to decide how to aggregate info for every
# target/compound combo.
# I think it makes sense to try to filter out mutant data, but keep post-translational
# modified targets (especially phosphorylation), because those are likely to be
# broadly relevant. if somebody is interested in specific mutations they will
# have to query the Kinomescan directly and shoulnd't rely on TAS.

hmsl_kinomescan_cleaned <- hmsl_kinomescan %>%
  filter(
    # Filter mutant annotation e.g. S154A
    !grepl("[^A-Za-z]+[A-Z][0-9]+([A-Z]|del|Del)[^A-Za-z]", description),
    !grepl("domain", description),
    !grepl("Dom", description),
    !grepl("inhibited", description)
  )

hmsl_kinomescan_mapped <- hmsl_kinomescan_cleaned %>%
  genebabel::join_hgnc("gene_symbol", c("symbol", "alias_symbol", "prev_symbol"), c("entrez_id", "name")) %>%
  left_join(
    all_cmpds_eq_classes %>%
      select(id, eq_class),
    by = c("hms_id" = "id")
  ) %>%
  rename(gene_id = entrez_id, pref_name_cmpd = pref_name, pref_name_target = name) %>%
  drop_na(gene_id)

write_csv(
  hmsl_kinomescan_mapped,
  file.path(dir_release, "biochemicaldata_single_dose_inhouse.csv.gz")
)

hmsl_kinomescan_q1 <- hmsl_kinomescan_mapped %>%
  as.data.table() %>%
  .[
    ,
    .(percent_control_Q1 = quantile(percent_control, 0.25, names = FALSE)),
    by = .(eq_class, gene_id, cmpd_conc_nM)
    ] %>%
  as_tibble()

write_csv(
  hmsl_kinomescan_q1,
  file.path(dir_release, "biochemicaldata_single_dose_inhouse_Q1.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

aggregation_activity <- Activity(
  name = "Aggregate affinity data",
  used = c(
    "syn20693885",
    "syn20693825",
    "syn20692433",
    "syn20692432"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/02_aggregate_drug_affinities.R"
)

list(
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl.csv.gz"),
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl_Q1.csv.gz"),
  file.path(dir_release, "biochemicaldata_single_dose_inhouse.csv.gz"),
  file.path(dir_release, "biochemicaldata_single_dose_inhouse_Q1.csv.gz")
) %>%
  map(
    . %>%
      File(parent = syn_release) %>%
      synStore(activity = aggregation_activity)
  )
