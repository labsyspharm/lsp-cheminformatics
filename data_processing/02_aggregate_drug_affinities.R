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

all_cmpds_eq_classes <- syn("syn20830516") %>%
  read_rds()

biochem_alldata <- syn("syn20693825") %>%
  read_csv(col_types = "iiiicccdciccccicccccii")

biochem_neat <- all_cmpds_eq_classes %>%
  mutate(
    data = map(
      data,
      ~biochem_alldata %>%
        rename(pref_name_target = pref_name) %>%
        filter(tax_id == 9606, !is.na(entrez_gene_id)) %>%
        left_join(
          .x %>%
            select(id, eq_class),
          by = c("chembl_id_compound" = "id")
        )
    )
  )

doseresponse_inhouse <- syn("syn20692433") %>%
  read_csv(col_types = "idcccccicccc")

doseresponse_inhouse_neat <- all_cmpds_eq_classes %>%
  mutate(
    data = map(
      data,
      ~doseresponse_inhouse %>%
        mutate(hms_id = paste0("HMSL", hms_id)) %>%
        rename(entrez_gene_id = gene_id) %>%
        left_join(
          .x %>%
            select(id, eq_class),
          by = c("hms_id" = "id")
        )
    )
  )

pheno_data <- syn("syn20693827") %>%
  read_csv(col_types = "iiiicccdciccccd")

pheno_data_neat <- all_cmpds_eq_classes %>%
  mutate(
    data = map(
      data,
      ~pheno_data %>%
        left_join(
          .x %>%
            select(id, eq_class),
          by = c("chembl_id_compound" = "id")
        )
    )
  )

# Aggregate dose-response data -------------------------------------------------
###############################################################################T

biochem_rowbind <- biochem_neat %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
        reference_id = chembl_id_doc,
        reference_type = "chembl_id",
        file_url = paste0("https://www.ebi.ac.uk/chembl/document_report_card/",chembl_id_doc),
        value = standard_value,
        value_unit = standard_units
      ) %>%
        select(
          lspci_id = eq_class,
          value, value_unit = standard_units, value_type = standard_type,
          value_relation = standard_relation, description_assay = description,
          uniprot_id, entrez_gene_id,
          reference_id, reference_type, file_url
        )
    )
  )


doseresponse_inhouse_rowbind <- doseresponse_inhouse_neat %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
          reference_id = synapse_id,
          reference_type = "synapse_id"
        ) %>%
        select(
          lspci_id = eq_class,
          value, value_unit, value_type,
          value_relation, description_assay = description,
          uniprot_id, entrez_gene_id,
          reference_id, reference_type, file_url
        )
    )
  )

complete_table <- biochem_rowbind %>%
  rename(biochem = data) %>%
  left_join(
    doseresponse_inhouse_rowbind %>%
      rename(inhouse = data)
  ) %>%
  mutate(
    # Call to distinct important, since some assays can be recorded multiple times
    # for the same eq_class now, when multiple forms of the same drug where mapped
    # to the same eq_class and an assay was stored in the db for all forms
    data = map2(biochem, inhouse, ~bind_rows(.x, .y) %>% distinct())
  )

write_rds(
  complete_table %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl.rds"),
  compress = "gz"
)

#test: unit
complete_table$value_unit%>%unique()

# Using data.table here for speed
calculate_q1 <- function(data) {
  data %>%
    data.table::as.data.table() %>%
    .[
      ,
      .(Q1 = round(quantile(value, 0.25, names = FALSE), 2)),
      by = .(eq_class, entrez_gene_id)
      ] %>%
    as_tibble() %>%
    mutate(binding = if_else(Q1 < 10000, 1L, 0L))
}

complete_table_Q1 <- complete_table %>%
  mutate(
    data = map(
      data,
      calculate_q1
    )
  )

write_rds(
  complete_table_Q1 %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl_Q1.rds"),
  compress = "gz"
)

# Aggregate single-dose data ---------------------------------------------------
###############################################################################T


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

hmsl_kinomescan_mapped <- all_cmpds_eq_classes %>%
  mutate(
    data = map(
      data,
      ~hmsl_kinomescan_cleaned %>%
        genebabel::join_hgnc(
          "gene_symbol",
          c("symbol", "alias_symbol", "prev_symbol"),
          c("entrez_id", "name", "uniprot_ids")
        ) %>%
        # I checked, no gene_symbol maps to multiple uniprot, so this is safe
        mutate(
          uniprot_id = as.character(uniprot_ids)
        ) %>%
        select(-uniprot_ids) %>%
        left_join(
          .x %>%
            select(id, eq_class),
          by = c("hms_id" = "id")
        ) %>%
        rename(
          entrez_gene_id = entrez_id,
          pref_name_cmpd = pref_name,
          pref_name_target = name,
          lspci_id = eq_class
        ) %>%
        drop_na(entrez_gene_id)
    )
  )

write_rds(
  hmsl_kinomescan_mapped,
  file.path(dir_release, "biochemicaldata_single_dose_inhouse.rds"),
  compress = "gz"
)

hmsl_kinomescan_q1 <- hmsl_kinomescan_mapped %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        as.data.table() %>%
        .[
          ,
          .(percent_control_Q1 = quantile(percent_control, 0.25, names = FALSE)),
          by = .(eq_class, entrez_gene_id, cmpd_conc_nM)
          ] %>%
        as_tibble()
    )
  )


write_rds(
  hmsl_kinomescan_q1,
  file.path(dir_release, "biochemicaldata_single_dose_inhouse_Q1.rds"),
  compress = "gz"
)

# Aggregate phenotypic data ----------------------------------------------------
###############################################################################T

pheno_data_formatted <- pheno_data_neat %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
          reference_id = chembl_id_doc,
          reference_type = "chembl_id",
          file_url = paste0("https://www.ebi.ac.uk/chembl/document_report_card/", chembl_id_doc)
        ) %>%
        select(
          lspci_id = eq_class,
          value = standard_value, value_unit = standard_units, value_type = standard_type,
          value_relation = standard_relation, description_assay = description,
          reference_id, reference_type, file_url
        )
    )
  )

write_rds(
  pheno_data_formatted,
  file.path(dir_release, "pheno_data.rds"),
  compress = "gz"
)

pheno_data_q1 <- pheno_data_neat %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        as.data.table() %>%
        .[
          ,
          .(
            standard_value_Q1 = quantile(standard_value, 0.25, names = FALSE),
            log10_value_Q1 = quantile(log10_value, 0.25, names = FALSE)
          ),
          keyby = .(eq_class, assay_id)
          ] %>%
        as_tibble()
    )
  )

write_rds(
  pheno_data_q1,
  file.path(dir_release, "pheno_data_Q1.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

aggregation_activity <- Activity(
  name = "Aggregate affinity data",
  used = c(
    "syn20692432",
    "syn20692433",
    "syn20693825",
    "syn20693827",
    "syn20830516"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/02_aggregate_drug_affinities.R"
)

syn_aggregate <- Folder("aggregate_data", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl.rds"),
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl_Q1.rds"),
  file.path(dir_release, "biochemicaldata_single_dose_inhouse.rds"),
  file.path(dir_release, "biochemicaldata_single_dose_inhouse_Q1.rds"),
  file.path(dir_release, "pheno_data.rds"),
  file.path(dir_release, "pheno_data_Q1.rds")
) %>%
  synStoreMany(parent = syn_aggregate, activity = aggregation_activity)
