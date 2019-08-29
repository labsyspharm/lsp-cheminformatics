#this script prepares allbiochem data for BAFP clustering

library(tidyverse)
library(httr)
library(data.table)

#############################################################################################################T
# set directories & import files----------
#############################################################################################################T
dir_chembl <- "~/repo/tas_vectors/chembl24_1"
dir_klaegerdata<-"~/repo/tas_vectors"
dir_doseresponse_inhouse<-"~/repo/tas_vectors"
dir_cheminformatics_files<-"~/repo/tas_vectors/chembl24_1"
dir_hmsl_data = "~/repo/tas_vectors"

setwd(dir_chembl)
hms_lincs_chembl_mapping <- read_csv(
  "hms_lincs_chembl_mapping.csv.gz",
  col_types = "ccccc"
)

all_cmpds <- read_csv(
  "all_compounds_chembl24_1_parent_mapped.csv.gz",
  col_types = "icciciiccic"
)

setwd(dir_klaegerdata)
klaegerdata <- read_csv(
  "Klaeger_reformatted.csv",
  col_types = "ccicdcccccc"
)
klaeger_neat <- klaegerdata %>%
  filter(chembl_id != "")

setwd(dir_chembl)
#fragment_table<-read.csv("BAFP_chemblID_fragmentID_20190707.csv")
biochem_alldata <- read_csv(
  "biochemicaldata_allcmpds_chembl24_1.csv.gz",
  col_types = "iiiiccdciccccccicicc"
)

biochem_neat <- biochem_alldata %>%
  rename(pref_name_target = pref_name) %>%
  filter(tax_id == 9606, !is.na(gene_id)) %>%
  left_join(
    all_cmpds %>%
      distinct(molregno, pref_name_cmpd = pref_name, chembl_id),
    by = "molregno"
  )

setwd(dir_doseresponse_inhouse)
doseresponse_inhouse <- read_csv(
  "all_doseresponse_summaryfile_20190417_vF.csv",
  col_types = "idcccccicccc"
)
doseresponse_inhouse_neat <- doseresponse_inhouse %>%
  mutate(hms_id = paste0("HMSL", hms_id)) %>%
  left_join(
    hms_lincs_chembl_mapping %>%
      distinct(hms_id, chembl_id),
    by = "hms_id"
  ) %>%
  left_join(
    all_cmpds %>%
      distinct(molregno,  pref_name_cmpd = pref_name, chembl_id),
    by = "chembl_id"
  ) %>%
  unique()

#############################################################################################################T
# bind files and make Q1 file----------
#############################################################################################################T

names(klaeger_neat)
names(biochem_neat)
names(doseresponse_inhouse_neat)
klaeger_rowbind<-klaeger_neat%>%
  mutate(chembl_id_cmpd=chembl_id)%>%
  dplyr::select(chembl_id_cmpd, value, value_unit, gene_id, reference_id, reference_type, file_url)
biochem_rowbind<-biochem_neat%>%
  mutate(reference_id=chembl_id_doc,
         reference_type="chembl_id",
         file_url=paste0("https://www.ebi.ac.uk/chembl/document_report_card/",chembl_id_doc),
         chembl_id_cmpd=chembl_id,
         value=standard_value,
         value_unit=standard_units)%>%
  dplyr::select(chembl_id_cmpd, value, value_unit, gene_id, reference_id, reference_type, file_url)
doseresponse_inhouse_rowbind<-doseresponse_inhouse_neat%>%
  mutate(chembl_id_cmpd=chembl_id,
         reference_id=synapse_id,
         reference_type="synapse_id")%>%
  dplyr::select(chembl_id_cmpd,value,value_unit, gene_id, reference_id,reference_type,file_url)

complete_table <- bind_rows(
  klaeger_rowbind,
  biochem_rowbind,
  doseresponse_inhouse_rowbind
) %>%
  as_tibble() %>%
  distinct() %>%
  mutate(value = round(value, 2)) %>%
  drop_na(chembl_id_cmpd, gene_id) %>%
  filter(value > 0)

setwd(dir_cheminformatics_files)
write_csv(complete_table, "complete_inhouse_klaeger_chembl_biochem_20190708.csv.gz")

#test: unit
complete_table$value_unit%>%unique()

# Using data.table here for speed
complete_table_Q1 <- complete_table %>%
  data.table::as.data.table() %>%
  .[
    ,
    .(Q1 = round(quantile(value, 0.25, names = FALSE))),
    by = .(chembl_id_cmpd, gene_id)
  ] %>%
  as_tibble()

write_csv(complete_table_Q1, "Q1_table_alldata_20190708.csv.gz")

complete_table_binding <- complete_table_Q1 %>%
  mutate(binding = ifelse(Q1 < 10000, 1, 0))

write_csv(complete_table_binding, "binding_table_alldata_20190708.csv.gz")


# Also calculate Q1 values for kinomescan data from HMS LINCS for which no
# complete dose response curve is available
setwd(dir_hmsl_data)
hmsl_kinomescan <- read_csv(
  "all_discoverX_kinomescan_20190619.csv",
  col_types = "cccdcdcc"
)

# check how common the situation is where a single gene was tested in different
# variants (mutants, post-translational modification, etc.)
hmsl_kinomescan %>%
  count(hms_id, gene_symbol, cmpd_conc_nM) %>%
  count(n)
# # A tibble: 19 x 2
# n    nn
# <int> <int>
#   1     1 64128
# 2     2  2749
# 3     3   311
# 4     4   199
# 5     5    32
# 6     6    43
# 7     7    56
# 8     8    85
# 9     9   132
# 10    10   171
# 11    11    71
# 12    12   145
# 13    15   145
# 14    17     1
# 15    18     2
# 16    19     1
# 17    20     2
# 18    24     2
# 19    30     2
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
    hms_lincs_chembl_mapping %>%
      distinct(hms_id, chembl_id),
    by = "hms_id"
  ) %>%
  rename(gene_id = entrez_id, pref_name_cmpd = pref_name, pref_name_target = name) %>%
  drop_na(chembl_id, gene_id)

hmsl_kinomescan_q1 <- hmsl_kinomescan_mapped %>%
  as.data.table() %>%
  .[
    ,
    .(percent_control_Q1 = quantile(percent_control, 0.25, names = FALSE)),
    by = .(chembl_id, gene_id, cmpd_conc_nM)
    ] %>%
  as_tibble()

write_csv(hmsl_kinomescan_q1, "hmsl_kinomescan_q1_20190619.csv.gz")

