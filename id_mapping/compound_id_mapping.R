library(tidyverse)
library(httr)
library(RPostgres)
library(bit64)

dir_chembl <- file.path("~", "repo", "tas_vectors", "chembl24_1")

## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_24",
                 host = "localhost", port = 5433,
                 user = "chembl_public")



##################################################################################################################T
# get table w/ all molecule structures, indicate if cmpd has different parental, create temp_id  ------------
##################################################################################################################T

# Chembl has salts sometimes separate from free bases, parent_molregno points to
# free base
# Trying to reconcile internal LINCS database, which only has a single ID for
# each compound, regardless of salt, with Chembl

setwd(dir_chembl)

# step 1--> get basic info on molecules
all_cmpds <- dbGetQuery(
  con,
  paste0(
    "select distinct dict.molregno, dict.pref_name, dict.chembl_id, dict.max_phase, syn.synonyms,
    hi.parent_molregno, hi.active_molregno, struct.standard_inchi, struct.standard_inchi_key
   from molecule_dictionary as dict
   left join molecule_synonyms as syn on dict.molregno = syn.molregno
   left join molecule_hierarchy as hi on dict.molregno = hi.molregno
   left join compound_structures as struct on dict.molregno=struct.molregno"
  )
) %>%
  as_tibble() %>%
  mutate(parental_flag = ifelse(molregno != parent_molregno, 1L, 0L))

write_csv(all_cmpds, "all_compounds_chembl24_1.csv.gz")
# Again, molregno, parent_molregno and active_molregno are integer64...


# Trying to map compounds with an annotated parent compound to its parents
# ChemblID
all_cmpds %>% filter(parental_flag==1) %>% dim()
all_cmpds %>% filter(parental_flag==0 | is.na(parental_flag)) %>% dim()

all_cmpds_parent_mapped <- all_cmpds %>%
  left_join(
    all_cmpds %>%
      distinct(parent_chembl_id = chembl_id, molregno) %>%
      drop_na(),
    by = c("parent_molregno" = "molregno")
  ) %>%
  mutate(
    parent_chembl_id = ifelse(
      parental_flag == 0 | is.na(parental_flag),
      chembl_id,
      parent_chembl_id
    )
  )

all_cmpds_parent_mapped %>% filter(is.na(parent_chembl_id)) %>% dim()
write_csv(all_cmpds_parent_mapped, "all_compounds_chembl24_1_parent_mapped.csv.gz")

#################################################################################T
# Mapping LINCS IDs to Chembl -----------------------------
#################################################################################T

# Attempt to get hms lincs ID map automagically from HMS LINCS website
hms_lincs_compounds_json <- GET(
  modify_url(
    "http://lincs.hms.harvard.edu",
    path = "/db/api/v1/smallmolecule/",
    query = list("limit" = "100000")
  ),
  accept_json()
)

hms_lincs_compounds <- content(hms_lincs_compounds_json, "text") %>%
  jsonlite::fromJSON() %>%
  .$objects %>%
  as_tibble() %>%
  mutate(hms_id = paste0("HMSL", str_match(hmsLincsID, "([0-9]+)-([0-9]+)")[, 2]))

# Disregard the annotated Chembl ID in the HMS LINCS database because it may
# point to the salt compound and we always want the free base ID
# Only use annotated LINCS Chembl ID if we can't find it using the inchi_key

# We have to first generate a list of possible tautomers of the molecules,
# the inchi listed is not necessarily the one used in Chembl

request_tautomers <- function(inchi, max_tautomers = 10) {
  res <- POST(
    "http://127.0.0.1:5000/query/tautomers",
    body = list(
      "inchi" = inchi
    ),
    encode = "multipart",
    accept_json()
  )
  content(res, as = "text") %>%
    jsonlite::fromJSON() %>%
    .$tautomers %>%
    as_tibble()
}

hms_lincs_compounds_tautomers <- hms_lincs_compounds %>%
  distinct(hms_id, smInchiParent) %>%
  filter(smInchiParent != "-") %>%
  mutate(
    tautomers = map(
      smInchiParent,
      request_tautomers
    )
  )

hms_lincs_compounds_tautomers_flat <- bind_rows(
  # First the original inchis from HMS LIINCS
  hms_lincs_compounds %>%
    distinct(hms_id, inchi = smInchiParent) %>%
    mutate(inchi_source = "HMS_LINCS"),
  # And the found tautomers
  hms_lincs_compounds_tautomers %>%
    mutate(
      tautomers = tautomers %>%
        map(mutate, inchi_source = "HMS_LINCS_TAUTOMERS") %>%
        map(distinct, inchi_source, inchi) %>%
        # Remove cases where the found tautomer is exactly the same as the query
        map2(smInchiParent, ~filter(.x, inchi != .y))
    ) %>%
    select(hms_id, tautomers) %>%
    unnest() %>%
    distinct()
) %>%
  arrange(hms_id)

hms_lincs_compounds_mapped <- hms_lincs_compounds_tautomers_flat %>%
  left_join(
    all_cmpds_parent_mapped %>%
      # Using parent_chembl_id because it points to the parent free base ID
      distinct(chembl_id = parent_chembl_id, standard_inchi),
    by = c("inchi" = "standard_inchi")
  )

# Any cases where for a single hms_id multiple Chembl IDs where found?
hms_lincs_compounds_mapped %>%
  drop_na(chembl_id) %>%
  count(hms_id) %>%
  count(n)
# # A tibble: 1 x 2
# n    nn
# <int> <int>
#   1     1   465

hms_lincs_compounds_mapped_merged <- bind_rows(
  # Coumpounds where we got the chembl_id by querying chembl with the inchi
  hms_lincs_compounds_mapped %>%
    mutate(chembl_id_source = "INCHI_MAPPING"),
  # Compounds where a Chembl ID was stored in HMS LINCS
  hms_lincs_compounds %>%
    distinct(hms_id, inchi = smInchiParent, chembl_id = chemblID) %>%
    drop_na(chembl_id) %>%
    mutate(
      chembl_id_source = "HMS_LINCS",
      chembl_id = paste0("CHEMBL", chembl_id),
      inchi_source = "HMS_LINCS"
    )
) %>%
  arrange(hms_id) %>%
  drop_na(chembl_id) %>%
  group_by(hms_id, inchi, chembl_id) %>%
  summarize(
    inchi_source = paste(unique(inchi_source), collapse = "/"),
    chembl_id_source = paste(unique(chembl_id_source), collapse = "/")
  ) %>%
  ungroup()

hms_lincs_compounds_mapped_merged %>%
  count(hms_id, chembl_id) %>%
  count(n)

write_csv(hms_lincs_compounds_mapped_merged, "hms_lincs_chembl_mapping.csv.gz")
