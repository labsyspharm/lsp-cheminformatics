library(tidyverse)
library(vroom)
library(data.table)
library(biomaRt)
library(bit64)
library(here)
library(synapser)
library(synExtra)
library(RPostgres)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# connect to chembl ------------------------------------------------------------
###############################################################################T


## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_25",
                 host = "localhost", port = 5432,
                 user = "chug")

# create target conversion table -----------------------------------------------
###############################################################################T

download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_25/chembl_uniprot_mapping.txt",
  file.path(dir_release, "chembl_uniprot_mapping.txt")
)

# step 1 --> import chembl generated file & convert to gene_ID

organism_biomart_mapping <- c(
  "Homo sapiens" = "hsapiens_gene_ensembl",
  "Rattus norvegicus" = "rnorvegicus_gene_ensembl",
  "Mus musculus" = "mmusculus_gene_ensembl"
)

chembl_info_all_targets <- dbGetQuery(
  con,
  paste0(
    "select dict.tid, dict.pref_name, dict.tax_id, dict.organism, dict.chembl_id
    from target_dictionary as dict"
  )
) %>%
  as_tibble() %>%
  filter(organism %in% names(organism_biomart_mapping))

map_uniprot_chembl <- read_tsv(
  file.path(dir_release, "chembl_uniprot_mapping.txt"), skip = 1,
  col_names = c("uniprot_id", "chembl_id", "protein_name", "target_type")
) %>%
  inner_join(chembl_info_all_targets, by = "chembl_id")

uniprot_to_entrez <- function(df, group) {
  mart <- useMart(
    "ENSEMBL_MART_ENSEMBL",
    organism_biomart_mapping[group$organism]
  )
  map_uniprot_geneID <- biomaRt::select(
    mart, unique(df$uniprot_id),
    c("uniprot_gn_id", "entrezgene_id"), "uniprot_gn_id"
  )
  df <- df %>%
    left_join(map_uniprot_geneID, by = c("uniprot_id" = "uniprot_gn_id"))
  if (group$organism == "Homo sapiens") {
    # https://github.com/datarail/genebabel
    # Additionally query the genebabel package for human genes
    df <- df %>%
      filter(is.na(entrezgene_id)) %>%
      dplyr::select(-entrezgene_id) %>%
      genebabel::join_hgnc(
        query_col = "uniprot_id",
        match_cols = "uniprot_ids",
        select_cols = "entrez_id"
      ) %>%
      mutate(entrez_id = as.integer(entrez_id)) %>%
      dplyr::rename(entrezgene_id = entrez_id) %>%
      bind_rows(filter(df, !is.na(entrezgene_id)))
  }
  df
}

map_uniprot_chembl <- map_uniprot_chembl %>%
  group_by(organism) %>%
  group_modify(uniprot_to_entrez) %>%
  dplyr::rename(gene_id = entrezgene_id) %>%
  ungroup()

map_uniprot_chembl %>%
  filter(is.na(gene_id)) %>%
  dplyr::count(organism)
# # A tibble: 3 x 2
# organism              n
# <chr>             <int>
#   1 Homo sapiens         19
# 2 Mus musculus        150
# 3 Rattus norvegicus   621

# Only 19 human uniprot IDs left without a matching Entrez ID, good enough...

map_uniprot_chembl %>%
  drop_na(gene_id) %>%
  dplyr::count(gene_id) %>%
  dplyr::count(n)
# # A tibble: 15 x 2
# n    nn
# <int> <int>
#   1     1  4051
# 2     2   685
# 3     3   211
# 4     4   100
# 5     5    41
# 6     6    22
# 7     7    17
# 8     8     4
# 9     9     6
# 10    10     7
# 11    11     2
# 12    12     3
# 13    13     1
# 14    16     1
# 15    17     1

# step 4 --> match with gene_info

download.file(
  "ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz",
  file.path(dir_release, "gene_info_20190829.gz")
)


# Using vroom here instead of loading the entire csv because it is downright massive
# and vroom is much faster
gene_info <- vroom(
  file.path(dir_release, "gene_info_20190829.gz"),
  delim = "\t",
  col_names = c(
    "tax_id", "gene_id", "ncbi_symbol", "locus_tag", "synonyms", "db_xrefs", "chromosome",
    "map_location", "description", "type_of_gene", "symbol", "name_ncbi", "nomenclature_status",
    "other_designations", "modification_date", "feature_type"
  ),
  col_types = "iic_c___c_cc____",
  skip = 1
)

map_chemblID_geneID <- map_uniprot_chembl %>%
  # Have to cast as integer64 because the Chembl postgresql DB stores some ids
  # (tid, tax_id) as 64bit integer
  left_join(
    mutate_at(gene_info, "tax_id", bit64::as.integer64),
    by = c("tax_id", "gene_id")
  ) %>%
  distinct()

write_csv(map_chemblID_geneID, file.path(dir_release, "target_map.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

target_wrangling_activity <- Activity(
  name = "Map and wrangle drug target data",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/06_target_id_mapping.R"
)

list(
  file.path(dir_release, "target_map.csv.gz")
) %>%
  map(
    . %>%
      File(parent = syn_release) %>%
      synStore(activity = target_wrangling_activity)
  )
