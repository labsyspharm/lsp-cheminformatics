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

allowed_target_types = c(
  "SELECTIVITY GROUP",
  "PROTEIN NUCLEIC-ACID COMPLEX",
  "PROTEIN FAMILY",
  "CHIMERIC PROTEIN",
  "PROTEIN COMPLEX",
  "SINGLE PROTEIN",
  "PROTEIN COMPLEX GROUP"
)

chembl_info_all_targets <- dbGetQuery(
  con,
  paste0(
    "SELECT dict.tid, dict.pref_name, dict.tax_id, dict.organism, dict.chembl_id, dict.target_type
    FROM target_dictionary AS dict"
  )
) %>%
  as_tibble() %>%
  filter(organism %in% names(organism_biomart_mapping), target_type %in% allowed_target_types)

chembl_info_all_targets %>%
  count(chembl_id, target_type) %>%
  count(n)
# # A tibble: 1 x 2
# n    nn
# <int> <int>
#   1     1  5842
# Only a single target type per chembl_id


map_uniprot_chembl <- read_tsv(
  file.path(dir_release, "chembl_uniprot_mapping.txt"), skip = 1,
  col_names = c("uniprot_id", "chembl_id", "pref_name", "target_type")
)


map_with_uniprot <- chembl_info_all_targets %>%
  left_join(
    map_uniprot_chembl %>%
      dplyr::select(chembl_id, uniprot_id),
    by = "chembl_id"
  )

download.file(
  "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz",
  file.path(dir_release, "uniprot_id_mapping.tsv.gz")
)

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

map_with_entrez <- map_with_uniprot %>%
  group_by(organism) %>%
  group_modify(uniprot_to_entrez) %>%
  dplyr::rename(entrez_gene_id = entrezgene_id) %>%
  ungroup()

map_with_entrez %>%
  filter(is.na(entrez_gene_id)) %>%
  dplyr::count(organism)
# # A tibble: 3 x 2
# organism              n
# <chr>             <int>
#   1 Homo sapiens         15
# 2 Mus musculus        149
# 3 Rattus norvegicus   620

# Only 15 human uniprot IDs left without a matching Entrez ID, good enough...

map_with_entrez %>%
  drop_na(entrez_gene_id) %>%
  dplyr::count(entrez_gene_id) %>%
  dplyr::count(n)
# # A tibble: 15 x 2
# n    nn
# <int> <int>
#   1     1  4078
# 2     2   653
# 3     3   200
# 4     4    93
# 5     5    36
# 6     6    20
# 7     7    16
# 8     8     3
# 9     9     5
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
    "tax_id", "entrez_gene_id", "entrez_symbol", "locus_tag", "entrez_synonyms", "db_xrefs", "chromosome",
    "map_location", "entrez_description", "entrez_type_of_gene", "symbol", "entrez_name", "nomenclature_status",
    "other_designations", "modification_date", "feature_type"
  ),
  col_types = "iic_c___cc_c____",
  skip = 1
)

map_chemblID_geneID <- map_with_entrez %>%
  # Have to cast as integer64 because the Chembl postgresql DB stores some ids
  # (tid, tax_id) as 64bit integer
  left_join(
    mutate_at(gene_info, "tax_id", bit64::as.integer64),
    by = c("tax_id", "entrez_gene_id")
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
