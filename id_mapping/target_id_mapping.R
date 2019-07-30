library(tidyverse)
library(vroom)
library(data.table)
library(biomaRt)

dir_chembl <- file.path("~", "repo", "tas_vectors", "chembl24_1")

##################################################################################################################T
# create target conversion table  ------------
##################################################################################################################T
dir.create(dir_chembl)
setwd(dir_chembl)
list.files()
download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_24_1/chembl_uniprot_mapping.txt",
  "chembl_uniprot_mapping.txt"
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
  "chembl_uniprot_mapping.txt", skip = 1,
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
# organism                n
# <chr>               <int>
# 1 Homo sapiens         18
# 2 Mus musculus        147
# 3 Rattus norvegicus   604

# Only 18 human uniprot IDs left without a matching Entrez ID, good enough...

map_uniprot_chembl %>%
  drop_na(gene_id) %>%
  dplyr::count(gene_id) %>%
  dplyr::count(n)
# # A tibble: 15 x 2
#       n    nn
#    <int> <int>
# 1     1  3860
# 2     2   658
# 3     3   192
# 4     4    97
# 5     5    35
# 6     6    23
# 7     7    16
# 8     8     3
# 9     9     5
# 10   10     5
# 11   11     3
# 12   12     1
# 13   13     1
# 14   14     1
# 15   16     1

# step 4 --> match with gene_info

## go to https://www.ncbi.nlm.nih.gov/gene
## click Download/FTP
## select gene_info table
## copy table to dir_chembl24_1

setwd(dir_chembl)
download.file("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz", "gene_info_20190725.gz")

#gene_info<-read.delim("gene_info_20190226")
#gene_info<-gene_info[c("X.tax_id","GeneID","Symbol","Synonyms","description","Full_name_from_nomenclature_authority")]
#dim(gene_info)
#write.csv(gene_info,file="gene_info_20190226_compact.csv",row.names = F)

# Using vroom here instead of loading the entire csv because it is downright massive
# and vroom is much faster
gene_info <- vroom(
  "gene_info_20190725.gz",
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

write_csv(map_chemblID_geneID, "map_targets_chemblID_genID_uniprot_tid.csv.gz")

