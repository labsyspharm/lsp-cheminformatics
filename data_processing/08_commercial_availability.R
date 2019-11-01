library(tidyverse)
library(synapser)
library(synExtra)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

compound_mapping <- syn("syn20934414") %>%
  read_rds()

# Fetch vendor info from ZINC --------------------------------------------------
###############################################################################T

vendor_libraries <- tribble(
  ~id, ~url, ~vendor_url,
  "targetmol", "http://files.docking.org/catalogs/40/", "https://www.targetmol.com/search?keyword=",
  "mce", "http://files.docking.org/catalogs/50/", "https://www.medchemexpress.com/search.html?q=",
  "selleck", "http://files.docking.org/catalogs/50/", "https://www.selleckchem.com/search.html?searchDTO.searchParam=",
  "tocris", "http://files.docking.org/catalogs/50/", "https://www.tocris.com/products/",
  "enamine", "http://files.docking.org/catalogs/50/", "https://www.enaminestore.com/catalog/",
  "enaminebb", "http://files.docking.org/catalogs/50/", "https://www.enaminestore.com/catalog/"
) %>%
  mutate(
    codemap_url = paste0(url, id, "/", id, ".codemap.txt.gz"),
    info_url = paste0(url, id, "/", id, ".info.txt.gz")
  )

pwalk(
  vendor_libraries,
  function(codemap_url, info_url, id, ...) {
    download.file(codemap_url, file.path(tempdir(), paste0("codemap_", id, ".txt.gz")))
    download.file(info_url, file.path(tempdir(), paste0("info_", id, ".txt.gz")))
  }
)

vendor_tables <- vendor_libraries %>%
  transmute(
    id,
    codemap = map(id, ~read_delim(file.path(tempdir(), paste0("codemap_", .x, ".txt.gz")), delim = " ", col_names = c("zinc_id", "vendor_id"), col_types = "cc")),
    info = map(id, ~read_tsv(file.path(tempdir(), paste0("info_", .x, ".txt.gz")), col_names = c("vendor_id", "zinc_id", "inchi_key", "tranche", "notes"), col_types = "ccccc")),
  )

download.file("http://files.docking.org/catalogs/1/chembl23/chembl23.codemap.txt.gz", file.path(tempdir(), "codemap_chembl23.txt.gz"))
chembl_codemap <- read_delim(file.path(tempdir(), "codemap_chembl23.txt.gz"), delim = " ", col_names = c("zinc_id", "chembl_id"), col_types = "cc")

vendor_tables_chembl <- vendor_tables %>%
  select(vendor = id, info) %>%
  unnest(info) %>%
  select(-tranche) %>%
  left_join(chembl_codemap, by = "zinc_id")

vendor_tables_lscpi <- compound_mapping %>%
  mutate(
    data = map(
      data,
      ~ left_join(vendor_tables_chembl, select(.x, id, lspci_id), by = c("chembl_id" = "id")) %>%
        select(-inchi_key, -chembl_id, -zinc_id, -notes) %>%
        distinct() %>%
        left_join(select(vendor_libraries, id, vendor_url), by = c("vendor" = "id")) %>%
        mutate(
          vendor_url = paste0(vendor_url, vendor_id)
        )
    )
  )

write_rds(
  vendor_tables_lscpi,
  file.path(dir_release, "compound_commercial_info_zinc.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

wrangle_activity <- Activity(
  name = "Wrangle commercial availability of compounds from ZINC",
  used = c(
    "syn20934414",
    "http://files.docking.org/catalogs/1/chembl23/chembl23.codemap.txt.gz",
    vendor_libraries$info_url
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/08_commercial_availability.R"
)

syn_vendors <- Folder("vendors", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "compound_commercial_info_zinc.rds")
) %>%
  synStoreMany(parentId = syn_vendors, activity = wrangle_activity)
