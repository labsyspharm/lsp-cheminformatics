library(tidyverse)
library(data.table)
library(biomaRt)
library(bit64)
library(here)
library(synapser)
library(synExtra)

source(here("id_mapping", "chemoinformatics_funcs.R"))

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

dir_unichem <- file.path(dir_release, "unichem")

dir.create(dir_unichem, showWarnings = FALSE)

# Download Unichem xref tables -------------------------------------------------
###############################################################################T

unichem_ftp <- "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/oracleDumps/UDRI238/"

download.file(
  file.path(unichem_ftp, "UC_SOURCE.txt.gz"),
  file.path(dir_unichem, "unichem_sources.txt.gz")
)

system2(
  "wget",
  c(
    "--recursive",
    "--no-parent",
    "-q",
    "-R", "'index.html*'",
    "-nH",
    "-P", dir_unichem,
    "--cut-dirs=6",
    "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1"
  )
)

unichem_sources <- read_tsv(
  file.path(dir_unichem, "unichem_sources.txt.gz"),
  col_names = c("src_id", "src_name", "date_created", "base_url"),
  col_types = "ic___c_____c______________"
) %>%
  filter(
    src_id %in% c(3L, 4L, 7L, 9L, 10L, 14L, 20L, 22L)
  ) %>%
  mutate(xref_type = paste0(src_name, "_id_compound"))

Xref_total <- unichem_sources %>%
  mutate(
    xref_map = map(
      src_id,
      function(src_id) {
        read_tsv(
          file.path(dir_unichem, "src_id1", paste0("src1src", src_id, ".txt.gz")),
          col_names = c("chembl_id_compound", "Xref_id_compound"),
          col_types = "cc",
          skip = 1
        )
      }
    )
  ) %>%
  unnest(xref_map) %>%
  mutate(url = paste0(base_url, Xref_id_compound))

write_csv(Xref_total, file.path(dir_unichem, "xref_table_unichem.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

unichem_wrangling_activity <- Activity(
  name = "Wrangle Unichem cross-references",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/05_unichem_xref_mapping.R"
)

list(
  file.path(dir_unichem, "xref_table_unichem.csv.gz")
) %>%
  map(
    . %>%
      File(parent = syn_release) %>%
      synStore(activity = unichem_wrangling_activity)
  )
