library(tidyverse)
library(vroom)
library(data.table)
library(biomaRt)
library(bit64)

dir_UniChem<- file.path("~", "repo", "tas_vectors", "unichem")

dir.create(dir_UniChem, showWarnings = FALSE)
setwd(dir_UniChem)
list.files()
system2(
  "wget",
  c(
    "--recursive",
    "--no-parent",
    "-q",
    "-R", "'index.html*'",
    "-nH",
    "--cut-dirs=6",
    "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1"
  )
)
download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/oracleDumps/UDRI232/UC_SOURCE.txt.gz",
  "UniChem_SourceCode.txt.gz"
)

unichem_sources <- read_tsv(
  "UniChem_SourceCode.txt.gz",
  col_names = c("source_number", "source_name", "date_created", "base_url"),
  col_types = "ic___c_____c______________"
) %>%
  filter(
    source_number %in% c(3L, 4L, 7L, 9L, 10L, 14L, 20L, 22L)
  ) %>%
  mutate(xref_type = paste0(source_name, "_id_compound"))

Xref_total <- unichem_sources %>%
  mutate(
    xref_map = map(
      source_number,
      function(source_number) {
        read_tsv(
          file.path(dir_UniChem, "src_id1", paste0("src1src", source_number, ".txt.gz")),
          col_names = c("chembl_id_compound", "Xref_id_compound"),
          col_types = "cc",
          skip = 1
        )
      }
    )
  ) %>%
  unnest(xref_map) %>%
  mutate(url = paste0(base_url, Xref_id_compound))

write_csv(Xref_total, "xref_table_unichem_20190325.csv.gz")
