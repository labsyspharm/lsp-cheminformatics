library(tidyverse)
library(httr)
library(furrr)
library(RPostgres)
library(bit64)
library(jsonlite)
library(data.table)
library(here)
library(synapser)
library(synExtra)

source(here("id_mapping", "chemoinformatics_funcs.R"))

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


# Load compound data -----------------------------------------------------------
###############################################################################T

compound_sources <- tribble(
  ~source_name, ~synapse_id,
  "chembl", "syn20692439",
  "hmsl", "syn20692442"
)

compounds_canonical <- compound_sources %>%
  mutate(
    canonical = map(
      synapse_id,
      . %>%
        syn() %>%
        read_csv(col_types = "ccc")
    ) %>%
      map(
        rename_all,
        str_replace,
        "chembl_id|hms_id", "id"
      )
  ) %>%
  select(-synapse_id) %>%
  unnest(canonical)

write_csv(
  compounds_canonical,
  file.path(dir_release, "all_compounds_canonical.csv.gz")
)


# Create molecular fingerprints ------------------------------------------------
###############################################################################T

fingerprinting_args <- tribble(
  ~fp_name, ~fp_type, ~fp_args,
  "morgan_chiral", "morgan", c("--useChirality", "1"),
  "morgan_normal", "morgan", c("--useChirality", "0")
)

plan(multisession, workers = 10)
cmpd_fingerprints <- compounds_canonical %>%
  select(name = id, compound = inchi) %>%
  # Some HMSL compounds don't report inchi or smiles or anything, should have removed
  # them earlier, removing them here now
  drop_na() %>%
  slice(sample(nrow(.))) %>%
  # slice(1:100) %>%
  split((seq(nrow(.)) - 1) %/% 10000) %>%
  enframe("index", "compounds") %>%
  mutate(names = map(compounds, "name"), compounds = map(compounds, "compound")) %>%
  crossing(fingerprinting_args) %>%
  mutate(
    fp_res = future_pmap(
      select(., compounds, names, fp_type, fp_args),
      safely(get_fingerprints),
      .progress = TRUE
    )
  )

write_rds(cmpd_fingerprints, file.path(dir_release, "all_compounds_fingerprints_raw.rds"))
# cmpd_fingerprints <- read_rds(file.path(dir_release, "all_compounds_fingerprints_raw.rds"))

merge_fp_res <- function(fp_res) {
  fingerprint_db_lines <- map(fp_res, pluck, "result", "fingerprint_db") %>%
  # split into lines
  map(stringi::stri_split_lines1)

  fingerprint_fps_combined <- c(
    # Keep first file including header
    fingerprint_db_lines[[1]],
    # Remaining lines from all files without headers
    fingerprint_db_lines[-1] %>%
      # In all following files, discard  header (lines starting with #)
      map(~.x[!str_starts(.x, fixed("#"))]) %>%
      # merge all remaining lines
      {do.call(c, .)}
  )

  list(
    skipped = do.call(c, map(fp_res, pluck, "result", "skipped")),
    fingerprint_db = fingerprint_fps_combined,
    error = do.call(c, map(fp_res, pluck, "error"))
  )
}

cmpd_fingerprints_all <- cmpd_fingerprints %>%
  group_by(fp_name, fp_type) %>%
  summarize(
    result = list(merge_fp_res(fp_res))
  ) %>%
  ungroup() %>%
  unnest_wider(result)

# Check if any jobs errored out
cmpd_fingerprints %>%
  pull("error") %>%
  map(is.null) %>%
  as.logical() %>%
  all()

skipped_cmpds <- cmpd_fingerprints_all %>%
  unnest(skipped) %>%
  select(-fingerprint_db)

write_lines(skipped_cmpds, file.path(dir_release, "all_compounds_fingerprints_skipped.txt"))

write_rds(
  cmpd_fingerprints_all,
  file.path(dir_release, "all_compounds_fingerprints.rds")
)
cmpd_fingerprints_all %>%
  mutate(
    fn = file.path(
      dir_release,
      paste0(paste("all_compounds_fingerprints_fp_type", fp_type, "fp_name", fp_name, sep = "_"), ".fps")
    )
  ) %>%
  {walk2(.$fingerprint_db, .$fn, write_lines)}

# Store to synapse -------------------------------------------------------------
###############################################################################T

fingerprint_activity <- Activity(
  name = "Calculate molecular fingerprints",
  used = c(
    "syn20692439",
    "syn20692442"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_fingerprints.R"
)

synExtra::synPluck(syn_release, "fingerprints")
fp_folder <- Folder("fingerprints", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_fingerprints.rds"),
  Sys.glob(file.path(dir_release, "all_compounds_fingerprints*fps")),
  file.path(dir_release, "all_compounds_canonical.csv.gz")
) %>%
  synStoreMany(parent = fp_folder, activity = fingerprint_activity)
