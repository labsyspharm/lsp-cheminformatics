library(tidyverse)
library(httr)
library(furrr)
library(RPostgres)
library(bit64)
library(jsonlite)
library(data.table)

source("chemoinformatics_funcs.R")
# source("../id_mapping/chemoinformatics_funcs.R")

dir_chembl <- file.path("~", "repo", "tas_vectors", "chembl25")

## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_25",
                 host = "localhost", port = 5432,
                 user = "chug")

# Fetch Chembl compound data ---------------------------------------------------
###############################################################################T

# Chembl has salts sometimes separate from free bases, parent_molregno points to
# free base
# Trying to reconcile internal LINCS database, which only has a single ID for
# each compound, regardless of salt, with Chembl

dir.create(dir_chembl, showWarnings = FALSE)
setwd(dir_chembl)

# step 1--> get basic info on molecules
all_cmpds <- dbGetQuery(
  con,
  "SELECT DISTINCT dict.molregno, dict.pref_name, dict.chembl_id, dict.max_phase,
    hi.parent_molregno, hi.active_molregno, struct.standard_inchi, struct.standard_inchi_key,
    struct.canonical_smiles, dict.inorganic_flag, ass.n_assays
  FROM molecule_dictionary AS dict
  LEFT JOIN molecule_hierarchy AS hi ON dict.molregno = hi.molregno
  LEFT JOIN compound_structures AS struct ON dict.molregno=struct.molregno
  LEFT JOIN (
    SELECT COUNT(activity_id) AS n_assays, molregno FROM activities GROUP BY molregno
  ) AS ass ON dict.molregno = ass.molregno"
) %>%
  as_tibble() %>%
  mutate(parental_flag = ifelse(molregno != parent_molregno, 1L, 0L))

all_synonyms <- dbGetQuery(
  con,
  "SELECT DISTINCT dict.molregno, syn.synonyms
  FROM molecule_dictionary AS dict
  LEFT JOIN molecule_synonyms AS syn ON DICT.molregno = syn.molregno"
) %>%
  drop_na() %>%
  as.data.table() %>%
  .[
    ,
    .(synonyms = list(synonyms)),
    keyby = molregno
  ]

all_cmpds <- all_cmpds %>%
  left_join(
    all_synonyms,
    by = "molregno"
  )

write_rds(
  all_cmpds,
  "all_compounds_chembl25.rds",
  compress = "gz"
)

all_cmpds <- all_cmpds %>%
  mutate(
    synonyms = map_chr(synonyms, paste, collapse = "; ")
  )

write_csv(
  all_cmpds,
  "all_compounds_chembl25.csv.gz"
)
# all_cmpds <- read_csv(
#   "all_compounds_chembl25.csv.gz",
#   col_types = "icciiiccciiic"
# )
# Again, molregno, parent_molregno and active_molregno are integer64...

unique_cmpds <- all_cmpds %>%
  distinct(compound = standard_inchi, name = chembl_id) %>%
  mutate(compound = trimws(compound)) %>%
  drop_na()

# Find canonical tautomer ------------------------------------------------------
###############################################################################T

library(ssh)
library(uuid)

session <- ssh_connect("ch305@o2.hms.harvard.edu", keyfile = "~/.ssh/id_ecdsa")
remote_wd <- file.path("/n", "scratch2", "ch305", paste0("chembl_25_canonicalize_", UUIDgenerate()))
# remote_wd <- file.path("/n", "scratch2", "ch305", "chembl_25_canonicalize_a21b056d-a655-486f-af70-15ad2b09db8c")
ssh_exec_wait(session, paste0("mkdir -p ", remote_wd))

local_temp_wd <- tempdir()

set.seed(42)
cmpd_split <- unique_cmpds %>%
  slice(sample(nrow(.))) %>%
  split((seq(nrow(.)) - 1) %/% 10000) %>%
  enframe("index", "compounds") %>%
  mutate(
    path = file.path(
      local_temp_wd, paste0("chembl_compounds_", sprintf("%03d", as.integer(index)), ".csv")
    )
  )

walk2(
  cmpd_split$compounds, cmpd_split$path,
  write_csv
)
scp_upload(session, cmpd_split$path, remote_wd)

# Processing 1000 compounds takes around 8 min on the cluster, so for 10000
# compounds in every file scheduling a 2 h job
ssh_exec_cmd <- paste(
  "sbatch", "-c", 1, "--time=05:00:00", "-p short", "--mem=6000",
  "-D", remote_wd,
  "--wrap",
  "'",
  paste(
    "source", "~/miniconda3/etc/profile.d/conda.sh", ";",
    "conda", "activate", "rdkit", ";",
    "tautomers", "canonicalize",
    "--compound-col", "compound",
    "--compound-encoding", "inchi",
    "--standardize",
    file.path(remote_wd, basename(cmpd_split$path)),
    file.path(remote_wd, paste0(tools::file_path_sans_ext(basename(cmpd_split$path)), "_canonical.csv"))
  ),
  "'"
)
ssh_exec_cmd %>% paste(collapse = "\n") %>% clipr::write_clip()

# submitting the jobs using ssh_exec_wait results in segfaul for whatever reason...
# ssh_exec_wait(session, ssh_exec_cmd)

# Checking which jobs completed sucessfully
job_status <- ssh_exec_internal(session, paste("ls", "-lh", remote_wd)) %>%
  pluck("stdout") %>%
  rawToChar() %>%
  read_delim(delim = " ", trim_ws = TRUE, skip = 1, col_names = FALSE)

sucessful <- job_status %>%
  filter(str_detect(X9, "canonical.csv")) %>%
  mutate(index = as.integer(str_extract(X9, "[0-9]+")))

cmpd_split_repeat <- cmpd_split %>%
  filter(!(index %in% sucessful$index))

ssh_exec_cmd_repeat <- paste(
  "sbatch", "-c", 1, "--time=4:00:00", "-p short", "--mem=20000",
  "-D", remote_wd,
  "--wrap",
  "'",
  paste(
    "source", "~/miniconda3/etc/profile.d/conda.sh", ";",
    "conda", "activate", "rdkit", ";",
    "tautomers", "canonicalize",
    "--compound-col", "compound",
    "--compound-encoding", "inchi",
    file.path(remote_wd, basename(cmpd_split_repeat$path)),
    file.path(remote_wd, paste0(tools::file_path_sans_ext(basename(cmpd_split_repeat$path)), "_canonical.csv"))
  ),
  "'"
)

ssh_exec_cmd_repeat %>% paste(collapse = "\n") %>% clipr::write_clip()

canonical_path <- file.path("chembl_compounds", "canonical") %>%
  normalizePath()
dir.create(canonical_path, showWarnings = FALSE, recursive = TRUE)

scp_download(session, file.path(remote_wd, "*_canonical.csv"), to = canonical_path)

canonical <- list.files(
    canonical_path,
    pattern = "chembl_compounds_[0-9]+_canonical\\.csv",
    full.names = TRUE
  ) %>%
  map(read_csv, col_types = "ccc") %>%
  bind_rows() %>%
  left_join(
    all_cmpds %>%
      select(chembl_id, chembl_inchi = standard_inchi),
    by = c("query" = "chembl_inchi")
  )

canonical_all <- canonical %>%
  distinct() %>%
  left_join(
    all_cmpds %>%
      select(chembl_id, inorganic_flag),
    by = "chembl_id"
  ) %>%
  mutate(
    inchi = case_when(
      # Whe canonicalization failed, use original inchi
      is.na(inchi) ~ query,
      # Using non-canonicalized inchi for compounds that are classified as non-organic
      inorganic_flag == 1 ~ query,
      TRUE ~ inchi
    )
  ) %>%
  select(
    chembl_id,
    original_inchi = query,
    inchi
  )

write_csv(canonical_all, file.path(canonical_path, "all_canonical.csv.gz"))
# canonical_all <- read_csv(file.path(canonical_path, "all_canonical.csv.gz"))

# Create molecular fingerprints ------------------------------------------------
###############################################################################T

plan(multisession, workers = 8)
cmpd_fingerprints <- canonical_all %>%
  select(name = chembl_id, compound = inchi) %>%
  split((seq(nrow(.)) - 1) %/% 10000) %>%
  enframe("index", "compounds") %>%
  mutate(
    fingerprint_res = future_map(
      compounds, get_fingerprints, .progress = TRUE
    )
  )
write_rds(cmpd_fingerprints, "chembl_fingerprints.rds")
# cmpd_fingerprints <- read_rds("chembl_fingerprints.rds")

fingerprint_fps <- cmpd_fingerprints$fingerprint_res %>%
  map("fps_file") %>%
  map(rawToChar) %>%
  # split into lines
  map(stringi::stri_split_lines1)

# Have to remove header of all fps files except first one

fingerprint_fps_combined <- c(
  fingerprint_fps[[1]],
  fingerprint_fps[-1] %>%
    # discard lines starting with #
    map(~.x[!str_starts(.x, fixed("#"))]) %>%
    {do.call(c, .)}
)

writeLines(fingerprint_fps_combined, "all_compounds_fingerprints.fps")

# Calculate similarity between Chembl compounds --------------------------------
###############################################################################T


cmdp_similarity <- cmpd_fingerprints %>%
  mutate(
    similarity_res = future_map(
      fingerprint_res,
      function(fp) {
        fp_out <- tempfile(fileext = ".fps")
        writeBin(fp$fps_file, fp_out)
        scan_fingerprint_matches(
          "all_compounds_fingerprints.fps",
          query = fp_out,
          threshold = 0.99
        )
      },
      .progress = TRUE
    )
  )

write_rds(cmdp_similarity, "all_compounds_similarity_raw.rds")
# cmdp_similarity <- read_rds("all_compounds_similarity_raw.rds")

similarity_df <- cmdp_similarity$similarity_res %>%
  map(as_tibble) %>%
  bind_rows() %>%
  mutate(score = as.double(score))

write_csv(similarity_df, "all_compounds_similarity.csv.gz")

cmpds_identical <- similarity_df %>%
  filter(match != query, score > 0.999999) %>%
  select(-score)




# Defining equivalence classes -------------------------------------------------
###############################################################################T


library(igraph)

# Making graph of compounds where edges represent identical compounds
# Extracting the "components" of the graph, isolated subgraphs, they
# represent equivalence clusters of identical compounds

cmpds_identical_graph <- cmpds_identical %>%
  mutate(
    chembl_id1 = ifelse(match > query, query, match),
    chembl_id2 = ifelse(match > query, match, query)
  ) %>%
  distinct(chembl_id1, chembl_id2) %>%
  as.matrix() %>%
  graph_from_edgelist(directed = FALSE)

cmpd_eq_classes <- components(cmpds_identical_graph) %>%
  pluck("membership") %>%
  enframe(value = "eq_class", name = "chembl_id") %>%
  mutate(eq_class = as.integer(eq_class))

all_cmpds_with_parent <- all_cmpds %>%
  filter(molregno != parent_molregno) %>%
  left_join(
    all_cmpds %>%
      select(molregno, parent_chembl_id = chembl_id, parent_standard_inchi = standard_inchi),
    by = c("parent_molregno" = "molregno")
  ) %>%
  select(molregno, chembl_id, parent_chembl_id, standard_inchi, parent_standard_inchi) %>%
  left_join(
    cmpd_eq_classes,
    by = "chembl_id"
  ) %>%
  left_join(
    cmpd_eq_classes %>%
      rename(parent_eq_class = eq_class),
    by = "chembl_id"
  )

all_cmpds_with_parent %>%
  filter(eq_class != parent_eq_class)
# Nothing, so we found all of the parent molecule relationships using our strategy
# of canonicalization

cmpd_eq_classes_solo <- unique_cmpds %>%
  distinct(chembl_id = name) %>%
  filter(!(chembl_id %in% cmpd_eq_classes$chembl_id)) %>%
  mutate(
    eq_class = seq(max(cmpd_eq_classes$eq_class) + 1, max(cmpd_eq_classes$eq_class) + n())
  )

cmpd_eq_classes_all <- bind_rows(
    cmpd_eq_classes,
    cmpd_eq_classes_solo
) %>%
  left_join(
    all_cmpds %>%
      distinct(molregno, chembl_id, chembl_inchi = standard_inchi, pref_name, max_phase, parent_molregno, active_molregno, n_assays),
    by = "chembl_id"
  ) %>%
  left_join(
    canonical_all %>%
      distinct(chembl_id, canonical_inchi = inchi),
    by = "chembl_id"
  )
write_csv(
  cmpd_eq_classes_all,
  "equivalence_classes.csv.gz"
)
# cmpd_eq_classes_all <- read_csv("equivalence_classes.csv.gz", col_types = "ciicciiiic")

fingerprints_tbl <- fingerprint_fps_combined %>%
  read_tsv(col_names = c("fingerprint", "chembl_id"), skip = 6)

fingerprints_eq_classes <- cmpd_eq_classes_all %>%
  select(chembl_id, eq_class) %>%
  inner_join(
    fingerprints_tbl,
    by = "chembl_id"
  ) %>%
  distinct(fingerprint, eq_class)

# FPS Header
write_lines(fingerprint_fps_combined[1:6], "eq_classes_fingerprints.fps")
# Content
write_tsv(
  fingerprints_eq_classes,
  "eq_classes_fingerprints.fps",
  col_names = FALSE,
  append = TRUE,
  quote_escape = "none"
)

# Finding canonical member of equivalence class --------------------------------
###############################################################################T

# Find canonical member of equivalence class
# 1. If present, choose member annotated as parent by Chembl
# 2. Choose member with highest clinical phase
# 3. Choose member with annotated pref_name
# 4. Choose member with highest n_assays
# 5. Choose member with lowest chembl_id
cmpd_eq_classes_canonical <- cmpd_eq_classes_all %>%
  mutate(
    parent_present = if_else(parent_molregno == molregno, NA_integer_, as.integer(parent_molregno)),
    molregno = as.integer(molregno),
    n_assays = as.integer(n_assays)
  ) %>%
  as.data.table() %>%
  .[
    ,
    `:=`(
      annotated_as_parent = as.integer(molregno) %in% parent_present,
      annotated_pref_name = !is.na(pref_name),
      annotated_assays = replace_na(n_assays, 0)
    ),
    keyby = eq_class
  ] %>%
  .[
    order(
      eq_class,
      -annotated_as_parent,
      -max_phase,
      -annotated_pref_name,
      -annotated_assays,
      -chembl_id
    ),
    head(.SD, 1),
    keyby = eq_class
  ]
write_csv(
  cmpd_eq_classes_canonical,
  "canonical_equivalence_classes.csv.gz"
)


