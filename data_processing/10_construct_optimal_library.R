library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Load files -------------------------------------------------------------------
###############################################################################T

selectivity <- syn("syn20836653") %>%
  read_rds()

canonical_fps <- syn("syn21042105") %>%
  read_rds()

# Prepare chemical similarity  files -------------------------------------------
###############################################################################T

fps_headers <- list(
morgan_normal = "#FPS1
#num_bits=2048
#type=RDKit-Morgan/1 radius=2 fpSize=2048 useFeatures=0 useChirality=0 useBondTypes=1
#software=RDKit/2018.09.3 chemfp/1.5
",
morgan_chiral = "#FPS1
#num_bits=2048
#type=RDKit-Morgan/1 radius=2 fpSize=2048 useFeatures=0 useChirality=1 useBondTypes=1
#software=RDKit/2018.09.3 chemfp/1.5
",
topological_normal = "#FPS1
#num_bits=2048
#type=RDKit-Fingerprint/2 minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=2 useHs=1
#software=RDKit/2018.09.3 chemfp/1.5
"
)

fps_selective <- canonical_fps %>%
  left_join(rename(selectivity, selectivity_df = data)) %>%
  mutate(
    data = pmap(
      list(data, fp_name, selectivity_df),
      ~filter(
        ..1,
        fp_name == ..2,
        lspci_id %in% filter(..3, selectivity_class %in% c("most_selective", "semi_selective"))$lspci_id
      ) %>%
        select(fingerprint, lspci_id)
    ),
    selectivity_fn = file.path(dir_release, paste0("fps_selective_", fp_name, ".fps"))
  )

pwalk(
  fps_selective,
  function(data, selectivity_fn, fp_name, ...) {
    write_file(
      fps_headers[[fp_name]],
      selectivity_fn
    )
    write_tsv(
      data,
      selectivity_fn,
      append = TRUE
    )
  }
)


source(here("id_mapping", "chemoinformatics_funcs.R"))
chemical_sim_selective <- fps_selective %>%
  mutate(
    data = map(
      selectivity_fn,
      scan_fingerprint_matches, threshold = 0.01
    ) %>%
      # Use only upper triangle of symmetric matrix
      map(mutate_at, vars("match", "query"), as.integer) %>%
      map(mutate_at, vars("score"), as.double) %>%
      map(filter, query < match) %>%
      map(select, lspci_id_1 = query, lspci_id_2 = match, tanimoto_similarity = score)
  )

# Construct optimal library ----------------------------------------------------
###############################################################################T

find_optimal_compounds <- function(
  selectivity, chemical_similarity
) {
  message(".")
  chemical_similarity %>%
    filter(
      lspci_id_1 %in% selectivity$lspci_id,
      lspci_id_2 %in% selectivity$lspci_id
    ) %>%
    left_join(
      select(
        selectivity,
        lspci_id,
        tool_score,
        selectivity_class
      ) %>%
        rename_all(paste0, "_1"),
      by = "lspci_id_1"
    ) %>%
    left_join(
      select(
        selectivity,
        lspci_id,
        tool_score,
        selectivity_class
      ) %>%
        rename_all(paste0, "_2"),
      by = "lspci_id_2"
    ) %>%
    mutate(tool_score_sum = tool_score_1 + tool_score_2) %>%
    arrange(
      as.integer(selectivity_class_1) + as.integer(selectivity_class_2),
      desc(tool_score_sum),
      tanimoto_similarity
    )
}

find_optimal_compounds_all_targets <- function(
  selectivity, chemical_sim, maximum_tanimoto_similarity = 0.2,
  selectivity_classes = c("most_selective", "semi_selective", "poly_selective")
) {
  sel_by_gene <- selectivity %>%
    filter(selectivity_class %in% selectivity_classes) %>%
    mutate(selectivity_class = factor(selectivity_class, levels = selectivity_classes)) %>%
    group_nest(gene_id)
  chemical_sim_filtered <- chemical_sim %>%
    filter(tanimoto_similarity <= maximum_tanimoto_similarity)
  cmpds <- sel_by_gene %>%
    transmute(
      gene_id,
      all_pairs = map(
        data,
        find_optimal_compounds, chemical_sim_filtered
      ),
      best_pair = map(
        all_pairs,
        slice, 1
      )
    ) %>%
    unnest(best_pair)
  cmpds
}

optimal_libraries <- selectivity %>%
  rename(selectivity_df = data) %>%
  left_join(
    select(
      chemical_sim_selective,
      fp_type, fp_name, chemical_similarity_df = data
    ),
    by = c("fp_type", "fp_name")
  ) %>%
  mutate(
    data = pmap(
      list(selectivity_df, chemical_similarity_df),
      find_optimal_compounds_all_targets
    )
  )

write_rds(
  optimal_libraries %>%
    select(fp_type, fp_name, data),
  file.path(dir_release, "optimal_library.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

library_activity <- Activity(
  name = "Construct optimal compound library",
  used = c(
    "syn20836653",
    "syn21042105"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/10_construct_optimal_library.R"
)

syn_library_folder <- Folder("compound_library", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "optimal_library.rds")
) %>%
  synStoreMany(parentId = syn_library_folder, activity = library_activity)

