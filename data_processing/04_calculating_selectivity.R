## this script calculates selectivity scores from biochemical data obtained with scripts '02_collecting_data_chembl.R'


library(tidyverse)
library(data.table)
library(here)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))



`%nin%` <- Negate(`%in%`)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


activities <- syn("syn20830825") %>%
  read_rds()

# set toolscore function -------------------------------------------------------
###############################################################################T

calc_toolscore<-function(data, lspci_id, gene_id)
{
  example_subset<-data[data$lspci_id==lspci_id,]

  ontarget<-example_subset[example_subset$gene_id == gene_id,]
  offtarget<-example_subset[example_subset$gene_id !=  gene_id,]

  chembl_active_data<-as.numeric(as.character(ontarget$standard_value))
  chembl_active_low_nM<-chembl_active_data[chembl_active_data<=100]
  if (length(chembl_active_low_nM) > 1) {
    chembl_strength<-7
  }else if (length(chembl_active_data[chembl_active_data<=1000]) > 4) {
    chembl_strength<-4
  }else {
    chembl_strength<-1
  }

  strength<-chembl_strength + 1 # bonus for multiple sources

  ontarget_IC50<-as.numeric(as.character(ontarget[!is.na(ontarget$standard_value),]$standard_value))
  ontarget_IC50_N<-length(ontarget_IC50)
  ontarget_IC50_Q1<-quantile(ontarget_IC50, probs=0.25, na.rm=TRUE)
  offtarget_IC50<-as.numeric(as.character(offtarget[!is.na(offtarget$standard_value),]$standard_value))
  offtarget_IC50_N<-length(offtarget_IC50)
  offtarget_IC50_Q1<-quantile(offtarget_IC50, probs=0.25, na.rm=TRUE)
  Q1_IC50_diff<-log10(offtarget_IC50_Q1)-log10(ontarget_IC50_Q1);
  if (length(ontarget_IC50)==0 || length(offtarget_IC50)==0) {
    wilcox_pval<-10;
  } else {
    w<-wilcox.test(ontarget_IC50, offtarget_IC50, alternative='less');
    if (w$p.value>=1E-15){
      wilcox_pval<-w$p.value;
    } else {
      wilcox_pval<-1E-15;
    }
  }

  IC50<-c(ontarget_IC50, offtarget_IC50)
  IC50_Q1<-quantile(IC50, probs=0.25, na.rm=TRUE)
  investigation_bias<-length(ontarget_IC50)/length(IC50)

  selectivity<-(Q1_IC50_diff/3 + (1-investigation_bias) - log10(wilcox_pval)/15)/3
  tool_score<- strength * selectivity
  #names(tool_score)<-'tool score'
  return (list(tool_score[[1]],strength,selectivity[[1]],investigation_bias,wilcox_pval,
               Q1_IC50_diff[[1]],offtarget_IC50_Q1[[1]],ontarget_IC50_Q1[[1]],ontarget_IC50_N,offtarget_IC50_N))
}

iterate_targets <- function(c.data) {
  toolscore<-list()
  toolscore_index<-0
  # c.data<-as.data.frame(lspci_id_list[[index_cmpd]])
  c.lspci_id<-unique(c.data$lspci_id)
  c.targets<-unique(c.data$gene_id)
  for(index_target in 1:length(c.targets)){
    c.df<-list()
    c.gene_id<-c.targets[index_target]
    c.df$lspci_id<-c.lspci_id
    c.df$gene_id<-c.gene_id
    c.return<-calc_toolscore(c.data, c.lspci_id, c.gene_id)
    if(!is.na(c.return[[3]])){
      toolscore_index<-toolscore_index+1
      c.df$tool_score<-c.return[[1]]
      c.df$strength<-c.return[[2]]
      c.df$selectivity<-c.return[[3]]
      c.df$investigation_bias<-c.return[[4]]
      c.df$wilcox_pval<-c.return[[5]]
      c.df$IC50_diff<-c.return[[6]]
      c.df$ontarget_IC50_Q1<-c.return[[8]]
      c.df$offtarget_IC50_Q1<-c.return[[7]]
      c.df$ontarget_IC50_N<-c.return[[9]]
      c.df$offtarget_IC50_N<-c.return[[10]]
      c.df$N_total<-c.return[[9]]+c.return[[10]]
      toolscore[[toolscore_index]]<-as.data.frame(c.df)
    }
  }
  # print(paste0(index_cmpd,"-",length(lspci_id_list)))
  # print(paste("Done", c.lspci_id))
  toolscore
}


# calculate toolscores ---------------------------------------------------------
###############################################################################T



activities_lspci_id_geneid <- activities %>%
  unnest(data) %>%
  rename(lspci_id = eq_class, gene_id = entrez_gene_id, standard_value = value)

# lspci_id_list<-dlply(activities_lspci_id_geneid,.(lspci_id),c)

plan(multisession(workers = 8))
toolscore.b <- activities_lspci_id_geneid %>%
  group_nest(fp_name, fp_type, lspci_id, keep = TRUE) %>%
  mutate(
    result = future_map(
      data,
      iterate_targets,
      .progress = TRUE
    )
  )

write_rds(
  toolscore.b,
  file.path(dir_release, "tool_score_raw.rds"),
  compress = "gz"
)

# toolscore.b <- read_rds(file.path(dir_release, "tool_score_raw.rds"))

toolscore_all <- toolscore.b %>%
  select(-data) %>%
  as.data.table() %>%
  .[
    ,
    .(result = map(
      result,
      compose(bind_rows, as_tibble, .dir = "forward")
    ) %>%
      bind_rows() %>%
      list()),
    keyby = .(fp_name, fp_type)
  ] %>%
  as_tibble()


write_rds(
  toolscore_all,
  file.path(dir_release, "selectivity.rds"),
  compress = "gz"
)

# toolscore_all <- read_rds(file.path(dir_release, "selectivity.rds"))

# Calculating selectivity classes ----------------------------------------------
###############################################################################T

calculate_most_selective <- function(df) {
  df %>%
    filter(
      tool_score >= 5,
      IC50_diff >= 2,
      wilcox_pval <= 0.1,
      strength == 8,
      investigation_bias <= 0.2
    )
}

calculate_semi_selective <- function(df, most_selective) {
  df %>%
    filter(
      IC50_diff >= 1,
      wilcox_pval <= 0.1,
      strength >= 5,
      investigation_bias <= 0.2
    ) %>%
    anti_join(
      most_selective,
      by = c("lspci_id", "gene_id")
    )
}

calculate_poly_selective <- function(df, most_selective, semi_selective) {
  genes_already_covered <- union(
    most_selective$gene_id,
    semi_selective$gene_id
  )
  df %>%
    filter(
      ontarget_IC50_N > 1,
      investigation_bias <= 0.2,
      IC50_diff >= 0,
      ontarget_IC50_Q1 < 9000,
      gene_id %nin% genes_already_covered
    ) %>%
    anti_join(
      bind_rows(most_selective, semi_selective),
      by = c("lspci_id", "gene_id")
    ) %>%
    arrange(desc(tool_score)) %>%
    group_by(gene_id) %>%
    slice(1) %>%
    ungroup()
}

calculate_unknown_selective <- function(df, most_selective, semi_selective, poly_selective) {
  df %>%
    filter(
      IC50_diff >= 0,
      ontarget_IC50_Q1 < 9000
    ) %>%
    anti_join(
      bind_rows(most_selective, semi_selective, poly_selective),
      by = c("lspci_id", "gene_id")
    )
}

selectivity_classes <- toolscore_all %>%
  mutate(
    most_selective = map(
      result,
      calculate_most_selective
    )
  ) %>%
  mutate(
    semi_selective = map2(
      result, most_selective,
      calculate_semi_selective
    )
  ) %>%
  mutate(
    poly_selective = pmap(
      list(result, most_selective, semi_selective),
      calculate_poly_selective
    )
  ) %>%
  mutate(
    unknown_selective = pmap(
      list(result, most_selective, semi_selective, poly_selective),
      calculate_unknown_selective
    )
  )

write_rds(
  selectivity_classes %>%
    select(-result),
  file.path(dir_release, "selectivity_classes.rds"),
  compress = "gz"
)

selectivity_classes <- read_rds(file.path(dir_release, "selectivity_classes.rds"))

# Combine classes in single table ----------------------------------------------
###############################################################################T

canonical_table <- syn("syn20835543") %>%
  read_rds()

target_table <- syn("syn20693721") %>%
  read_csv(col_types = "ciciccciccccc")

selectivity_classes_long <- selectivity_classes %>%
  gather("selectivity_class", "selectivity", ends_with("_selective")) %>%
  group_by(fp_name, fp_type) %>%
  summarize(
    selectivity = selectivity %>%
      set_names(selectivity_class) %>%
      bind_rows(.id = "selectivity_class") %>%
      list()
  ) %>%
  left_join(
    canonical_table %>%
      rename(canonical_table = data) %>%
      mutate(
        canonical_table = map(canonical_table, select, lspci_id, chembl_id, hms_id, pref_name)
      )
  ) %>%
  mutate(
    selectivity = map2(
      selectivity, canonical_table,
      left_join
    ) %>%
      map(
        left_join,
        target_table %>%
          distinct(gene_id = entrez_gene_id, gene_symbol = entrez_symbol)
      ) %>%
      map(
        select,
        selectivity_class,
        lspci_id, chembl_id, hms_id, pref_name, gene_id, gene_symbol, tool_score, strength, selectivity,
        investigation_bias, wilcox_pval, IC50_diff, ontarget_IC50_Q1, offtarget_IC50_Q1,
        ontarget_IC50_N , offtarget_IC50_N, N_total
      )
  )

pwalk(
  selectivity_classes_long,
  function(fp_name, selectivity, ...) {
    write_csv(selectivity, file.path(dir_release, paste0("selectivity_classes_", fp_name, ".csv.gz")))
  }
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

selectivity_calc_activity <- Activity(
  name = "Calculate compound-target selectivity and selectivity classes",
  used = c(
    "syn20693721",
    "syn20830825",
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/04_calculating_selectivity.R"
)

syn_selectivity_folder <- Folder("selectivity", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  # file.path(dir_release, "selectivity.rds"),
  # file.path(dir_release, "selectivity_classes.rds"),
  Sys.glob(file.path(dir_release, "selectivity_classes_*.csv.gz"))
) %>%
  synStoreMany(parentId = syn_selectivity_folder, activity = selectivity_calc_activity)

