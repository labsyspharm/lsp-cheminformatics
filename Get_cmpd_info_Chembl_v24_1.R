# get infor for smallmoleculesuite 2.0 from chembl


library(tidyverse)
library(data.table)
#options(scipen = 99999999)
##################################################################################################################T
# connect to chembl v. 24_1  ------------
##################################################################################################################T

library(RPostgres)
## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_24",
                 host = "localhost", port = 5433,
                 user = "chembl_public")

##################################################################################################################T
# set directories  ------------
##################################################################################################################T
dir_chembl <- file.path("~", "repo", "tas_vectors", "chembl24_1")
# dir_chembl23<-"/Users/nienke/Dropbox (Personal)/Harvard/Novartis_Internship/atHMS/Chembl_V23"
dir_UniChem<- file.path("~", "repo", "tas_vectors", "unichem")

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

#write.csv(map_uniprot_chembl, file = "map_uniprotID_chemblID.csv",row.names = F)

chembl_info_all_targets <- dbGetQuery(
  con,
  paste0(
    "select dict.tid, dict.pref_name, dict.tax_id, dict.organism, dict.chembl_id
    from target_dictionary as dict"
  )
) %>%
  as_tibble() %>%
  filter(organism %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus"))

map_uniprot_chembl <- read_tsv(
  "chembl_uniprot_mapping.txt", skip = 1,
  col_names = c("uniprot_id", "chembl_id", "protein_name", "target_type")
) %>%
  inner_join(chembl_info_all_targets, by = "chembl_id")

library(biomaRt)
organism_mart_mapping <- c(
  "Homo sapiens" = "hsapiens_gene_ensembl",
  "Rattus norvegicus" = "rnorvegicus_gene_ensembl",
  "Mus musculus" = "mmusculus_gene_ensembl"
)


uniprot_to_entrez <- function(df, group) {
  mart <- useMart(
    "ENSEMBL_MART_ENSEMBL",
    organism_mart_mapping[group$organism]
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
      genebabel::join_hgnc(query_col = "uniprot_id", match_cols = "uniprot_ids", select_cols = "entrez_id") %>%
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
library(vroom)
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
  left_join(mutate_at(gene_info, "tax_id", bit64::as.integer64), by = c("tax_id", "gene_id")) %>%
  distinct()

write_csv(map_chemblID_geneID, "map_targets_chemblID_genID_uniprot_tid.csv.gz")

##################################################################################################################T
# get table w/ all molecule structures, indicate if cmpd has different parental, create temp_id  ------------
##################################################################################################################T

# Chembl has salts sometimes separate from free bases, parent_molregno points to
# free base
# Trying to reconcile internal LINCS database, which only has a single ID for
# each compound, regardless of salt, with Chembl

# step 1--> get basic info on molecules
all_cmpds <- dbGetQuery(
  con,
  paste0(
  "select distinct dict.molregno, dict.pref_name, dict.chembl_id, dict.max_phase, syn.synonyms,
    hi.parent_molregno, hi.active_molregno, struct.standard_inchi
   from molecule_dictionary as dict
   left join molecule_synonyms as syn on dict.molregno = syn.molregno
   left join molecule_hierarchy as hi on dict.molregno = hi.molregno
   left join compound_structures as struct on dict.molregno=struct.molregno"
  )
) %>%
  as_tibble() %>%
  mutate(parental_flag = ifelse(molregno != parent_molregno, 1, 0))
write.csv(all_cmpds, file="all_compounds_chembl24_1.csv",row.names = F)
# Again, molregno, parent_molregno and active_molregno are integer64...

    #step 1.1 create tempID201903_XXX for cmpds w/ parental cmpd
all_cmpds%>%filter(parental_flag==1)%>%dim()
all_cmpds%>%filter(parental_flag==0 | is.na(parental_flag)==T)%>%dim()

all_cmpds <- all_cmpds %>%
 # mutate(molregno_sub=gsub(1,"a",molregno))%>%
  mutate(
    temp_id = ifelse(
      parental_flag == 0 | is.na(parental_flag),
      paste0("TempID201903_", molregno),
      paste0("TempID201903_", parent_molregno)
    )
  )

all_cmpds%>%filter(parental_flag==1)%>%dplyr::select(temp_id)%>%unique%>%dim()#93901
all_cmpds%>%filter(parental_flag==1)%>%dplyr::select(molregno)%>%unique%>%dim()#96503
all_cmpds%>%filter(parental_flag==0)%>%dplyr::select(temp_id)%>%unique%>%dim()#1653230
all_cmpds%>%filter(parental_flag==0)%>%dplyr::select(molregno)%>%unique%>%dim()#1653230
all_cmpds%>%dplyr::select(temp_id)%>%unique()%>%dim()#1732317

tail(all_cmpds%>%filter(parental_flag==1))
all_cmpds[all_cmpds$temp_id=='TempID201903_01782648',]%>%View()

write_csv(all_cmpds, "all_compounds_chembl24_1_tempID1.csv")


################################################################################################################################################################################################T
# get biochemical data all compounds ------------
################################################################################################################################################################################################T
setwd(dir_chembl24_1)
all_cmpds<- read.csv("all_compounds_chembl24_1_tempID1.csv")


## redo biochem assays


biochem_test<-dbGetQuery(con, paste0("select *
                        from activities as ACT
                        left join assays as A
                        on ACT.assay_id = A.assay_id
                        where ACT.molregno in (",toString(unique(all_cmpds$molregno)),")
                        and ACT.standard_value is not null
                        and A.assay_type = 'B'
                        and A.relationship_type in ('D', 'H', 'M', 'U')
                        and ACT.standard_units = 'nM'
                        and ACT.standard_type in ('IC50','Ki','EC50','Kd','IC90','CC50','ID50','AC50','Inhibition','MIC','Potency','Activity','ED50')
                        and A.assay_cell_type is NULL
                        and A.bao_format not in ('BAO_0000221', 'BAO_0000219','BAO_0000218')"))


BAO_format <- dbGetQuery(con, paste0("select *
                        from bioassay_ontology
                                   "))#where bao_id in ('",paste(biochem_test$bao_format%>%unique,collapse="','"),"')
View(BAO_format)






activities_biochem_1<-dbGetQuery(con, paste0("select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno,ACT.standard_relation, ACT.standard_type,
                                             ACT.standard_value,ACT.standard_units,
                                             A.tid,
                                             A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc
                                             from activities as ACT
                                             left join assays as A
                                             on ACT.assay_id = A.assay_id
                                             left join DOCS
                                             on DOCS.doc_id=A.doc_id
                                             left join bioassay_ontology as BAO
                                             on A.bao_format=BAO.bao_id
                                             where ACT.molregno in (",toString(unique(all_cmpds$molregno)),")
                                             and ACT.standard_value is not null
                                             and A.assay_type = 'B'
                                             and A.relationship_type in ('D', 'H', 'M', 'U')
                                             and ACT.standard_units = 'nM'
                                             and ACT.standard_type in ('IC50','Ki','EC50','Kd','IC90','CC50','ID50','AC50','Inhibition','MIC','Potency','Activity','ED50')
                                             and A.assay_cell_type is NULL
                                             and A.bao_format not in ('BAO_0000221', 'BAO_0000219','BAO_0000218')
                                             "))
View(activities_biochem_1)

# Data from "Navigating the Kinome" paper is annotated using F (functional) assay type
# whereas most other data B (binding)
activities_biochem_2<-dbGetQuery(con, paste0("select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno,ACT.standard_relation, ACT.standard_type,
                                             ACT.standard_value,ACT.standard_units,
                                             A.tid,
                                             A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc
                                             from activities as ACT
                                             left join assays as A
                                             on ACT.assay_id = A.assay_id
                                             left join DOCS
                                             on DOCS.doc_id=A.doc_id
                                             left join bioassay_ontology as BAO
                                             on A.bao_format=BAO.bao_id
                                             where ACT.molregno in (",toString(unique(all_cmpds$molregno)),")
                                             and ACT.standard_value is not null
                                             and A.assay_type = 'F'
                                             and A.description like '%Navigating the Kinome%'
                                             "))

activities_biochem<-list(
  activities_biochem_1,
  activities_biochem_2
)
activities_biochem <- data.table::rbindlist(activities_biochem)

map_chemblID_geneID_formerge <- map_chemblID_geneID %>%
  dplyr::select(chembl_id, symbol, gene_id, uniprot_id, tax_id, tid, target_type, pref_name)%>%
  rename(chembl_id_target = chembl_id)%>%
  unique()

activities_biochem_geneid <-
  left_join(activities_biochem, map_chemblID_geneID_formerge, by="tid")%>%
  as_tibble()# multplicity of table length due to many-to-many gene_id to uniprot_id mapping


setwd(dir_chembl)
write.csv(activities_biochem_geneid,file = "biochemicaldata_allcmpds_chembl24_1.csv",row.names = F)


################################################################################################################################################################################################T
# get phenotypic data all compounds ------------
################################################################################################################################################################################################T
standard_units <- dbGetQuery(con, paste0("
                                       select distinct(standard_units)
                                       from activities as ACT
                                       where ACT.molregno in (",toString(all_cmpds$molregno),")"))

standard_units_ok<-c('M','mol/L','nM','nmol/L',
                     'nmol.L-1','pM','pmol/L','pmol/ml','um',
                     'uM','umol/L','umol/ml','umol/uL')



# assay_freq<-dbGetQuery(con, paste0(" select assay_id, count(molregno)
#                                    from activities
#                                    where standard_units in ('",paste(standard_units_ok,collapse="','"),"')
#                                    group by assay_id
#                                    "))

units_per_assay<-dbGetQuery(con, paste0("select assay_id, count(distinct(standard_units)) as count_units,
                                        count(molregno) as count_molregno
                                        from activities
                                        where standard_units in ('",paste(standard_units_ok, collapse="','"),"')
                                        group by assay_id
                                        "))
assays_qualified<-units_per_assay%>%filter(count_units==1)%>%filter(count_molregno>2)

activities_1<-dbGetQuery(con, paste0("select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno,ACT.standard_relation, ACT.standard_type,
                                             ACT.standard_value,ACT.standard_units,
                                             A.tid,
                                             A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc
                                             from activities as ACT
                                             left join assays as A
                                             on ACT.assay_id = A.assay_id
                                             left join DOCS
                                             on DOCS.doc_id=A.doc_id
                                             left join bioassay_ontology as BAO
                                             on A.bao_format=BAO.bao_id
                                     where ACT.standard_value is not null
                                     and A.relationship_type like 'N'
                                     and A.assay_id in (",toString(assays_qualified$assay_id),")
                                     and ACT.standard_units in ('",paste(standard_units_ok,collapse="','"),"')
                                     "))

# BAO list is a list of assays that are curated by Nienke to be interpretable
# https://bioportal.bioontology.org
# This is to check if we missed any of these assays with the BAO ids listed
# in the query
activities_2<-dbGetQuery(con, paste0("select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno,ACT.standard_relation, ACT.standard_type,
                                             ACT.standard_value,ACT.standard_units,
                                             A.tid,
                                             A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc
                                             from activities as ACT
                                             left join assays as A
                                             on ACT.assay_id = A.assay_id
                                             left join DOCS
                                             on DOCS.doc_id=A.doc_id
                                             left join bioassay_ontology as BAO
                                             on A.bao_format=BAO.bao_id
                                     where ACT.standard_value is not null
                                     and A.assay_id in (",toString(assays_qualified$assay_id),")
                                     and A.assay_id not in (",toString(activities_1$assay_id),")
                                     and ACT.standard_units in ('",paste(standard_units_ok,collapse="','"),"')
                                     and BAO.label in ('BAO_0000006','BAO_0000090','BAO_0000094','BAO_0000096','BAO_0000218','BAO_0000219','BAO_0000221')
                                     "))


dim(activities_2)

pheno_activities<-activities_1

pheno_activities$log10_value<-log10(pheno_activities$standard_value)
pheno_activities<-pheno_activities%>%drop_na(log10_value)%>%filter(standard_units %in% standard_units_ok)

dim(pheno_activities)

setwd(dir_chembl)
write.csv(file="phenotypic_assaydata_allcmpds_chemblv24_1.csv",pheno_activities, row.names=F)
saveRDS(pheno_activities,file = "phenotypic_assaydata_allcmpds_chemblv24_1.RDS")

##################################################################################################################################################################################################T
# get clinical &  info ------------
##################################################################################################################################################################################################T

approval_info<-dbGetQuery(con, paste0("select pref_name, chembl_id as chembl_id_compound, MOLDICT.molregno, max_phase, first_approval, oral, parenteral,
                              topical, black_box_warning,first_in_class, prodrug, indication_class, withdrawn_flag, withdrawn_year, withdrawn_country,
                              withdrawn_reason, max_phase_for_ind, mesh_id, mesh_heading, efo_id, efo_term, ref_type, ref_id, ref_url
                              from molecule_dictionary as MOLDICT
                              left join drug_indication as DRUGIND on MOLDICT.molregno = DRUGIND. molregno
                              left join indication_refs as INDREF on DRUGIND.drugind_id=INDREF.drugind_id
                              where max_phase >0"))
View(approval_info)

setwd(dir_chembl)
write.csv(approval_info, file="approval_info_phase1to4cmpds_ChemblV24_1.csv",row.names = F)

# dbGetQuery(con, paste0("select * from drug_indication as DRUGIND
#                         left join indication_refs as INDREF on DRUGIND.drugind_id=INDREF.drugind_id"))%>%View


##################################################################################################################T
# get crossref info ------------
##################################################################################################################T
## goal: get temp_id to crossref doc
#step1: import crossref docs

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

setwd(dir_UniChem)
write.csv(Xref_total,file="xref_table_unichem_20190325.csv",row.names = F)

# columns in Xref table: chembl_id|external_id|external_id_type|external_url|chembl_url











Unichem_sources


## has to be done via unichem
cross_reference_table

dbGetQuery(con, paste0("select *
                       from compound_records as COMPREC
                       left join docs on COMPREC.src_id=DOCS.src_id
                       where COMPREC.src_id = 20
                       limit 10"))%>%View()

#rec.src_id,rec.src_compound_id, rec.cidx, rec.compound_name, rec.compound_key,


# molregno, chembl_id, inchi,synonyms


#molecule_dictionary:
#molecule_synonyms: synonyms
#compound_records: SRC_ID, src_compound_id, CIDX, compound_name, compound_key
#molecule_hierarchy: parent_molregno, active_molregno

##################################################################################################################T
# get table w/ all 'same compound, different salt' IDs  ------------
##################################################################################################################T
#molecule hierarchy: parent_molregno, active_molregno


##################################################################################################################T
# give temp ID to each set   ------------
##################################################################################################################T



