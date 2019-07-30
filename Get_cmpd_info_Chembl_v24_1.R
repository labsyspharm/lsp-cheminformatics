# get infor for smallmoleculesuite 2.0 from chembl


library(tidyverse)
library(vroom)
library(data.table)
library(biomaRt)
library(bit64)
# options(scipen = 99999999)

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

# Get target gene mapping
map_chembl_tid_gene_id <- read_csv(
  "map_targets_chemblID_genID_uniprot_tid.csv.gz",
  col_types = "ccccciciiccccc"
)


################################################################################################################################################################################################T
# get biochemical data all compounds ------------
################################################################################################################################################################################################T
setwd(dir_chembl)
all_cmpds <- read_csv(
  "all_compounds_chembl24_1_parent_mapped.csv.gz",
  col_types = "icciciicci"
)

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
                        and A.bao_format not in ('BAO_0000221', 'BAO_0000219','BAO_0000218') LIMIT 100"))


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

activities_biochem_geneid <- activities_biochem %>%
  left_join(
    map_chembl_tid_gene_id %>%
      distinct(chembl_id_target = chembl_id, symbol, gene_id, uniprot_id, tax_id, tid, target_type, pref_name) %>%
      mutate_at(vars(tid), as.integer64),
    by = "tid"
  ) %>%
  as_tibble()

write_csv(activities_biochem_geneid, "biochemicaldata_allcmpds_chembl24_1.csv.gz")


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
# in the query, result should be empty
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
write_csv(pheno_activities, "phenotypic_assaydata_allcmpds_chemblv24_1.csv.gz")


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
write_csv(approval_info, "approval_info_phase1to4cmpds_ChemblV24_1.csv.gz")

# dbGetQuery(con, paste0("select * from drug_indication as DRUGIND
#                         left join indication_refs as INDREF on DRUGIND.drugind_id=INDREF.drugind_id"))%>%View


##################################################################################################################T
# get crossref info ------------
##################################################################################################################T
## goal: get temp_id to crossref doc
#step1: import crossref docs


# columns in Xref table: chembl_id|external_id|external_id_type|external_url|chembl_url











# Unichem_sources


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



