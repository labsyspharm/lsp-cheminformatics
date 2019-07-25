# get infor for smallmoleculesuite 2.0 from chembl



library(plyr)
library(tidyverse)
library(data.table)
#options(scipen = 99999999)
##################################################################################################################T
# connect to chembl v. 24_1  ------------
##################################################################################################################T

require("RPostgreSQL")
## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")
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
# dir_UniChem<-"/Users/nienke/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/Xrefs_UniChem"

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

# map_uniprot_chembl %>%
#   filter(is.na(gene_id)) %>%
#   dplyr::count(organism)
# # A tibble: 3 x 2
# organism                n
# <chr>               <int>
# 1 Homo sapiens         18
# 2 Mus musculus        147
# 3 Rattus norvegicus   604

# Only 18 human uniprot IDs left without a matching Entrez ID, good enough...

# map_uniprot_chembl %>%
#   drop_na(gene_id) %>%
#   dplyr::count(gene_id) %>%
#   dplyr::count(n)
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
gene_info<-read.csv("gene_info_20190226_compact.csv")

names(chembl_info_unmapped_targets)
setwd(dir_chembl24_1)
gene_info_1<-read.csv("gene_info_20190226_addition_to_chemblinfo.csv")
#gene_info_1<-gene_info%>%
#  filter(Symbol %in% chembl_info_unmapped_targets$component_synonym)
target_not_gene_ifo<-chembl_info_unmapped_targets[! chembl_info_unmapped_targets$component_synonym %in% gene_info_1$Symbol,]
dim(target_not_gene_ifo)
head(target_not_gene_ifo)
tid_not_gene_info<-target_not_gene_ifo[!target_not_gene_ifo$tid %in% chembl_info_unmapped_targets$tid,]
dim(tid_not_gene_info) # success when 0,X
names(gene_info_1)<-c("tax_id","gene_id","symbol","synonyms","description","name_ncbi")
#setwd(dir_chembl24_1)
#write.csv(gene_info_1,file="gene_info_20190226_addition_to_chemblinfo.csv",row.names = F)


# step 5 --> merge gene_info to chembl_info_unmapped_targets on symbol&tax_id
unmapped_chembl_info<-chembl_info_unmapped_targets%>%
  merge(gene_info_1, by.x=c("tax_id","component_synonym"),by.y = c("tax_id","symbol"))%>%
  unique()
dim(chembl_info_unmapped_targets%>%unique())
dim(unmapped_chembl_info)
head(unmapped_chembl_info)

# step 6 --> merge doc step 5 with map chembl_id to gene_id
map_chemblID_geneID%>%head()
unique(map_chemblID_geneID$mapping_chembl2entrez_performed_by)
map_chemblID_geneID_2<-merge(map_chemblID_geneID,unmapped_chembl_info, by=c("chembl_id","target_type"),all.x = T )
names(map_chemblID_geneID_2)[c(5,7,15)]<-c("geneID_chemblv24","geneID_nmoretChemblv23","geneID_nmoretChemblv24")
names(map_chemblID_geneID_2)
map_chemblID_geneID_2[is.na(map_chemblID_geneID_2$geneID_nmoretChemblv24)==F,]$mapping_chembl2entrez_performed_by<-"nmoret_ncbi_merge (v24)"
map_chemblID_geneID_2[is.na(map_chemblID_geneID_2$geneID_chemblv24)==F,]$mapping_chembl2entrez_performed_by<-"nmoret_uniprot"
unique(map_chemblID_geneID_2$mapping_chembl2entrez_performed_by)

dim(map_chemblID_geneID_2%>%filter((is.na(tax_id)==F & is.na(mapping_chembl2entrez_performed_by)==T)))# zero lines
map_chemblID_geneID_3<-map_chemblID_geneID_2%>%filter(is.na(mapping_chembl2entrez_performed_by)==F)
unique(map_chemblID_geneID_3$mapping_chembl2entrez_performed_by)
dim(map_chemblID_geneID_3)

#step 7 --> prettyfy table
# reorder columns
names(map_chemblID_geneID_3)
map_chemblID_geneID_3<-map_chemblID_geneID_3%>%
  mutate(symbol=component_synonym)%>%
  select(chembl_id, symbol, uniprot_id, geneID_chemblv24, geneID_nmoretChemblv23,geneID_nmoretChemblv24,
            mapping_chembl2entrez_performed_by, mapping_uniprot2chemblID_performed_by,
            protein_name, gene_name,pref_name, name_ncbi, description, synonyms, organism, tax_id, tid, target_type) ## symbol is empty very often

#make tall& skinny version of gene id
map_chemblID_geneID_4<-map_chemblID_geneID_3%>%gather("geneID_source","gene_id",geneID_chemblv24:geneID_nmoretChemblv24,na.rm = T)%>%
  select(chembl_id, symbol, uniprot_id, gene_id, geneID_source,
         mapping_chembl2entrez_performed_by, mapping_uniprot2chemblID_performed_by,
         protein_name, gene_name,pref_name, name_ncbi, description, synonyms, organism, tax_id, tid, target_type)%>%
  arrange(chembl_id)
head(map_chemblID_geneID_4)

# make sure all rows have gene symbol, all rows tid, target_type, pref_name & organism
map_chemblID_geneID_5<-map_chemblID_geneID_4%>%
  merge(gene_info,by.x="gene_id", by.y="GeneID")%>%
  mutate(symbol=Symbol,
         name_ncbi=Full_name_from_nomenclature_authority,
         synonyms=Synonyms,
         tax_id=X.tax_id,
         description.x=description.y,
         description=description.x)%>%
  select(chembl_id,symbol,gene_id,uniprot_id,geneID_source,mapping_chembl2entrez_performed_by,mapping_uniprot2chemblID_performed_by,
         protein_name,gene_name,name_ncbi,description, synonyms, tax_id)%>% #tid, target_type, pref_name & organism in seperate query
  unique()


chembl_info_all_targets<-
    dbGetQuery(con, paste0("select dict.tid, dict.target_type, dict.pref_name, dict.organism, dict.chembl_id
                            from target_dictionary as dict
                            left join target_components as comp on comp.tid=dict.tid
                            and dict.chembl_id in ('",paste(map_chemblID_geneID_5$chembl_id%>%as.character(),collapse="','"),"')
                            "))%>%unique()

map_chemblID_geneID_5<-data.table(map_chemblID_geneID_5)
chembl_info_all_targets<-data.table(chembl_info_all_targets)

map_chemblID_geneID_6<-map_chemblID_geneID_5%>%
  merge(chembl_info_all_targets, by= "chembl_id")%>%
  mutate(tid=as.character(tid))

class(map_chemblID_geneID_6$tid)

setwd(dir_chembl24_1)
#write.csv(map_chemblID_geneID_6%>%unique(),file="map_targets_chemblID_genID_uniprot_tid.csv",row.names = F)
map_chemblID_geneID_6<-read.csv("map_targets_chemblID_genID_uniprot_tid.csv", stringsAsFactors = T)
map_chemblID_geneID_6[map_chemblID_geneID_6$tid=="1e+05",]$tid<-"100000"

##!! map to famplex id


##################################################################################################################T
# get table w/ all molecule structures, indicate if cmpd has different parental, create temp_id  ------------
##################################################################################################################T

# Chembl has salts sometimes separate from free bases, parent_molregno points to
# free base
# Trying to reconcile internal LINCS database, which only has a single ID for
# each compound, regardless of salt, with Chembl

# step 1--> get basic info on molecules
all_cmpds<-dbGetQuery(con, paste0("select dict.molregno, dict.pref_name, dict.chembl_id, dict.max_phase, syn.synonyms,
                                   hi.parent_molregno, hi.active_molregno, struct.standard_inchi
                                   from molecule_dictionary as dict
                                   left join molecule_synonyms as syn on dict.molregno = syn.molregno
                                   left join molecule_hierarchy as hi on dict.molregno = hi.molregno
                                   left join compound_structures as struct on dict.molregno=struct.molregno"))%>%
  unique()%>%
  mutate(parental_flag=ifelse(molregno!=parent_molregno,1,0))
#setwd(dir_chembl24_1)
#write.csv(all_cmpds, file="all_compounds_chembl24_1.csv",row.names = F)


    #step 1.1 create tempID201903_XXX for cmpds w/ parental cmpd
all_cmpds%>%filter(parental_flag==1)%>%dim()
all_cmpds%>%filter(parental_flag==0 | is.na(parental_flag)==T)%>%dim()
all_cmpds<-all_cmpds%>%
 # mutate(molregno_sub=gsub(1,"a",molregno))%>%
  mutate(temp_id=ifelse(parental_flag==0|is.na(parental_flag)==T,paste0("TempID201903_",formatC(format(molregno,scientific=F)%>%as.numeric(),format = "d",width=8,flag = "0")),
                        paste0("TempID201903_",formatC(format(parent_molregno,scientific=F)%>%as.numeric(),format = "d",width=8,flag = "0"))))

all_cmpds%>%filter(parental_flag==1)%>%select(temp_id)%>%unique%>%dim()#93901
all_cmpds%>%filter(parental_flag==1)%>%select(molregno)%>%unique%>%dim()#96503
all_cmpds%>%filter(parental_flag==0)%>%select(temp_id)%>%unique%>%dim()#1653230
all_cmpds%>%filter(parental_flag==0)%>%select(molregno)%>%unique%>%dim()#1653230
all_cmpds%>%select(temp_id)%>%unique()%>%dim()#1732317

tail(all_cmpds%>%filter(parental_flag==1))
all_cmpds[all_cmpds$temp_id=='TempID201903_01782648',]%>%View()

setwd(dir_chembl24_1)
write.csv(all_cmpds, file="all_compounds_chembl24_1_tempID1.csv",row.names = F)


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
                        and A.bao_format not in ('BAO_0000221', 'BAO_0000219','BAO_0000218')
                        LIMIT 100"))
View(biochem_test)


BAO_format<-dbGetQuery(con, paste0("select *
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
activities_biochem<-rbindlist(activities_biochem)%>%data.table(.)
activities_biochem$tid<-as.character(activities_biochem$tid)

map_chemblID_geneID_formerge<-map_chemblID_geneID_6%>%
  select(chembl_id,symbol,gene_id,uniprot_id,tax_id,tid,target_type,pref_name)%>%
  mutate(chembl_id_target=chembl_id)%>%
  select(-chembl_id)%>%
  unique()

activities_biochem_geneid<-
  left_join(activities_biochem,map_chemblID_geneID_formerge,by="tid")%>%
  as.data.frame(.)# multplicity of table length due to many-to-many gene_id to uniprot_id mapping


setwd(dir_chembl24_1)
write.csv(activities_biochem_geneid,file = "biochemicaldata_allcmpds_chembl24_1.csv",row.names = F)


################################################################################################################################################################################################T
# get phenotypic data all compounds ------------
################################################################################################################################################################################################T
standard_units<-dbGetQuery(con, paste0("
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

units_per_assay<-dbGetQuery(con, paste0(" select assay_id, count(distinct(standard_units)) as count_units,
                                        count(molregno) as count_molregno
                                        from activities
                                        where standard_units in ('",paste(standard_units_ok,collapse="','"),"')
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
                                     limit 100
                                     "))


dim(activities_2)

pheno_activities<-activities_1

pheno_activities$log10_value<-log10(pheno_activities$standard_value)
pheno_activities<-pheno_activities%>%filter(log10_value != "NaN")%>%filter(standard_units %in% standard_units_ok)

dim(pheno_activities)

setwd(dir_chembl24_1)
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
                              where max_phase >0 LIMIT 100"))
View(approval_info)

setwd(dir_chembl24_1)
write.csv(approval_info, file="approval_info_phase1to4cmpds_ChemblV24_1.csv",row.names = F)

dbGetQuery(con, paste0("select * from drug_indication as DRUGIND
                        left join indication_refs as INDREF on DRUGIND.drugind_id=INDREF.drugind_id
                          limit 1000 "))%>%View


##################################################################################################################T
# get crossref info ------------
##################################################################################################################T
## goal: get temp_id to crossref doc
#step1: import crossref docs

setwd(dir_UniChem)
list.files()

Unichem_sources<-read.delim("UniChem_SourceCode")

chembl2pdb<-read.delim("src1src3.txt")
names(chembl2pdb)<-c("chembl_id_compound", "Xref_id_compound")
chembl2pdb<-chembl2pdb%>%
  mutate(url=paste0(Unichem_sources$base_url[2],Xref_id_compound),
         Xref_type="pdb_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2gtopdb<-read.delim("src1src4.txt")
names(chembl2gtopdb)<-c("chembl_id_compound", "Xref_id_compound")
chembl2gtopdb<-chembl2gtopdb%>%
  mutate(url=paste0(Unichem_sources$base_url[3],Xref_id_compound),
         Xref_type="gtopdb_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2chebi<-read.delim("src1src7.txt")
names(chembl2chebi)<-c("chembl_id_compound", "Xref_id_compound")
chembl2chebi<-chembl2chebi%>%
  mutate(url=paste0(Unichem_sources$base_url[4],Xref_id_compound),
         Xref_type="chebi_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2zinc<-read.delim("src1src9.txt")
names(chembl2zinc)<-c("chembl_id_compound", "Xref_id_compound")
chembl2zinc<-chembl2zinc%>%
  mutate(url=paste0(Unichem_sources$base_url[5],Xref_id_compound),
         Xref_type="zinc_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2emolecules<-read.delim("src1src10.txt")
names(chembl2emolecules)<-c("chembl_id_compound","Xref_id_compound")
chembl2emolecules<-chembl2emolecules%>%
  mutate(url=paste0(Unichem_sources$base_url[6],Xref_id_compound),
         Xref_type="emolecules_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2fdasrs<-read.delim("src1src14.txt")
names(chembl2fdasrs)<-c("chembl_id_compound", "Xref_id_compound")
chembl2fdasrs<-chembl2fdasrs%>%
  mutate(url=paste0(Unichem_sources$base_url[7],Xref_id_compound),
         Xref_type="fdasrs_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2selleck<-read.delim("src1src20.txt")
names(chembl2selleck)<-c("chembl_id_compound", "Xref_id_compound")
chembl2selleck<-chembl2selleck%>%
  mutate(url=paste0(Unichem_sources$base_url[8],Xref_id_compound,".html"),
         Xref_type="selleck_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

chembl2pubchem<-read.delim("src1src22.txt")
names(chembl2pubchem)<-c("chembl_id_compound", "Xref_id_compound")
chembl2pubchem<-chembl2pubchem%>%
  mutate(url=paste0(Unichem_sources$base_url[9],Xref_id_compound),
         Xref_type="pubchem_id_compound",
         Xref_id_compound=as.character(Xref_id_compound))

head(chembl2chebi)
head(chembl2pdb)
head(chembl2fdasrs)
head(chembl2gtopdb)
head(chembl2pubchem)
head(chembl2selleck)
head(chembl2zinc)
head(chembl2emolecules)

list_Xref_total<-list()
list_Xref_total[[1]]<-chembl2chebi
list_Xref_total[[2]]<-chembl2emolecules
list_Xref_total[[3]]<-chembl2fdasrs
list_Xref_total[[4]]<-chembl2gtopdb
list_Xref_total[[5]]<-chembl2pdb
list_Xref_total[[6]]<-chembl2pubchem
list_Xref_total[[7]]<-chembl2selleck
list_Xref_total[[8]]<-chembl2zinc
Xref_total<-bind_rows(list_Xref_total)%>%arrange(chembl_id_compound)

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



