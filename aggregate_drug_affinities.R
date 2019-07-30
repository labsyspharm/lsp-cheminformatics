#this script prepares allbiochem data for BAFP clustering

library(plyr);library(dplyr)
library(tidyr)
options(scipen=999)

#############################################################################################################T
# set directories & import files----------
#############################################################################################################T
dir_chembl_v24<-"/Users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/ChEMBL_24_1"
dir_klaegerdata<-"/Users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/InHouse_Data/Klaeger_etAl_2017"
dir_doseresponse_inhouse<-"/Users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/InHouse_Data"
dir_RT<-"/Users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/Reagenttracker"
dir_cheminformatics_files<-"/Users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files"

setwd(dir_RT)
map_RT_hmsID<-read.csv("rt_chembl_matches_20190617.txt")
names(map_RT_hmsID)[1]<-"hms_id"

setwd(dir_klaegerdata)
klaegerdata<-read.csv("Klaeger_reformatted.csv",stringsAsFactors = F)
klaeger_neat<-klaegerdata%>%filter(chembl_id != "")

setwd(dir_chembl_v24)
#fragment_table<-read.csv("BAFP_chemblID_fragmentID_20190707.csv")
biochem_alldata<-read.csv("biochemicaldata_allcmpds_chembl24_1.csv")
cmpd_info<-read.csv("all_compounds_chembl24_1.csv",stringsAsFactors = F)
names(cmpd_info)[2]<-"pref_name_cmpd"
biochem_neat<-biochem_alldata%>%filter(tax_id==9606, is.na(gene_id)==FALSE)%>%
  merge(.,cmpd_info[c("molregno","pref_name_cmpd","chembl_id")], by="molregno")%>%
  unique()

setwd(dir_doseresponse_inhouse)
doseresponse_inhouse<-read.csv("all_doseresponse_summaryfile_20190417_vF.csv",stringsAsFactors = F)
doseresponse_inhouse_neat<-doseresponse_inhouse%>%
  mutate(hms_id=sub("1","HMSL1",hms_id))%>%
  merge(.,map_RT_hmsID, by="hms_id")

#############################################################################################################T
# bind files and make Q1 file----------
#############################################################################################################T

names(klaeger_neat)
names(biochem_neat)
names(doseresponse_inhouse_neat)
klaeger_rowbind<-klaeger_neat%>%
  mutate(chembl_id_cmpd=chembl_id)%>%
  select(chembl_id_cmpd, value, value_unit, gene_id, reference_id, reference_type, file_url)
biochem_rowbind<-biochem_neat%>%
  mutate(reference_id=chembl_id_doc,
         reference_type="chembl_id",
         file_url=paste0("https://www.ebi.ac.uk/chembl/document_report_card/",chembl_id_doc),
         chembl_id_cmpd=chembl_id,
         value=standard_value,
         value_unit=standard_units)%>%
  select(chembl_id_cmpd, value, value_unit, gene_id, reference_id, reference_type, file_url)
doseresponse_inhouse_rowbind<-doseresponse_inhouse_neat%>%
  mutate(chembl_id_cmpd=chembl_id,
         reference_id=synapse_id,
         reference_type="synapse_id")%>%
  select(chembl_id_cmpd,value,value_unit, gene_id, reference_id,reference_type,file_url)

complete_table<-bind_rows(klaeger_rowbind,biochem_rowbind,doseresponse_inhouse_rowbind)%>%
  unique()%>%
  mutate(value=round(value,2))%>%
  filter(value>0)%>%
  filter(is.na(value)==F)
setwd(dir_cheminformatics_files)
saveRDS(complete_table,"complete_inhouse_klaeger_chembl_biochem_20190708.RDS")

#test: unit
complete_table$value_unit%>%unique()

options(scipen=0)
complete_table_Q1<-complete_table%>%
  group_by(chembl_id_cmpd,gene_id)%>%
  summarise(Q1=round(quantile(value)[[2]]))
setwd(dir_cheminformatics_files)
saveRDS(complete_table_Q1,"Q1_table_alldata_20190708.RDS")
write.csv(complete_table_Q1,"Q1_table_alldata_20190708.csv",row.names = F)

complete_table_binding<-complete_table_Q1%>%
  mutate(binding=ifelse(Q1<10000,1,0))
setwd(dir_cheminformatics_files)
saveRDS(complete_table_binding, "binding_table_alldata_20190708.RDS")
write.csv(complete_table_binding,"binding_table_alldata_20190708.csv",row.names = F)
