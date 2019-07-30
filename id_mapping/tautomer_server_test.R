library(tidyverse)
library(httr)

x <- POST(
  "http://127.0.0.1:5000/query/tautomers",
  body = list(
    "smiles" = "C[C@@H](C1=C(C(=O)C2=C(O1)C=CC(=C2)F)C3=CC(=CC=C3)F)N4C5=C(C(=N4)C6=CC(=C(C=C6)OC(C)C)F)C(=NC=N5)N"
  ),
  encode = "multipart",
  accept_json()
)

content(x, as = "text") %>%
  jsonlite::fromJSON()
# $request
# $request$inchi_key
# [1] "IUVCFHHAEHNCFT-INIZCTEOSA-N"
#
# $request$smiles
# [1] "C[C@@H](C1=C(C(=O)C2=C(O1)C=CC(=C2)F)C3=CC(=CC=C3)F)N4C5=C(C(=N4)C6=CC(=C(C=C6)OC(C)C)F)C(=NC=N5)N"
#
#
# $tautomers
# inchi
# 1 InChI=1S/C31H24F3N5O3/c1-15(2)41-24-9-7-18(12-22(24)34)27-26-30(35)36-14-37-31(26)39(38-27)16(3)29-25(17-5-4-6-19(32)11-17)28(40)21-13-20(33)8-10-23(21)42-29/h4-16H,1-3H3,(H2,35,36,37)/t16-/m0/s1
# 2 InChI=1S/C31H24F3N5O3/c1-15(2)41-24-9-7-18(12-22(24)34)27-26-30(35)36-14-37-31(26)39(38-27)16(3)29-25(17-5-4-6-19(32)11-17)28(40)21-13-20(33)8-10-23(21)42-29/h4-16H,1-3H3,(H2,35,36,37)/t16-/m0/s1
# 3 InChI=1S/C31H24F3N5O3/c1-15(2)41-24-9-7-18(12-22(24)34)27-26-30(35)36-14-37-31(26)39(38-27)16(3)29-25(17-5-4-6-19(32)11-17)28(40)21-13-20(33)8-10-23(21)42-29/h4-16H,1-3H3,(H2,35,36,37)/t16-/m0/s1
# 4         InChI=1S/C31H24F3N5O3/c1-15(2)41-24-9-7-18(12-22(24)34)27-26-30(35)36-14-37-31(26)39(38-27)16(3)29-25(17-5-4-6-19(32)11-17)28(40)21-13-20(33)8-10-23(21)42-29/h4-16,35,38H,1-3H3/t16-/m0/s1
# inchi_key                                                                              smiles
# 1 IUVCFHHAEHNCFT-INIZCTEOSA-N      CC(C)Oc1ccc(-c2nn([C@@H](C)c3oc4ccc(F)cc4c(=O)c3-c3cccc(F)c3)c3ncnc(N)c23)cc1F
# 2 IUVCFHHAEHNCFT-INIZCTEOSA-N  CC(C)Oc1ccc(-c2nn([C@@H](C)c3oc4ccc(F)cc4c(=O)c3-c3cccc(F)c3)c3nc[nH]c(=N)c23)cc1F
# 3 IUVCFHHAEHNCFT-INIZCTEOSA-N  CC(C)Oc1ccc(-c2nn([C@@H](C)c3oc4ccc(F)cc4c(=O)c3-c3cccc(F)c3)c3[nH]cnc(=N)c23)cc1F
# 4 WPKFOPKHYJMHNR-INIZCTEOSA-N CC(C)Oc1ccc(-c2[nH]n([C@@H](C)c3oc4ccc(F)cc4c(=O)c3-c3cccc(F)c3)c3ncnc(=N)c2-3)cc1F



x <- POST(
  "http://127.0.0.1:5000/query/convert",
  body = list(
    "in" = "smiles",
    "out" = "inchi",
    "value" = "C[C@@H](C1=C(C(=O)C2=C(O1)C=CC(=C2)F)C3=CC(=CC=C3)F)N4C5=C(C(=N4)C6=CC(=C(C=C6)OC(C)C)F)C(=NC=N5)N"
  ),
  encode = "multipart",
  accept_json()
)

content(x, "text") %>%
  jsonlite::fromJSON()
# $converted
# $converted$inchi
# [1] "InChI=1S/C31H24F3N5O3/c1-15(2)41-24-9-7-18(12-22(24)34)27-26-30(35)36-14-37-31(26)39(38-27)16(3)29-25(17-5-4-6-19(32)11-17)28(40)21-13-20(33)8-10-23(21)42-29/h4-16H,1-3H3,(H2,35,36,37)/t16-/m0/s1"
#
#
# $request
# $request$`in`
# [1] "smiles"
#
# $request$out
# [1] "inchi"
#
# $request$value
# [1] "C[C@@H](C1=C(C(=O)C2=C(O1)C=CC(=C2)F)C3=CC(=CC=C3)F)N4C5=C(C(=N4)C6=CC(=C(C=C6)OC(C)C)F)C(=NC=N5)N"

