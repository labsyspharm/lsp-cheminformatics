library(tidyverse)
library(httr)
library(jsonlite)

cmpds <- read_csv("../chembl24_1/hms_lincs_chembl_mapping.csv.gz") %>%
  select(compound = inchi, name = chembl_id) %>%
  head(10) %>%
  bind_rows(
    head(., 2) %>%
      mutate(name = paste0("copy_", name))
  )

request_json <- list(
  request = list(encoding = "inchi"),
  compounds = cmpds
) %>%
  toJSON(dataframe = "columns", auto_unbox = TRUE)

x <- POST(
  "http://127.0.0.1:8000/fingerprints/fingerprint_db",
  body = request_json,
  content_type_json(),
  write_disk("test.fps", overwrite = TRUE)
)


y <- POST(
  "http://127.0.0.1:8000/fingerprints/simsearch",
  body = list(
    "threshold" = 0.7,
    "fingerprint_db" = upload_file("test.fps")
  ),
  encode = "multipart",
  accept_json()
)

y_parsed <- content(y, as = "text") %>%
  fromJSON() %>%
  as_tibble()

