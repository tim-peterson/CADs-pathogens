



install.packages("dbparser")
## load dbparser package
library(dbparser)
library(dplyr)
library(ggplot2)
library(XML)

## parse data from XML and save it to memory
read_drugbank_xml_db("/Users/timpeterson/Downloads/full database.xml")

## load drugs data
drugs <- drug()


run_all_parsers(
  save_table = FALSE,
  save_csv = TRUE,
  csv_path = "/Volumes/GoogleDrive/My Drive/DATA/DrugBank",
  override_csv = FALSE,
  database_connection = NULL
)
