# Script to standardize taxonomic names for the SIREN project

# Reads in a species list generated from multiple databases and returns a non-redundant,
# standardized list in which all taxonomic names are standardized based on one database

# Written by Chris Hempel (christopher.hempel@kaust.edu.sa) on 31 Oct 2023

library(bdc) # For taxonomic standardization
library(duckdb) # Required for bdc
library(openxlsx) # To read in Excel file
library(dplyr)
library(dplyr)

# Edit these options as needed
## Database to use for taxonomy (for all options, check out https://brunobrr.github.io/bdc/articles/taxonomy.html)
db <- "gbif" # Note: GBIF seemed to perform best. NCBI just accepted everything
## Excel file containing species list
species_list_excel_file <- "/Users/christopherhempel/Google Drive/KAUST/SIREN project/Red Sea species list.xlsx"
## Sheet in Excel file containing species list
species_list_sheet <- "Species cleaned + DBs"

# Read in df
df <- read.xlsx(species_list_excel_file, sheet = species_list_sheet)
species_list <- df$Species
# Standardize taxonomy
query_names <- bdc_query_names_taxadb(
  sci_name            = species_list,
  replace_synonyms    = FALSE, # replace synonyms by accepted names?
  suggest_names       = TRUE, # try to find a candidate name for misspelled names?
  db                  = db, # taxonomic database
  parallel            = TRUE, # should parallel processing be used?
  ncores              = 2, # number of cores to be used in the parallelization process
)

query_names$database <- df$Database

# Extract species names that were not found in the database
species_df_not_found <- query_names[query_names$notes == "notFound", c("original_search", "database")]
names(species_df_not_found)[names(species_df_not_found) == "original_search"] <- "species_not_found"

# Generate search info
num_total <- length(species_list)
num_accepted <- sum(grepl("accepted", query_names$notes))
num_misspelled <- sum(grepl("wasMisspelled", query_names$notes))
num_synonym <- sum(grepl("replaceSynonym", query_names$notes))
num_notfound <- length(rownames(species_df_not_found))
cat("Total number of species in list:", num_total, "\n")
cat("Number of accepted species:", num_accepted, "\n")
cat("Number of misspelled species:", num_misspelled, "\n")
cat("Number of synonym species:", num_synonym, "\n")
cat("Number of species not found:", num_notfound, "\n")

# Create a new data frame for accepted species with tax info
species_df_found <- query_names[grepl("accepted", query_names$notes), c("kingdom", "phylum", "class", "order", "family", "genus", "scientificName", "database", "notes")]
# Rename the scientificName and notes column
names(species_df_found)[names(species_df_found) == "scientificName"] <- "species"
names(species_df_found)[names(species_df_found) == "notes"] <- "validation_status"
# Count and remove duplicate species, keep database info, and collect number info
num_dupl <- sum(duplicated(species_df_found$species))
# Split and unnest the 'database' column
species_df_found <- species_df_found %>%
  separate_rows(database, sep = ",\\s*") %>%
  mutate(database = trimws(database))
# Group by 'species', combine unique databases, and retain other columns
species_df_found <- species_df_found %>%
  group_by(species) %>%
  summarise_all(function(x) paste(unique(x), collapse = ", ")) %>%
  ungroup()
cat("Number of duplicate species:", num_dupl, "\n")
num_total_after_harmonization <- length(rownames(species_df_found))
cat("Number of unique species after harmonization:", num_total_after_harmonization, "\n")
                                        
# Save the info
search_info <- as.data.frame(t(as.data.frame(c(num_total, num_accepted, num_misspelled, num_synonym, num_notfound, num_dupl, num_total_after_harmonization))))
colnames(search_info) <- c("num_total", "num_accepted", "num_misspelled", "num_synonym", "num_notfound", "num_dupl", "num_total_after_harmonization")
rownames(search_info) <- NULL

# Export dfs
species_df_found <- species_df_found %>%
  select(kingdom, phylum, class, order, family, genus, species, database, validation_status) %>%
  arrange(kingdom, phylum, class, order, family, genus, species, database,validation_status)
write.csv(species_df_found, "/Users/christopherhempel/Google Drive/KAUST/SIREN project/red_sea_species_list_standardized.csv", row.names = FALSE)
write.csv(species_df_not_found, "/Users/christopherhempel/Google Drive/KAUST/SIREN project/red_sea_species_list_not_found_in_gbif.csv", row.names = FALSE)
write.csv(search_info, "/Users/christopherhempel/Google Drive/KAUST/SIREN project/red_sea_species_list_standardization_info.csv", row.names = FALSE)

