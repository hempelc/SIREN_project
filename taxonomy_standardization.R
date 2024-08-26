# Script to standardize taxonomic names for the SIREN project

# Reads in a species list generated from multiple databases and returns a non-redundant,
# standardized list in which all taxonomic names are standardized based on one database

# Written by Chris Hempel (christopher.hempel@kaust.edu.sa) on 31 Oct 2023

library(bdc) # For taxonomic standardization
library(duckdb) # Required for bdc
library(dplyr)
library(tidyr)
library(stringr)

# Edit these options as needed
## Database to use for taxonomy (for all options, check out https://brunobrr.github.io/bdc/articles/taxonomy.html)
db <- "gbif" # Note: Tried GBIF and NCBI. GBIF seemed to perform best. NCBI just accepted everything
## Excel file containing species list
species_list_file <- "red_sea_species_list_manually_cleaned_with_coi_counts.csv"

# Read in df
df <- read.csv(species_list_file)
species_list <- df$Species

# Standardize taxonomy
query_names <- bdc_query_names_taxadb(
  sci_name            = species_list,
  replace_synonyms    = TRUE, # replace synonyms by accepted names?
  suggest_names       = TRUE, # try to find a candidate name for misspelled names?
  db                  = db, # taxonomic database
  parallel            = TRUE, # should parallel processing be used?
  ncores              = 2, # number of cores to be used in the parallelization process
)

# Replace NaN with literal NA
query_names <- as.data.frame(apply(query_names, 2, function(x) ifelse(is.na(x), "NA", x)))

# Keep previous info
query_names$database <- df$Database
query_names$COI_sequences_on_GenBank <- df$Number.of.COI.sequences.on.GenBank

# Rename columns
names(query_names)[names(query_names) == "scientificName"] <- "species"
names(query_names)[names(query_names) == "notes"] <- "validation_status"

# Export relevant info of output as is ("raw")
raw_df <- query_names[c("original_search", "kingdom", "phylum", "class", "order", "family", "genus", "species", "database", "COI_sequences_on_GenBank", "validation_status")]
write.csv(raw_df, "red_sea_species_list_standardized_raw.csv", row.names = FALSE)

# Extract species names that were not found in the database
species_df_not_accepted <- query_names[!grepl("accepted", query_names$validation_status), c("original_search", "database", "COI_sequences_on_GenBank", "validation_status")]
names(species_df_not_accepted)[names(species_df_not_accepted) == "original_search"] <- "species_not_accepted"
write.csv(species_df_not_accepted, "red_sea_species_list_not_accepted_in_gbif.csv", row.names = FALSE)

# Create a new data frame for accepted species with tax and coi count info
species_df_accepted <- query_names[grepl("accepted", query_names$validation_status), c("kingdom", "phylum", "class", "order", "family", "genus", "species", "database", "COI_sequences_on_GenBank")]

# Capture duplication info
num_dupl <- sum(duplicated(species_df_accepted$species))

# Remove duplicate species, merge database info, and sum up sequences
species_df_accepted <- species_df_accepted %>%
  group_by(kingdom, phylum, class, order, family, genus, species) %>%
  summarise(
    database = paste(unique(unlist(strsplit(database, ", "))), collapse = ", "),
    COI_sequences_on_GenBank = sum(COI_sequences_on_GenBank)
  ) %>%
  ungroup()
# Export df
write.csv(species_df_accepted, "red_sea_species_list_standardized_synonyms_merged.csv", row.names = FALSE)

# Generate search info
## Replace validation statuses that have synonym in middle
query_names$validation_status <- gsub("homotypic synonym", "homotypicSynonym", query_names$validation_status)
query_names$validation_status <- gsub("heterotypic synonym", "heterotypicSynonym", query_names$validation_status)
query_names$validation_status <- gsub("proparte synonym", "proparteSynonym", query_names$validation_status)

## Generate a list of all encountered validation statuses
validation_status_list <- unique(unlist(strsplit(gsub(" ", "", as.character(query_names$validation_status)), "|", fixed = TRUE)))


## Create a new dataframe to store the results
info_df <- data.frame(
  observation = character(),
  count = numeric(),
  stringsAsFactors = FALSE
)
## Iterate over validation_status_list and count occurrences
for (validation_status in validation_status_list) {
  count <- sum(str_detect(query_names$validation_status, validation_status), na.rm = TRUE)
  info_df <- rbind(info_df, data.frame(observation = validation_status, count = count))
}

## Capture additional information, add to info_df, and save
num_total <- length(species_list)
num_notaccepted <- length(rownames(species_df_not_accepted))
num_total_after_standardization <- length(rownames(species_df_accepted))

info_df <- rbind(info_df, data.frame(observation = "number of species before standardization", count = num_total))
info_df <- rbind(info_df, data.frame(observation = "number of species not accepted in GBIF", count = num_notaccepted))
info_df <- rbind(info_df, data.frame(observation = "number of duplicated species after standardization", count = num_dupl))
info_df <- rbind(info_df, data.frame(observation = "number of unique, accepted species after standardization", count = num_total_after_standardization))

write.csv(info_df, "red_sea_species_list_standardization_info.csv", row.names = FALSE)

