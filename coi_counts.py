import pandas as pd
from Bio import Entrez
from time import sleep

# Set your email address for the Entrez API
Entrez.email = "christopher.hempel@kaust.edu.sa"
infile = "/Users/christopherhempel/Google Drive/KAUST/SIREN project/red_sea_species_list_standardized.csv"
outfile = "/Users/christopherhempel/Google Drive/KAUST/SIREN project/red_sea_species_list_standardized_with_coi_counts.csv"


# Function to get counts of COI sequences from species list
def count_coi_sequences(species_list):
    counts_dict = {}
    for species in species_list:
        # Search for COI sequences in GenBank for the current species
        query = f'"{species}"[Organism] AND (COI[All Fields] OR COX1[All Fields] OR "cytochrome oxidase subunit 1"[All Fields])'
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
        record = Entrez.read(handle)
        counts_dict[species] = int(record["Count"])
    return counts_dict


# Import df and clean
df = pd.read_csv(infile)

# Grab species and determine number of COI sequences on GenBank (=counts)
species_list = df["species"].drop_duplicates()
# get counts of COI sequences from species list
counts_dict = {}
while True:
    try:
        for species in species_list:
            # Search for COI sequences in GenBank for the current species
            if species in counts_dict:
                continue
            query = f'"{species}"[Organism] AND (COI[All Fields] OR COX1[All Fields] OR "cytochrome oxidase subunit 1"[All Fields])'
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
            record = Entrez.read(handle)
            counts_dict[species] = int(record["Count"])
    except:
        continue
    if len(counts_dict) == 6766:
        break


# Convert count dictionary to DataFrame
counts_df = pd.DataFrame.from_dict(
    counts_dict, orient="index", columns=["Number of COI sequences on GenBank"]
)
counts_df["species"] = counts_df.index
counts_df = counts_df.reset_index(drop=True)

# Merge info and export
s = pd.merge(df, counts_df, how="outer", on="species")
df_with_counts.to_csv(outfile, index=False)
