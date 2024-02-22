import pandas as pd
from Bio import Entrez

# Set your email address for the Entrez API
Entrez.email = "christopher.hempel@kaust.edu.sa"
infile = "/Users/christopherhempel/Google Drive/KAUST/SIREN project data/red_sea_species_list_manually_cleaned.csv"
outfile = "/Users/christopherhempel/Google Drive/KAUST/SIREN project data/red_sea_species_list_manually_cleaned_with_coi_counts.csv"

# Import df and clean
df = pd.read_csv(infile)

# Get counts of COI sequences from species list
counts_dict = {}
while True:
    try:
        for species in df["Species"]:
            # Search for COI sequences in GenBank for the current species
            if species in counts_dict:
                continue
            query = f'"{species}"[Organism] AND (COI[All Fields] OR COX1[All Fields] OR "cytochrome oxidase subunit 1"[All Fields])'
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
            record = Entrez.read(handle)
            counts_dict[species] = int(record["Count"])
    except Exception:
        continue
    if len(counts_dict) == len(df):
        break


# Convert count dictionary to DataFrame
counts_df = pd.DataFrame.from_dict(
    counts_dict, orient="index", columns=["Number of COI sequences on GenBank"]
)
counts_df["Species"] = counts_df.index
counts_df = counts_df.reset_index(drop=True)

# Merge info and export
df_with_counts = pd.merge(df, counts_df, how="outer", on="Species")
df_with_counts.to_csv(outfile, index=False)
