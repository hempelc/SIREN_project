import pandas as pd
from Bio import Entrez

# Set your email address for the Entrez API
Entrez.email = <email-address>
infile = "red_sea_species_list_manually_cleaned.csv"
outfile = "red_sea_species_list_manually_cleaned_with_coi_counts.csv"

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

counts_dict2 = {}
while True:
    try:
        for species in df["Species"]:
            # Search for COI sequences in GenBank for the current species
            if species in counts_dict2:
                continue
            query = f'"{species}"[Organism] AND (COI[All Fields] OR CO1[All Fields] OR COXI[All Fields] OR COX1[All Fields] OR "cytochrome oxidase subunit 1"[All Fields])'
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
            record = Entrez.read(handle)
            counts_dict2[species] = int(record["Count"])
    except Exception:
        continue
    if len(counts_dict2) == len(df):
        break


# Convert count dictionary to DataFrame
counts_df2 = pd.DataFrame.from_dict(
    counts_dict2, orient="index", columns=["Number of COI sequences on GenBank -v2"]
)
counts_df2["Species"] = counts_df2.index
counts_df2 = counts_df2.reset_index(drop=True)


# Merge info and export
df_with_counts = pd.merge(df, counts_df, how="outer", on="Species")
df_with_counts = pd.merge(df_with_counts, counts_df2, how="outer", on="Species")

df_with_counts_diff = df_with_counts[
    df_with_counts["Number of COI sequences on GenBank"]
    != df_with_counts["Number of COI sequences on GenBank -v2"]
]

df_with_counts_diff.to_csv(outfile, index=False)
