import pandas as pd
import plotly.express as px
import os

# File with counts per species as .csv
infile = "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/Red Sea species list processing - chris - NCBI + BOLD counts/red_sea_species_list_standardized_synonyms_merged_unaccepted_merged.csv"

df = pd.read_csv(infile)

# Step 1: Convert 'BOLD' and 'GenBank' columns to binary
df['BOLD_binary'] = df['COI_sequences_on_BOLD'].apply(lambda x: 1 if x >= 1 else 0)
df['GenBank_binary'] = df['COI_sequences_on_GenBank'].apply(lambda x: 1 if x >= 1 else 0)

# Step 2: Count how many rows have a 1 in each column
count_bold = df['BOLD_binary'].sum()
count_genbank = df['GenBank_binary'].sum()

# Step 3: Count how many rows have a 1 in either 'BOLD' or 'GenBank'
count_either = df[['BOLD_binary', 'GenBank_binary']].max(axis=1).sum()

# Step 4: Count how many rows have a 1 in 'GenBank' but not in 'BOLD' and vice versa
count_genbank_only = ((df['GenBank_binary'] == 1) & (df['BOLD_binary'] == 0)).sum()
count_bold_only = ((df['BOLD_binary'] == 1) & (df['GenBank_binary'] == 0)).sum()

# Create a new DataFrame with the results
result_df = pd.DataFrame({
    'Count_all_taxa': len(df),
    'Count_BOLD': [count_bold],
    'Count_GenBank': [count_genbank],
    'Count_Either': [count_either],
    'Count_GenBank_Only': [count_genbank_only],
    'Count_BOLD_Only': [count_bold_only]
}).transpose()

result_df.to_csv("/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/Red Sea species list processing - chris - NCBI + BOLD counts/count_summary.csv", header=False)