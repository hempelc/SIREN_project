import pandas as pd
import plotly.express as px

# File with counts per species as .csv
infile = "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/Red Sea species list processing - chris - NCBI + BOLD counts/red_sea_species_list_standardized_synonyms_merged_unaccepted_merged_just_animals_no_aves_no_insecta.csv"

# Import data
df = pd.read_csv(infile)

# Filter the DataFrame based on the kingdom filter
df = df[df["kingdom"]=="Animalia"]

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

# Make figure
# Format df and manually add missing info
result_df = result_df.reset_index()
result_df.columns = ["metric", "count"]
result_df['metric'] = ["Animals*", "BOLD", "GenBank", "GenBank + BOLD", "GenBank exclusive", "BOLD exclusive"]
result_df.loc[len(result_df)] = {"metric": "All", "count": 7882}
result_df.loc[len(result_df)] = {"metric": "Standardized", "count": 6897}
result_df = result_df.reindex([6,7,0,2,1,3,4,5])

# Plot
fig = px.bar(result_df, x="metric", y="count", text_auto=True, color_discrete_sequence=["#b8860b"], labels={'count':'Number of species'})
fig.update_layout(xaxis_title=None)
fig.write_image("/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/Red Sea species list processing - chris - NCBI + BOLD counts/count_summary.png", height=400)
