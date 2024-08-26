import pandas as pd

df1 = pd.read_csv(
    "red_sea_species_list_standardized_synonyms_merged_unaccepted_merged_just_animals_no_aves_no_insecta_new_BOLD_species.csv"
)
df2 = pd.read_csv(
    "Dirk_lat_long_rec_search_noLand_allSpecies_standardized.csv"
)

# Convert count columns to binary
df1["COI_sequence_available_on_GenBank"] = (
    df1["COI_sequences_on_GenBank"] >= 1
).astype(int)
df1["COI_sequence_available_on_BOLD"] = (df1["COI_sequences_on_BOLD"] >= 1).astype(int)

# Drop the original columns
df1.drop(
    columns=["COI_sequences_on_GenBank", "COI_sequences_on_BOLD", "Notes"], inplace=True
)

# Turn all values in the 'values' column to positive
df2["depth"] = df2["depth"].abs()

# Initialize the new columns in df1
df1["lat_lon_available"] = 0
df1["lat_values"] = ""
df1["lon_values"] = ""
df1["depth_available"] = 0
df1["depth_values"] = ""

# Update the columns based on the presence of species in df2
for i, species in enumerate(df1["species"]):
    if species in df2["species"].values:
        df1.at[i, "lat_lon_available"] = 1
        lat_values = df2.loc[df2["species"] == species, "lat"].unique()
        long_values = df2.loc[df2["species"] == species, "long"].unique()
        depth_values = df2.loc[df2["species"] == species, "depth"]

        if not depth_values.isna().all():
            depth_values = depth_values.dropna().unique()
            df1.at[i, "depth_values"] = ",".join(map(str, depth_values))
            df1.at[i, "depth_available"] = 1

        df1.at[i, "lat_values"] = ",".join(map(str, lat_values))
        df1.at[i, "lon_values"] = ",".join(map(str, long_values))
    else:
        df1.at[i, "lat_values"] = ""
        df1.at[i, "lon_values"] = ""
        df1.at[i, "depth_values"] = ""

df1.to_csv(
    "red_sea_species_list_standardized_with_coi_counts_depth_lat_lon.csv",
    index=False,
)
