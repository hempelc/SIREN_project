import pandas as pd

file1 = "MZGdata-coi.csv"
file3 = (
    "red_sea_species_list_standardized_with_coi_counts_depth_lat_lon_search-updated.csv"
)

df1 = pd.read_csv(
    file1, usecols=[1, 16, 17], header=None, names=["species", "lon", "lat"]
)
df3 = pd.read_csv(file3)

numspecimens = df1["species"].isin(df3["species"]).sum()
print(f"Number of MZGdb specimens in our species list: {numspecimens}")

numspecies = df1["species"][df1["species"].isin(df3["species"])].nunique()
print(f"Number of MZGdb species in our species list: {numspecies}")

# Unique species in MZGdb:
species1 = set(df1["species"])
len(species1)

# Unique species in BOLD:
species3 = set(df3["species"])

# Number and percentage of species in our list found in MZGdb
num_notfound = len([x for x in species3 if x not in species1])
len(species3) - num_notfound
(len(species3) - num_notfound) / len(species3)

# Species found in MZGdb without COI seq in BOLD or NCBI
spec_found = [x for x in species3 if x in species1]
filtered_df = df3[df3["species"].isin(spec_found)]
uncovered_species_df = filtered_df[
    (filtered_df["COI_sequence_available_on_GenBank"] == 0)
    & (filtered_df["COI_sequence_available_on_BOLD"] == 0)
]
uncovered_phylum_counts = uncovered_species_df["phylum"].value_counts()

# Add lat lon info from MZGdb to existing table
df1["lat"] = df1["lat"].str.replace(":q!", "", regex=False).astype(float)
with_lat_lon = df1[df1["lon"] != -999]
just_species_from_list = with_lat_lon[with_lat_lon["species"].isin(spec_found)]
## Num species and specimens with lat lon
len(set(just_species_from_list["species"]))
len(just_species_from_list)
## Just within Red Sea coordinates
mask = (
    (just_species_from_list["lat"] >= 12.5)
    & (just_species_from_list["lat"] <= 30)
    & (just_species_from_list["lon"] >= 32)
    & (just_species_from_list["lon"] <= 44)
)
just_red_sea_lat_lon = just_species_from_list[mask]
## Num species with lat lon in Red Sea
len(set(just_red_sea_lat_lon["species"]))

## Add lat lon info
### Iterate over rows in the second DataFrame
for index, row in just_red_sea_lat_lon.iterrows():
    species = row["species"]
    lat = row["lat"]
    lon = row["lon"]

    ### Find the corresponding species in the first DataFrame
    mask = df3["species"] == species

    ### If lat_lon_available is 0, set it to 1
    df3.loc[mask & (df3["lat_lon_available"] == 0), "lat_lon_available"] = 1

    ### For lat_values
    existing_lat_values = (
        df3.loc[mask, "lat_values"]
        .fillna("")
        .astype(str)
        .apply(lambda x: x.split(",") if x else [])
    )
    updated_lat_values = existing_lat_values.apply(
        lambda x: ",".join(sorted(set(x + [str(lat)])))
    )

    ### For lon_values
    existing_lon_values = (
        df3.loc[mask, "lon_values"]
        .fillna("")
        .astype(str)
        .apply(lambda x: x.split(",") if x else [])
    )
    updated_lon_values = existing_lon_values.apply(
        lambda x: ",".join(sorted(set(x + [str(lon)])))
    )

    ### Update the DataFrame with unique values
    df3.loc[mask, "lat_values"] = updated_lat_values
    df3.loc[mask, "lon_values"] = updated_lon_values


### Remove leading commas in lat_values and lon_values if they were empty
df3["lat_values"] = df3["lat_values"].str.lstrip(",")
df3["lon_values"] = df3["lon_values"].str.lstrip(",")

df3.to_csv(
    "red_sea_species_list_standardized_with_coi_counts_depth_lat_lon_search-updated_with_MGZcoords.csv",
    index=False,
)

# Expand lat lon columns
# Create a copy of the relevant columns for lat/lon
df_lat_lon = df3[["lat_values", "lon_values"]].copy()

# Ensure lat_values and lon_values are lists
df_lat_lon["lat_values"] = df_lat_lon["lat_values"].str.split(",")
df_lat_lon["lon_values"] = df_lat_lon["lon_values"].str.split(",")

# Explode both lat_values and lon_values into separate rows
df_lat_lon = df_lat_lon.explode("lat_values").explode("lon_values")

# Drop any rows where lat_values or lon_values are empty or NaN
df_lat_lon = df_lat_lon.dropna()

# Drop duplicates
df_lat_lon = df_lat_lon.drop_duplicates()

# Reset index
df_lat_lon = df_lat_lon.reset_index(drop=True)
df_lat_lon.to_csv(
    "lat_lon_for_density_map.csv",
    index=False,
)
