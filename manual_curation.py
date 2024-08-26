# Manually curate DB species list
import pandas as pd
import re

# Import data and drop empty rows
df = pd.read_excel(
    "Red Sea species list database mining.xlsx",
    sheet_name="Non-std. sp. + DBs not cleaned",
).dropna()


# Function to clean up values in Species column
def clean_species(value):
    # Remove text in brackets
    value = re.sub(r"\([^)]*\)", "", value).strip()
    # Remove entries with 'cf.' and 'sp.'
    if "cf." in value or "sp." in value:
        return None
    # Remove entries ending with " sp" and " cf"
    if value.endswith(" sp") or value.endswith(" cf"):
        return None
    # Cut subspecies, variance, and forma names down to the species name by
    # removing all but the first 2 words
    if len(value.split()) > 2:
        return " ".join(value.split()[:2])
    # Remove entries that are only one word
    if len(value.split()) == 1:
        return None
    return value


# Function to check if a string contains characters other than A-Za-z space "-"
def contains_non_alpha(value):
    return bool(re.search(r"[^A-Za-z -]", value))


# Apply cleaning function
df["Species"] = df["Species"].apply(clean_species)

# Drop rows with None in 'Species'
df = df.dropna(subset=["Species"])

# Manually check for species names that dont contain A-Z, a-z, "-" and space
## Identify rows with characters other than A-Za-z
df["Contains_Non_Alpha"] = df["Species"].apply(contains_non_alpha)

## Manually inspect
df[df["Contains_Non_Alpha"]]

## Manually correct a few specific cases
df["Species"] = df["Species"].str.replace("Actäa", "Actaea")
df["Species"] = df["Species"].str.replace("Mülleria", "Muelleria")
df["Species"] = df["Species"].str.replace("okhaënsis", "okhaensis")

# Drop duplicate entries that are a result of cutting off subspecies names etc
df = df.drop_duplicates()

# Aggregate species and databases
grouped = df.groupby("Species")["Database"].agg(list).reset_index()
grouped["Database_all"] = grouped["Database"].apply(lambda x: ", ".join(x))
result_df = df.drop_duplicates(subset="Species").merge(
    grouped, on="Species", how="left"
)

# Clean and save df
result_df = result_df.drop(["Database_x", "Contains_Non_Alpha", "Database_y"], axis=1)
result_df = result_df.rename(columns={"Database_all": "Database"})

result_df.to_csv(
    "red_sea_species_list_manually_cleaned.csv",
    index=False,
)
