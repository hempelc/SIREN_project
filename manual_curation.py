# Manually curate DB species list
import pandas as pd
import re

df = pd.read_excel(
    "/Users/christopherhempel/Google Drive/KAUST/SIREN project/Red Sea species list.xlsx",
    sheet_name="All species non-stand. + DBs",
).dropna()


# Function to clean up values in Species column
def clean_species(value):
    # Remove text in brackets
    value = re.sub(r"\([^)]*\)", "", value).strip()
    # Remove entries with 'cf.' and 'sp.'
    if "cf." in value or "sp." in value:
        return None
    # Remove entries ending with " sp"
    if value.endswith(" sp"):
        return None
    # Cut subspecies, variance, and forma names down to the species name by
    # removing all but the first 2 words
    if len(value.split()) > 2:
        return " ".join(value.split()[:2])
    # Remove entries that are only one word
    if len(value.split()) == 1:
        return None
    return value


# Apply cleaning function
df["Species"] = df["Species"].apply(clean_species)

# Drop rows with None in 'Column1'
df = df.dropna(subset=["Species"])


# Manually check for species names that dont contain A-Z, a-z, "-" and space
# Function to check if a string contains characters other than A-Za-z space "-"
def contains_non_alpha(value):
    return bool(re.search(r"[^A-Za-z -]", value))


# Identify rows with characters other than A-Za-z
df["Contains_Non_Alpha"] = df["Species"].apply(contains_non_alpha)

# Manually inspect
df[df["Contains_Non_Alpha"]]

# Manually correct a few specific cases
df["Species"] = df["Species"].str.replace("Actäa", "Actaea")
df["Species"] = df["Species"].str.replace("Mülleria", "Muelleria")
df["Species"] = df["Species"].str.replace("okhaënsis", "okhaensis")

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
    "/Users/christopherhempel/Google Drive/KAUST/SIREN project/Red Sea species list manually cleaned.csv",
    index=False,
)
