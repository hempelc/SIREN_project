import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import plotly.express as px

# Read the Excel file into a pandas DataFrame
file_path = "lat lon depth bold all specimens no aves no insects.csv"
df = pd.read_csv(file_path)

# Create a GeoDataFrame from the DataFrame by converting lat long columns to Point geometries
geometry = [Point(xy) for xy in zip(df["long"], df["lat"])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

# Download the world map from GeoPandas datasets
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

# Plot the world map
fig, ax = plt.subplots(figsize=(10, 6))
world.plot(ax=ax, color="lightgrey")

# Plot the points on the world map
gdf.plot(ax=ax, marker="o", color="red", markersize=2)

# Set plot title and show the plot
plt.show()


# Interactive version
df["hover_text"] = (
    "Sample ID: "
    + df["sampleid"].astype(str)
    + "<br>"
    + "Species: "
    + df["name"].astype(str)
)

# Create a scatter mapbox plot using Plotly Express
fig = px.scatter_mapbox(df, lat="lat", lon="long", hover_name="hover_text", zoom=0)

# Customize the map layout
fig.update_layout(mapbox_style="open-street-map", margin=dict(l=0, r=0, t=0, b=0))

# Show the plot
fig.show()


# CODE TO REMOVE LAND POINTS

# Perform a spatial join to identify points that fall on land
points_on_land = gpd.sjoin(gdf, world, how="inner", op="intersects")

# Get the indices of points on land
indices_on_land = points_on_land.index

# Filter out points on land from the original GeoDataFrame
gdf_filtered = gdf[~gdf.index.isin(indices_on_land)]

# Make hover name
gdf_filtered["hover_text"] = (
    "Sample ID: "
    + gdf_filtered["sampleid"].astype(str)
    + "<br>"
    + "Species: "
    + gdf_filtered["name"].astype(str)
)

# Create a scatter mapbox plot using Plotly Express for points not on land
fig = px.scatter_mapbox(
    gdf_filtered, lat="lat", lon="long", hover_name="hover_text", zoom=1
)

# Customize the map layout
fig.update_layout(mapbox_style="open-street-map", margin=dict(l=0, r=0, t=0, b=0))

# Show the plot
fig.show()

# Non-interactive version:

# Create a GeoDataFrame from the DataFrame by converting lat long columns to Point geometries
geometry = [Point(xy) for xy in zip(gdf_filtered["long"], gdf_filtered["lat"])]
gdf = gpd.GeoDataFrame(gdf_filtered, geometry=geometry, crs="EPSG:4326")

# Plot the world map
fig, ax = plt.subplots(figsize=(10, 6))
world.plot(ax=ax, color="lightgrey")

# Plot the points on the world map
gdf.plot(ax=ax, marker="o", color="red", markersize=2)

# Show the plot
plt.show()

# Just Red Sea
# Define the bounding box for the Red Sea region
red_sea_bbox = [32, 44, 12.5, 30]  # [minx, maxx, miny, maxy]

# Plot the points on the world map within the specified bounding box
fig, ax = plt.subplots(figsize=(10, 6))
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world.plot(ax=ax, color="lightgrey", edgecolor="black")
ax.set_xlim(red_sea_bbox[0], red_sea_bbox[1])
ax.set_ylim(red_sea_bbox[2], red_sea_bbox[3])

gdf.plot(ax=ax, marker="o", color="red", markersize=6)

# Show the plot
plt.show()

# Get info on number of specimens and species in Red Sea
mask = (
    (gdf["lat"] >= 12.5)
    & (gdf["lat"] <= 30)
    & (gdf["long"] >= 32)
    & (gdf["long"] <= 44)
)
filtered_df = gdf[mask][["sampleid", "lat", "long", "depth", "name"]]
print("Specimens in Red Sea: " + str(len(filtered_df)))
print("Species in Red Sea: " + str(len(filtered_df["name"].drop_duplicates())))
print(
    "Specimens in Red Sea with depth: "
    + str(len(filtered_df[filtered_df["depth"].notna()]))
)
print(
    "Species in Red Sea with depth: "
    + str(len(filtered_df[filtered_df["depth"].notna()]["name"].drop_duplicates()))
)

# Plot all complete information
non_animal_phyla = [
    "Rhodophyta",
    "Amoebozoa",
    "Pyrrophycophyta",
    "Chlorophyta",
    "Magnoliophyta",
    "Thaumarchaeota",
]
gdf_filtered_animals_noNaN = gdf_filtered[
    ~gdf_filtered["phylum"].isin(non_animal_phyla)
].dropna(subset=["species"])
gdf_filtered_animals_with_species = gdf_filtered_animals_noNaN[
    ~gdf_filtered_animals_noNaN["species"].str.contains("sp\.|cf\.")
]
gdf_filtered_animals_with_species.to_csv(
    "Dirk_lat_long_rec_search_noLand_allSpecies.csv"
)


# Create a GeoDataFrame from the DataFrame by converting lat long columns to Point geometries
geometry = [
    Point(xy)
    for xy in zip(
        gdf_filtered_animals_with_species["long"],
        gdf_filtered_animals_with_species["lat"],
    )
]
gdf = gpd.GeoDataFrame(
    gdf_filtered_animals_with_species, geometry=geometry, crs="EPSG:4326"
)

# Just Red Sea
# Define the bounding box for the Red Sea region
red_sea_bbox = [32, 44, 12.5, 30]  # [minx, maxx, miny, maxy]

# Plot the points on the world map within the specified bounding box
fig, ax = plt.subplots(figsize=(10, 6))
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world.plot(ax=ax, color="lightgrey", edgecolor="black")
ax.set_xlim(red_sea_bbox[0], red_sea_bbox[1])
ax.set_ylim(red_sea_bbox[2], red_sea_bbox[3])

gdf.plot(ax=ax, marker="o", color="red", markersize=6)

# Get info on species not in our list
species_not_in_our_data = gdf_filtered_animals_with_species[
    ~gdf_filtered_animals_with_species["species"].isin(df["name"])
]
species_not_in_our_data.to_csv(
    "Dirk_lat_long_rec_search_noLand_notInOurList.csv"
)
