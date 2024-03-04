import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

# Read the Excel file into a pandas DataFrame
file_path = "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/BOLD processing/lat lon bold.xlsx"
df = pd.read_excel(file_path)

# Create a GeoDataFrame from the DataFrame by converting lat long columns to Point geometries
geometry = [Point(xy) for xy in zip(df['long'], df['lat'])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

# Download the world map from GeoPandas datasets
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# Plot the world map
fig, ax = plt.subplots(figsize=(10, 6))
world.plot(ax=ax, color='lightgrey')

# Plot the points on the world map
gdf.plot(ax=ax, marker='o', color='red', markersize=2)

# Set plot title and show the plot
plt.show()


import pandas as pd
import plotly.express as px

# Read the Excel file into a pandas DataFrame
file_path = "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/BOLD processing/lat lon bold.xlsx"
df = pd.read_excel(file_path)

df['hover_text'] = "Sample ID: " + df['sampleid'].astype(str) + '<br>' + "Species: " + df['name'].astype(str)

# Create a scatter mapbox plot using Plotly Express
fig = px.scatter_mapbox(df,
                        lat='lat',
                        lon='long',
                        hover_name='hover_text',
                        zoom=0)

# Customize the map layout
fig.update_layout(mapbox_style="open-street-map",
                  margin=dict(l=0, r=0, t=0, b=0))

# Show the plot
fig.show()


# CODE TO REMOVE LAND POINTS

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import plotly.express as px

# Read the Excel file into a pandas DataFrame
file_path = "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/BOLD processing/lat lon bold.xlsx"
df = pd.read_excel(file_path)

# Create a GeoDataFrame from the DataFrame by converting lat long columns to Point geometries
geometry = [Point(xy) for xy in zip(df['long'], df['lat'])]
gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')

# Download the world map from GeoPandas datasets
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# Perform a spatial join to identify points that fall on land
points_on_land = gpd.sjoin(gdf, world, how='inner', op='intersects')

# Get the indices of points on land
indices_on_land = points_on_land.index

# Filter out points on land from the original GeoDataFrame
gdf_filtered = gdf[~gdf.index.isin(indices_on_land)]

# Make hover name
gdf_filtered['hover_text'] = "Sample ID: " + gdf_filtered['sampleid'].astype(str) + '<br>' + "Species: " + gdf_filtered['name'].astype(str)

# Create a scatter mapbox plot using Plotly Express for points not on land
fig = px.scatter_mapbox(gdf_filtered,
                        lat='lat',
                        lon='long',
                        hover_name='hover_text',
                        zoom=1)

# Customize the map layout
fig.update_layout(mapbox_style="open-street-map",
                  margin=dict(l=0, r=0, t=0, b=0))

# Show the plot
fig.show()

# Non-interactive version:

# Create a GeoDataFrame from the DataFrame by converting lat long columns to Point geometries
geometry = [Point(xy) for xy in zip(gdf_filtered['long'], gdf_filtered['lat'])]
gdf = gpd.GeoDataFrame(gdf_filtered, geometry=geometry, crs='EPSG:4326')

# Plot the world map
fig, ax = plt.subplots(figsize=(10, 6))
world.plot(ax=ax, color='lightgrey')

# Plot the points on the world map
gdf.plot(ax=ax, marker='o', color='red', markersize=2)

# Set plot title and show the plot
plt.show()
