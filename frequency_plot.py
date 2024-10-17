import pandas as pd
import plotly.express as px

# Read CSV file
df = pd.read_csv(
    "red_sea_species_list_standardized_synonyms_merged_unaccepted_merged_just_animals_no_aves_no_insecta_new_BOLD_species_search-updated.csv"
)

# Define custom breaks and labels
custom_breaks = [-0.5, 0.5, 5, 20, 50, 100, 500, 1105]
custom_labels = ["0", "1-5", "6-20", "21-50", "51-100", "101-500", ">500"]

# Create intervals for NCBI using pd.cut
df["intervals"] = pd.cut(
    df["COI_sequences_on_GenBank"],
    bins=custom_breaks,
    labels=custom_labels,
    include_lowest=True,
)

# Create a bar plot
fig = px.bar(
    df["intervals"].value_counts().sort_index(),
    labels={"value": "Number of Species", "index": "Number of sequences per species"},
    title="Frequency of COI sequences per species",
    color_discrete_sequence=["darkgoldenrod"],
    text_auto=True,
)

# Update plot layout
fig.update_layout(
    xaxis_title="Number of sequences per species",
    yaxis_title="Number of species",
    showlegend=False,
)

# Show the plot
fig.show()

# Save the plot as a PNG file
fig.write_image(
    "frequency_plot_NCBI.svg",
    height=400,
)

# Save the plot as a PNG file
fig.write_image(
    "frequency_plot_NCBI.png",
    height=400,
)
