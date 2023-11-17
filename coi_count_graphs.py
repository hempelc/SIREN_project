import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os

# FIle with counts per species as .csv
infile = "/Users/christopherhempel/Google Drive/KAUST/SIREN project data/red_sea_species_list_standardized_with_coi_counts.csv"
# Graph output directory
plot_outdir = (
    "/Users/christopherhempel/Google Drive/KAUST/SIREN project data/coi_counts_graphs"
)
# Min number of available COI seq for taxon to count as found
min_seqs = 1
# Maximum number of taxa to show in barplots
num_taxa_barplots = 50
# Determine ranks to graph
ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
# ranks = ["phylum"]
# Kingdom filter (all listed will be kept)
kingdoms = ["Animalia", "Archaea", "Bacteria", "Chromista", "Plantae"]
# kingdoms = ["Animalia"]


# Make plot directory
os.makedirs(plot_outdir, exist_ok=True)

# Import data
df = pd.read_csv(infile)

# Filter the DataFrame based on the kingdom filter
df = df[df["kingdom"].isin(kingdoms)]

# Turn to p/a based on min_seqs
df["Number of COI sequences on GenBank"] = df[
    "Number of COI sequences on GenBank"
].apply(lambda x: 1 if x >= min_seqs else 0)
df = df.rename(columns={"Number of COI sequences on GenBank": "COI sequence available"})

# Strip plots

## Loop over ranks
for rank in ranks:
    ## Group counts of ranks
    rank_series = df.groupby(rank)["COI sequence available"].sum()
    ## Make separate violin plot
    fig = px.strip(rank_series, width=300).update_traces(jitter=1)
    median = rank_series.median()
    if int(median) == median and isinstance(median, float):
        median = round(median)
    mean = round(rank_series.mean(), 1)
    fig.update_xaxes(showticklabels=False)
    fig.update_layout(
        xaxis_title=f"{rank} - median: {median} - mean: {mean}",
        yaxis_title=f"# of taxa with ≥ {min_seqs} COI seqs on GenBank",
    )
    fig.show()
    fig.write_image(
        os.path.join(plot_outdir, f"stripplot_{rank}_{min_seqs}_seqs.png"), width=300
    )

if "pecies" in ranks[-1]:
    species_series = df.groupby(ranks[-1])["COI sequence available"].sum()
    species_df = pd.DataFrame(
        {
            "spacer": ["a", "a"],
            "category": ["without", "with"],
            "count": [(species_series == 0).sum(), (species_series == 1).sum()],
        }
    )
    fig = px.bar(
        species_df,
        text_auto=True,
        x="spacer",
        y="count",
        color="category",
        title=f"Species with and without ≥ {min_seqs} COI seqs",
    )
    fig.update_xaxes(showticklabels=False)
    fig.update_layout(
        xaxis_title="",
        yaxis_title="",
    )
    fig.update_traces(insidetextanchor="middle")
    fig.show()
    fig.write_image(
        os.path.join(plot_outdir, f"stackedbarplot_species_{min_seqs}_seqs.png")
    )

# Bar plots
## Loop over ranks
for rank in ranks:
    ### Group counts of ranks and count how many taxa are grouped
    rank_df = df.groupby(rank).agg({"COI sequence available": ["sum", "count"]})
    rank_df.columns = [
        "COI sequence available",
        "Taxon count",
    ]
    ### Determine number of taxa with missing COI sequences
    rank_df["COI sequence not available"] = (
        rank_df["Taxon count"] - rank_df["COI sequence available"]
    )
    ### Determine relative number of counts
    rank_df["proportion"] = (
        rank_df["COI sequence available"] / rank_df["Taxon count"] * 100
    )

    # Save df
    rank_df.to_csv(
        "/Users/christopherhempel/Google Drive/KAUST/SIREN project data/red_sea_species_list_final_{rank}_{min_seqs}_seqs.csv"
    )

    ### Define barplot title based on number of taxa
    if len(rank_df) >= num_taxa_barplots:
        barplot_title_total = f"Total number of taxa with ≥ {min_seqs} COI seqs available on GenBank - {rank} - first {num_taxa_barplots} taxa"
        barplot_title_relative = f"Relative number of taxa with ≥ {min_seqs} COI seqs available on GenBank - {rank} - first {num_taxa_barplots} taxa"
    else:
        barplot_title_total = f"Total number of taxa with ≥ {min_seqs} COI seq available on GenBank - {rank}"
        barplot_title_relative = f"Relative number of taxa with ≥ {min_seqs} COI seqs available on GenBank - {rank}"

    ### Total barplots
    #### Sort and keep top X taxa, where X = num_taxa_barplots
    rank_df_total = rank_df["COI sequence available"].sort_values(ascending=False)[
        :num_taxa_barplots
    ]
    #### Make bar plots
    barplot = px.bar(rank_df_total)
    barplot.update_layout(
        xaxis_title="Taxa",
        yaxis_title="Counts",
        title={
            "text": barplot_title_total,
            "font": {"size": 15},
        },
        showlegend=False,
        width=max(300, 25 * len(rank_df_total)),
    )
    barplot.update_xaxes(tickangle=45)
    barplot.show()
    barplot.write_image(
        os.path.join(plot_outdir, f"barplot_total_{rank}_{min_seqs}_seqs.png")
    )

    ### Relative barplots
    #### Sort and keep top X taxa, where X = num_taxa_barplots
    rank_df_relative = rank_df["proportion"].sort_values(ascending=False)[
        :num_taxa_barplots
    ]
    #### Make bar plots
    barplot = px.bar(rank_df_relative)
    barplot.update_layout(
        xaxis_title="Taxa",
        yaxis_title="%",
        title={
            "text": barplot_title_relative,
            "font": {"size": 15},
        },
        showlegend=False,
        width=max(300, 25 * len(rank_df_relative)),
    )
    barplot.update_xaxes(tickangle=45)
    barplot.update_layout(yaxis_range=[0, 100])
    barplot.show()
    barplot.write_image(
        os.path.join(plot_outdir, f"barplot_relative_{rank}_{min_seqs}_seqs.png")
    )

    ### Stacked barplot
    stacked_df = rank_df.reset_index().sort_values("Taxon count", ascending=False)
    stacked_bar = px.bar(
        stacked_df,
        x=rank,
        y=["COI sequence available", "COI sequence not available"],
        title=f"Number of taxa on rank {rank} with and without ≥ {min_seqs} COI seqs available on GenBank",
    )

    stacked_bar.update_layout(
        xaxis_title="Taxa",
        yaxis_title="Sequence counts",
        width=max(300, 25 * len(stacked_df)),
        legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99),
    )
    stacked_bar.update_xaxes(tickangle=45)
    stacked_bar.show()
    stacked_bar.write_image(
        os.path.join(plot_outdir, f"stacked_barplot_{rank}_{min_seqs}_seqs.png")
    )

    ### Bubble plot
    bubble_df = rank_df.reset_index()
    bubble_plot = px.scatter(
        rank_df.reset_index(),
        x="COI sequence available",
        y="COI sequence not available",
        size="Taxon count",
        color="proportion",
        hover_name=rank,
        size_max=60,
        color_continuous_scale="Plasma_r",
        width=600,
        height=500,
        title=rank,
    )
    bubble_plot.update_traces(marker=dict(line=dict(width=0)))
    bubble_plot.show()
    bubble_plot.write_image(
        os.path.join(plot_outdir, f"bubbleplot_{rank}_{min_seqs}_seqs.png")
    )
